module paqueteFourier

import Base: +, -, *, /

import PyPlot

# Notar que los métodos de +, -, *, / y diff se exportan solos
export trigPoly, innerProduct, nodes, vals, tpfplot

# El objeto básico de este paquete es un trigPoly, cuyo campo coefs contiene
# los coeficientes de la representación de un polinomio trigonométrico como
# combinación lineal de 1, exp(-im θ), exp(im θ), exp(-2im θ), exp(2im θ), …,
# en ese orden.
# Junto a la definición de este 'type' Julia automáticamente deja disponible un
# constructor, que es una función que recibe un vector de números complejos y
# devuelve un trigPoly con ese vector en el campo coefs.
# A continuación justificaré mi elección de ordenamiento para los coeficientes
# con la ayuda de una tabla.  Denoto al orden elegido (0, -1, 1, -2, 2, …) por
# A, al orden visto en clase (-N/2, -N/2+1, …, -1, 0, 1, …, N/2-2, N/2) por B y
# al orden por defecto de la FFT (0, 1, …, N/2-1, -N/2, -N/2+1, …, -1) por C (B
# y C son ligeramente distintos si N es impar, pero la tabla sigue valiendo).
#                                 ┌───┬───┬───┐
#                                 │ A │ B │ C │
# ┌───────────────────────────────┼───┼───┼───┤
# │                          Suma │ ✓ │ ✗ │ ✗ │
# │Interpretación de coeficientes │ ✓ │ ✗ │ ~ │
# │     Vector de números de onda │ ~ │ ✓ │ ✗ │
# │           Interpolación (FFT) │ ✗ │ ~ │ ✓ │
# ├───────────────────────────────┴───┴───┴───┤
# │   ✓ = Fácil, ~ = Moderado, ✗ = Difícil    │
# └───────────────────────────────────────────┘
"""
    trigPoly(x)

Polinomio trigonométrico que es combinación lineal de 1, exp(-im θ), exp(im θ), exp(-2im θ), … con los coeficientes indicados en `x`.
"""
type trigPoly
	coefs::Vector{Complex128}
end

"""
    innerProduct(f,g)

Producto interior L² entre los trigPoly `f` y `g`.
"""
function innerProduct(f::trigPoly, g::trigPoly)
	ac = 0.0 + 0.0im
	for k = 1:min(length(f.coefs),length(g.coefs))
		ac += f.coefs[k] * conj(g.coefs[k])
	end
	2π*ac
end

# Sin cargar ningún paquete Julia ya tiene una función llamada norm; la nuestra
# será un método más
"""
    norm(f::trigPoly)

Norma L² de `f`.
"""
Base.norm(f::trigPoly) = sqrt(real(innerProduct(f,f)))

"""
    nodes(n)

Arreglo de `n` nodos de Fourier.
"""
function nodes(n::Integer)
	2π*(0:n-1)/n
end

# Función que transforma un vector del orden C al orden A; Julia reemplazará T
# con el tipo de número que contenga el vector con el cual se alimente a esta
# función
function reordering{T}(v::Vector{T})
	n = length(v)
	out = Vector{T}(n) # Reserva de memoria
	floorhalfn = n÷2
	for k = 1:floorhalfn
		out[2*k] = v[n+1-k]
	end
	for k = 1:n-floorhalfn
		out[2*k-1] = v[k]
	end
	out
end

# Función que transforma un vector del orden A al orden C
function reverse_reordering{T}(v::Vector{T})
	n = length(v)
	floorhalfn = n÷2
	out = Vector{T}(n)
	for k = 1:floorhalfn
		out[n+1-k] = v[2*k]
	end
	for k = 1:n-floorhalfn
		out[k] = v[2*k-1]
	end
	out
end

# Interpolador que asume que los valores recibidos corresponden a la evaluación
# de una función en los nodos de Fourier
function interp{T}(v::Vector{T})
	modes = fft(v)/length(v)
	trigPoly(reordering(modes))
end

"""
    trigPoly(f, n)

Polinomio trigonométrico que interpola a la función `f` en `n` nodos de Fourier.
"""
function trigPoly(f::Function, n::Integer)
	vals = map(f, nodes(n)) # Vector con evaluaciones de f
	interp(vals)
end

"""
    vals(f)

Evaluación del polinomio trigonométrico `f` en tantos nodos de Fourier como coeficientes almacena `f`.
"""
function vals(f::trigPoly)
	ifft(reverse_reordering(f.coefs*length(f.coefs)))
end

"""
    vals(f, x)

Evaluación del polinomio trigonométrico `f` en el arreglo o escalar real `x`.
"""
function vals(f::trigPoly, x)
	coefsplus = f.coefs[3:2:end]
	coefsminus = f.coefs[2:2:end]
	# Podría escribir exp(im*x), pero si x es escalar Julia emite una
	# advertencia (no un error) que invita a usar la versión vectorizada de
	# exp; usando map x puede ser escalar o vectorial indistintamente
	expplusix = map(exp, im*x)
	expminusix = conj(expplusix)
	# Arreglo con ceros o escalar cero según lo que sea x
	outplus = zero(expplusix)
	outminus = zero(expplusix)
	for val in reverse(coefsplus)
		outplus = expplusix.*(outplus + val)
	end
	for val in reverse(coefsminus)
		outminus = expminusix.*(outminus + val)
	end
	outplus + outminus + f.coefs[1]
end

# Evaluador por defecto de un trigPoly; permite evaluar un trigPoly f en x con
# la sintaxis clásica f(x)
(f::trigPoly)(x) = vals(f, x)

# Suma
function +(f::trigPoly, g::trigPoly)
	fl = length(f.coefs)
	gl = length(g.coefs)
	retl = max(fl, gl)
	retcoefs = zeros(Complex128, retl)
	retcoefs[1:fl] = f.coefs
	retcoefs[1:gl] += g.coefs
	trigPoly(retcoefs)
end

# Cambio de signo
-(f::trigPoly) = trigPoly(-f.coefs)

# Resta
-(f::trigPoly, g::trigPoly) = f + (-g)

# Suma, resta, multiplicación y división por escalar
function +(f::trigPoly, a::Number)
	retcoefs = f.coefs
	retcoefs[1] += a
	trigPoly(retcoefs)
end
-(f::trigPoly, a::Number) = f + (-a)
*(f::trigPoly, a::Number) = trigPoly(f.coefs*a)
/(f::trigPoly, a::Number) = trigPoly(f.coefs/a)
*(a::Number, f::trigPoly) = f*a
+(a::Number, f::trigPoly) = f+a
-(a::Number, f::trigPoly) = (-f) + a

# Multiplicación
function *(f::trigPoly, g::trigPoly)
	ndof = length(f.coefs) + length(g.coefs)
	paddedf = trigPoly([f.coefs; zeros(Complex128,ndof-length(f.coefs))])
	paddedg = trigPoly([g.coefs; zeros(Complex128,ndof-length(g.coefs))])
	productvals = vals(paddedf) .* vals(paddedg)
	interp(productvals)
end

"""
    diff(f::trigPoly)

Derivada de `f`.
"""
function Base.diff(f::trigPoly)
	retcoefs = Vector{Complex128}(length(f.coefs))
	for k = 1:length(f.coefs)
		#             ↓     modo     ↓
		retcoefs[k] = (-1)^(k+1)*(k÷2) * im * f.coefs[k]
	end
	trigPoly(retcoefs)
end

"""
    diff(f::trigPoly, k::Int64)

`k`-ésima derivada de `f`.
"""
function Base.diff(f::trigPoly, k::Int64)
	@assert k>=0
	if k == 0
		return f
	elseif k == 1
		return diff(f)
	else
		return diff(f, k-1)
	end
end

"""
    tpfplot(f)
    tpfplot(f1, f2, ...)
    tpfplot(f1, f2, ..., labels=[label1, label2, ...])

Grafica a un polinomio trigonométrico (trigPoly) o función (Function) `f` o a la colección de polinomios trigonométricos o funciones `f1`, `f2`, ...; el argumento opcional `labels`, de aparecer, debe especificar nombres para las curvas ingresadas como un arreglo de cadenas de tipo `String`.

La figura obtenida se puede guardar mediante la sintaxis

    fh = tpfplot(...)
    fh[:savefig]("some_filename.eps")
"""
function tpfplot(vf::Vararg{Union{Function,trigPoly}}; labels=[])
	fig, ax = PyPlot.subplots(1, 2, sharey=true)
	PyPlot.sca(ax[1])
	PyPlot.title("Real part")
	PyPlot.sca(ax[2])
	PyPlot.title("Imaginary part")
	handlesreal = []
	handlesimag = []
	defaultNamesFlag = isempty(labels)
	for i = 1:length(vf)
		c = ["b" "g" "r" "c" "m" "y" "k"][(i-1)%7+1] # Especificación de color
		if defaultNamesFlag
			push!(labels, "Curve $(i)")
		end
		if typeof(vf[i]) <: trigPoly
			nsamples = max(150, 6*length(vf[i].coefs))
			x = nodes(nsamples)
			vfix = vals(trigPoly([vf[i].coefs;zeros(Complex128,nsamples-length(vf[i].coefs))]))
			# Agrego el punto de abscisa 2π
			x = [x;2π]
			vfix = [vfix;vfix[1]]
			# Ahora los puntos con abscisas en la colección de
			# nodos de Fourier que corresponde a vf[i]
			xn = nodes(length(vf[i].coefs))
			vfin = vals(vf[i])
			PyPlot.sca(ax[1])
			h = PyPlot.plot(x, real(vfix), c*"-", xn, real(vfin), c*".")
			push!(handlesreal, h)
			PyPlot.sca(ax[2])
			h = PyPlot.plot(x, imag(vfix), c*"-", xn, imag(vfin), c*".")
			push!(handlesimag, h)
		else
			x = linspace(0,2π,2000)
			vfix = map(vf[i], x)
			PyPlot.sca(ax[1])
			h = PyPlot.plot(x, real(vfix), c*"-")
			push!(handlesreal, h)
			PyPlot.sca(ax[2])
			h = PyPlot.plot(x, imag(vfix), c*"-")
			push!(handlesimag, h)
		end
	end
	PyPlot.sca(ax[1])
	PyPlot.legend(map(h -> tuple(h...), handlesreal), tuple(labels...))
	PyPlot.sca(ax[2])
	PyPlot.legend(map(h -> tuple(h...), handlesimag), tuple(labels...))
	fig
end

end # module
