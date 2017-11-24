# paqueteFourier

Paquete para [Julia](https://julialang.org/) que aproxima, manipula y resuelve EDO con funciones representadas mediante series de Fourier finitas (polinomios trigonométricos).

## Instalación

Primero es necesario haber [instalado Julia](https://julialang.org/downloads/). Luego, desde la ventana de comandos de Julia, se debe ejecutar:

```julia-repl
julia> Pkg.clone("https://github.com/lfiguero/paqueteFourier.git")
```

## Ejemplo de utilización

```julia
using paqueteFourier
f = trigPoly([0.0, 0.0, 0.0, 0.5im, -0.5im]) # ¡Comas importan!
f(π/8)
vals(f)
nd = nodes(length(f.coefs))
f(nd)
g = trigPoly(t -> -1/2*cos(2t), 7)
innerProduct(f,g)
h = diff(g)
norm(f-h)
f*g
3*f - 2.0im*g
wave = t -> t^2*(2π-t)^2*exp((-1.3+3.9im)t)
w = trigPoly(wave, 18)
tpfplot(wave, w, labels=["wave", "I_18(wave)"])
```

## Contenidos

Este paquete provee el tipo `trigPoly` que representa un polinomio trigonométrico mediante sus coeficientes respecto a la base 1, exp(-im θ), exp(im θ), exp(-2im θ), exp(2im θ), etc. Para crear y manipular `trigPoly` este paquete provee las funciones descritas a continuación.

### trigPoly (constructores)

`trigPoly(x)`, donde `x` es un vector de números, devuelve un polinomio trigonométrico que es combinación lineal de 1, exp(-im θ), exp(im θ), exp(-2im θ), exp(2im θ), … con los coeficientes indicados en `x`.

`trigPoly(f, n)`, donde `f` es una función a valores complejos y `n` un entero no negativo, devuelve el polinomio trigonométrico que interpola a `f` en la colección de nodos de Fourier de `n` elementos.


## Cómo consultar ayuda del paquete

## Cómo examinar el código del paquete
