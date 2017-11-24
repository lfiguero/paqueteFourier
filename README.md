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

Si `f` es un `trigPoly`, el vector de coeficientes puede obtenerse mediante `f.coefs`.

### +, -, \*, /

Los `trigPoly` se pueden cambiar de signo y sumar, restar, multiplicar y dividir por escalar. Dos `trigPoly` se pueden sumar, restar y multiplicar entre sí.

### nodes, vals y evaluación por defecto

`nodes(n)`, donde `n` es un entero positivo, devuelve la colección de los nodos de Fourier en [0,2π) con `n` elementos.

`vals(f)`, donde `f` es un `trigPoly`, devuelve el arreglo con los valores de `f` en `nodes(length(f.coefs))`; esta operación cuesta O(length(f.coefs) log(length(f.coefs))) operaciones.

`vals(f, x)`, donde `f` es un `trigPoly` y `x` un escalar o un vector, devuelve a la evaluación de f en `x`; esta operación cuesta O(length(f.coefs) length(x)).

Si `f` es un `trigPoly` y `x` un escalar o un vector, la sintaxis de evaluación por defecto, `f(x)` equivale a `vals(f, x)`.

### diff

Si `f` es un `trigPoly`, `diff(f)` es un `trigPoly` con la derivada del polinomio trigonométrico representado en `f`.

### innerProduct y norm

Si `f` y `g` son `trigPoly`, `innerProduct(f,g)` devuelve el producto interior L² de los polinomios trigonométricos representados.

Si `f` es un `trigPoly`, `norm(f)` devuelve la norma L² del polinomio trigonométrico representado.

### tpfplot

`tpfplot(f)` grafica a un `trigPoly` o función `f`.

`tpfplot(f1, f2, ...)` grafica a una colección de `trigPoly` o funciones.

`tpfplot(f1, f2, ..., labels=[label1, label2, ...])` especifica etiquetas para las curvas.


## Cómo consultar ayuda

La ayuda de Julia se encuentra en [https://docs.julialang.org/](https://docs.julialang.org/); ahí se puede especificar una versión particular de Julia de ser necesario.

Si uno tipea `?` en la ventana de comandos de Julia, el *prompt* cambia de `julia>` a `help?>`. Entonces se puede tipear un comando y para obtener ayuda acerca de él, de haberla.


## Cómo examinar el código del paquete

El código fuente del paquete está en `src/paqueteFourier.jl` dentro de [la página de GitHub para este repositorio](https://github.com/lfiguero/paqueteFourier).

Notar que al instalar el paquete en una configuración local de Julia, automáticamente se *clona* el paquete (en Linux, por defecto, en `~/.julia/v0.6/paqueteFourier`; reemplazar el `0.6` con la versión de Julia instalada), por lo que también se podrían examinar los contenidos en forma local.
