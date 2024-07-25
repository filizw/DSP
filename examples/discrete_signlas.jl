include("../src/DSP.jl")

using CairoMakie

n = range(-10, 10)

## kronecker
x = DSP.kronecker.(n)
stem(n, x)
## heaviside
x = DSP.heaviside.(n)
stem(n, x)
##