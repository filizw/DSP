include("../src/DSP.jl")

using CairoMakie

N = 100
n = range(0, N-1)

## rect
x = DSP.rect(N)
lines(n, x)
## triang
x = DSP.triang(N)
lines(n, x)
## hanning
x = DSP.hanning(N)
lines(n, x)
## hamming
x = DSP.hamming(N)
lines(n, x)
## blackman
x = DSP.blackman(N)
lines(n, x)
##