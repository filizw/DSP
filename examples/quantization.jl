include("../src/DSP.jl")

using CairoMakie

t = range(0, 2, 1000)
x = cos.(2*π.*t)

## quantization
L = range(-1, 1, 10)
y = DSP.quantize.(x, L=L)

fig = Figure()
ax1 = Axis(fig[1,1])
ax2 = Axis(fig[2,1])
lines!(ax1, t, x)
lines!(ax2, t, y)
current_figure()
## SQNR
N = 4
SQNR = DSP.SQNR(N)
println("SQNR = $SQNR")
## SNR
t = range(0, 2, 1000)
noise = rand(length(t))
x = cos.(2*π.*t)
y = x .+ noise

SNR = DSP.SNR(DSP.power(x), DSP.power(noise))
println("SNR = $SNR")

lines(t, y)
##