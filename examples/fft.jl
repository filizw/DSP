include("../src/DSP.jl")

using CairoMakie

## fft radix2, ifft radix2
fs = 1000
N = 2^5
t = range(0, step = 1/fs, length = N)
x = sin.(2*π*100.0.*t) .+ 0.5*rand(N)

X = DSP.fft_radix2_dit_r(x)

x1 = DSP.ifft_radix2_dif_r(X)

err = maximum(abs.(x1 - x))
println("err = $err")

fig = Figure()
ax1 = Axis(fig[1,1])
ax2 = Axis(fig[2,1])
stem!(ax1, x)
stem!(ax2, abs.(X))
current_figure()
## fft, ifft
fs = 1000
N = 33
t = range(0, step = 1/fs, length = N)
x = sin.(2*π*100.0.*t) .+ 0.5*rand(N)

X = DSP.fft(x)

x1 = DSP.ifft(X)

err = maximum(abs.(x1 - x))
println("err = $err")

fig = Figure()
ax1 = Axis(fig[1,1])
ax2 = Axis(fig[2,1])
stem!(ax1, x)
stem!(ax2, abs.(X))
current_figure()
##