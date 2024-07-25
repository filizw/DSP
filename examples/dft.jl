include("../src/DSP.jl")

using CairoMakie

## dft, idft
fs = 1000
N = 20
t = range(0, step = 1/fs, length = N)
x = sin.(2*π*100.0.*t) .+ 0.5*rand(N)

X = DSP.dft(x)

x1 = DSP.idft(X)

err = maximum(abs.(x1 - x))
println("err = $err")

fig = Figure()
ax1 = Axis(fig[1,1])
ax2 = Axis(fig[2,1])
stem!(ax1, x)
stem!(ax2, abs.(X))
current_figure()
## rdft, irdft
fs = 1000
N = 20
t = range(0, step = 1/fs, length = N)
x = sin.(2*π*100.0.*t) .+ 0.5*rand(N)

X = DSP.rdft(x)

x1 = DSP.irdft(X, N)

err = maximum(abs.(x1 - x))
println("err = $err")

fig = Figure()
ax1 = Axis(fig[1,1])
ax2 = Axis(fig[2,1])
stem!(ax1, x)
stem!(ax2, abs.(X))
current_figure()
##