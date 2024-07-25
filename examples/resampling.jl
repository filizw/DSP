include("../src/DSP.jl")

using CairoMakie

fs = 10
N = 30
t = range(0, step = 1/fs, length = N)
x1 = cos.(2*π.*t)

x2 = DSP.resample(x1, 2, 3)

fig = Figure()
ax1 = Axis(fig[1,1])
ax2 = Axis(fig[2,1])
stem!(ax1, x1)
stem!(ax2, real.(x2))
current_figure()