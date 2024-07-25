include("../src/DSP.jl")

using CairoMakie

fs = 60
N = 20
n = range(0, N-1)

t = range(0, step = 1/fs, length = N)
x = cos.(2*π*10.0.*t)

N1 = 1000
t1 = range(0, t[end], N1)

y = zeros(N1)
for i in range(1, N1)
    y[i] = DSP.interpolate(t1[i], m = t, s = x)
end

fig = Figure()
ax1 = Axis(fig[1,1])
ax2 = Axis(fig[2,1])
stem!(ax1, n, x)
lines!(ax2, t1, y)
current_figure()