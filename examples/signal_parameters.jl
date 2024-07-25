include("../src/DSP.jl")

using CairoMakie

t = range(-3, 3, 1000)
x = DSP.triangular_wave.(t)

## parameters
lines(t, x)

mean = DSP.mean(x)
peak2peak = DSP.peak2peak(x)
energy = DSP.energy(x)
power = DSP.power(x)
rms = DSP.rms(x)

println("mean = $mean")
println("peak2peak = $peak2peak")
println("energy = $energy")
println("power = $power")
println("rms = $rms")
## running mean
running_mean = DSP.running_mean(x, 10)

fig = Figure()
ax1 = Axis(fig[1,1])
ax2 = Axis(fig[2,1])
lines!(ax1, t, x)
lines!(ax2, t, running_mean)
current_figure()
## running energy
running_energy = DSP.running_energy(x, 10)

fig = Figure()
ax1 = Axis(fig[1,1])
ax2 = Axis(fig[2,1])
lines!(ax1, t, x)
lines!(ax2, t, running_energy)
current_figure()
## running power
running_power = DSP.running_power(x, 10)

fig = Figure()
ax1 = Axis(fig[1,1])
ax2 = Axis(fig[2,1])
lines!(ax1, t, x)
lines!(ax2, t, running_power)
current_figure()
##