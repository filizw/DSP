include("../src/DSP.jl")

using CairoMakie

## conv, fast conv
x1 = rand(10)
x2 = [0, 0, 1, 0, 0]

y1 = DSP.conv(x1, x2)
y2 = DSP.conv(x1, x2)

err = maximum(abs.(y1 - y2))
println("err = $err")

fig = Figure()
ax1 = Axis(fig[1,1])
ax2 = Axis(fig[2,1])
ax3 = Axis(fig[3,1])
stem!(ax1, x1)
stem!(ax2, x2)
stem!(ax3, y1)
current_figure()
## overlap add, overlap save
x1 = rand(10)
x2 = [0, 0, 1, 0, 0]

L = 4
y0 = DSP.conv(x1, x2)
y1 = DSP.overlap_add(x1, x2, L)
y2 = DSP.overlap_save(x1, x2, L)

err1 = maximum(abs.(y0 - y1))
err2 = maximum(abs.(y0 - y2))
println("err1 = $err1, err2 = $err2")
## lti filter, filtfilt
b = DSP.firwin_lp_I(10, 0.3)
a = []
x = rand(100)

y1 = DSP.lti_filter(b, a, x)
y2 = DSP.filtfilt(b, a, x)

fig = Figure()
ax1 = Axis(fig[1,1])
ax2 = Axis(fig[2,1])
ax3 = Axis(fig[3,1])
lines!(ax1, x)
lines!(ax2, y1)
lines!(ax3, y2)
current_figure()
## lti amp, lti phase
b = DSP.firwin_lp_I(10, 0.3)
a = []

f = range(0, 1, 1000)

amp = zeros(1000)
phase = zeros(1000)

for i in range(1, 1000)
    amp[i] = DSP.lti_amp(f[i], b, a)
    phase[i] = DSP.lti_phase(f[i], b, a)
end

fig = Figure()
ax1 = Axis(fig[1,1])
ax2 = Axis(fig[2,1])
lines!(ax1, f, amp)
lines!(ax2, f, phase)
current_figure()
##