include("../src/DSP.jl")

using CairoMakie

## lp I
b = DSP.firwin_lp_I(10, 0.3)

N = 1000
f = range(0, 1, N)
amp = zeros(N)
phase = zeros(N)

for i in range(1, N)
    amp[i] = DSP.lti_amp(f[i], b, [])
    phase[i] = DSP.lti_phase(f[i], b, [])
end

fig = Figure()
ax1 = Axis(fig[1,1])
ax2 = Axis(fig[2,1])
lines!(ax1, f, amp)
lines!(ax2, f, phase)
current_figure()
## hp I
b = DSP.firwin_hp_I(10, 0.3)

N = 1000
f = range(0, 1, N)
amp = zeros(N)
phase = zeros(N)

for i in range(1, N)
    amp[i] = DSP.lti_amp(f[i], b, [])
    phase[i] = DSP.lti_phase(f[i], b, [])
end

fig = Figure()
ax1 = Axis(fig[1,1])
ax2 = Axis(fig[2,1])
lines!(ax1, f, amp)
lines!(ax2, f, phase)
current_figure()
## bp I
b = DSP.firwin_bp_I(10, 0.2, 0.3)

N = 1000
f = range(0, 1, N)
amp = zeros(N)
phase = zeros(N)

for i in range(1, N)
    amp[i] = DSP.lti_amp(f[i], b, [])
    phase[i] = DSP.lti_phase(f[i], b, [])
end

fig = Figure()
ax1 = Axis(fig[1,1])
ax2 = Axis(fig[2,1])
lines!(ax1, f, amp)
lines!(ax2, f, phase)
current_figure()
## bs I
b = DSP.firwin_bs_I(10, 0.2, 0.3)

N = 1000
f = range(0, 1, N)
amp = zeros(N)
phase = zeros(N)

for i in range(1, N)
    amp[i] = DSP.lti_amp(f[i], b, [])
    phase[i] = DSP.lti_phase(f[i], b, [])
end

fig = Figure()
ax1 = Axis(fig[1,1])
ax2 = Axis(fig[2,1])
lines!(ax1, f, amp)
lines!(ax2, f, phase)
current_figure()
## lp II
b = DSP.firwin_lp_II(11, 0.3)

N = 1000
f = range(0, 1, N)
amp = zeros(N)
phase = zeros(N)

for i in range(1, N)
    amp[i] = DSP.lti_amp(f[i], b, [])
    phase[i] = DSP.lti_phase(f[i], b, [])
end

fig = Figure()
ax1 = Axis(fig[1,1])
ax2 = Axis(fig[2,1])
lines!(ax1, f, amp)
lines!(ax2, f, phase)
current_figure()
## bp II
b = DSP.firwin_bp_II(10, 0.2, 0.3)

N = 1000
f = range(0, 1, N)
amp = zeros(N)
phase = zeros(N)

for i in range(1, N)
    amp[i] = DSP.lti_amp(f[i], b, [])
    phase[i] = DSP.lti_phase(f[i], b, [])
end

fig = Figure()
ax1 = Axis(fig[1,1])
ax2 = Axis(fig[2,1])
lines!(ax1, f, amp)
lines!(ax2, f, phase)
current_figure()
## diff
b = DSP.firwin_diff(10)

N = 1000
f = range(0, 1, N)
amp = zeros(N)
phase = zeros(N)

for i in range(1, N)
    amp[i] = DSP.lti_amp(f[i], b, [])
    phase[i] = DSP.lti_phase(f[i], b, [])
end

fig = Figure()
ax1 = Axis(fig[1,1])
ax2 = Axis(fig[2,1])
lines!(ax1, f, amp)
lines!(ax2, f, phase)
current_figure()
##