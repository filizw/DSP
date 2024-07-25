include("../src/DSP.jl")

using CairoMakie

## fftfreq, amplitude spectrum, power spectrum, psd
fs = 100
N = 100
t = range(0, step = 1/fs, length = N)
x = 3*sin.(2*π*3.0.*t) .+ cos.(2*π*11.0.*t)

f = DSP.fftfreq(N, fs)

fig = Figure()
ax1 = Axis(fig[1,1])
ax2 = Axis(fig[2,1])
ax3 = Axis(fig[3,1])
lines!(ax1, f, DSP.amplitude_spectrum(x))
lines!(ax2, f, DSP.power_spectrum(x))
lines!(ax3, f, DSP.psd(x, DSP.rect(N), fs))
current_figure()
## periodogram
fs = 100
N = 100
rand_bl = DSP.rand_siganl_bl(10, 20)
t = range(0, step = 1/fs, length = N)
x = rand_bl.(t)

f = DSP.fftfreq(N, fs)
lines(f, DSP.psd(x, DSP.rect(N), fs))
## rfftfreq, stft, istft
fs = 20
t1 = range(0, step = 1/fs, length = 100)
t2 = range(t1[end] + 1/fs, step = 1/fs, length = 100)

x1 = sin.(2*π*2.0.*t1)
x2 = cos.(2*π*9.0.*t2)

x = vcat(x1, x2)

window_size = 50
window = DSP.hamming(window_size)
L = 25

X = DSP.stft(x, window, L)

xx = DSP.istft(X, window, L)

err = maximum(abs.(xx - x))
println("err = $err")

f = DSP.rfftfreq(window_size, fs)
heatmap(t, f, abs.(X))
##