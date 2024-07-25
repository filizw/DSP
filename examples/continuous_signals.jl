include("../src/DSP.jl")

using CairoMakie

t = range(-2, 2, 1000)

## rectangular
x = DSP.rectangular.(t, T=2)
lines(t, x)
## triangle
x = DSP.triangle.(t, T=2)
lines(t, x)
## letter M
x = DSP.letter_M.(t, T=2)
lines(t, x)
## letter U
x = DSP.letter_U.(t, T=2)
lines(t,x)
## ramp wave
x = DSP.ramp_wave.(t)
lines(t, x)
## sawtooth wave
x = DSP.sawtooth_wave.(t)
lines(t, x)
## triangular wave
x = DSP.triangular_wave.(t)
lines(t, x)
## square wave
x = DSP.square_wave.(t)
lines(t, x)
## pulse wave
x = DSP.pulse_wave.(t, 0.5)
lines(t, x)
## impulse repeater
f = DSP.impulse_repeater(DSP.triangular_wave, 0, 0.5)
x = f.(t)
lines(t, x)
## ramp wave bl
t = range(-1.5, 1.5, 1000)
x = DSP.ramp_wave_bl.(t, band = 10)

lines(t, x)
## sawtooth wave bl
t = range(-1.5, 1.5, 1000)
x = DSP.sawtooth_wave_bl.(t, band = 10)

lines(t, x)
## triangular wave bl
t = range(-1.5, 1.5, 1000)
x = DSP.triangular_wave_bl.(t, band = 10)

lines(t, x)
## square wave bl
t = range(-1.5, 1.5, 1000)
x = DSP.square_wave_bl.(t, band = 10)

lines(t, x)
## pulse wave bl
t = range(-1.5, 1.5, 1000)
x = DSP.pulse_wave_bl.(t, D = 0.5, band = 10)

lines(t, x)
## impulse repeater bl
fb = DSP.impulse_repeater_bl(DSP.ramp_wave, 0, 0.5, 10)

t = range(0, 1, 1000)
x = fb.(t)

lines(t, x)
## rand signal bl
rand_bl = DSP.rand_siganl_bl(30, 50)
t = range(0, 1, 1000)
x = rand_bl.(t)

lines(t, x)
##