module DSP

using LinearAlgebra
using Primes

# Continuous signals
rectangular(t::Real; T::Real=1.0)::Real = abs(t) < T/2 ? 1 : 0
triangle(t::Real; T::Real=1.0)::Real = abs(t) < T/2 ? 1 - 2*abs(t)/T : 0
letter_M(t::Real; T::Real=1.0)::Real = abs(t) < T/2 ? 2*abs(t)/T : 0
letter_U(t::Real; T::Real=1.0)::Real = abs(t) < T/2 && abs(t) > 3*T/8 ? 1 : 0

ramp_wave(t::Real)::Real = 2*(t - floor(t - 0.5) - 1)
sawtooth_wave(t::Real)::Real = ramp_wave(-t)
triangular_wave(t::Real)::Real = 2*(triangle(t - floor(t + 0.25) - 0.25) - 0.5)
square_wave(t::Real)::Real = (t - floor(t)) < 0.5 ? 1 : -1
pulse_wave(t::Real, ρ::Real=0.2)::Real = rectangular(t - floor(t) - ρ/2, T = ρ)
impulse_repeater(g::Function, t1::Real, t2::Real)::Function = f(t::Real)::Real = g((t - t1) % (t2 - t1) + (t >= t1 ? t1 : t2))

function ramp_wave_bl(t::Real; A::Real=1.0, T::Real=1.0, band::Real=20.0)::Real
    y::Real = 0
    for k in range(1, floor(band*T))
        y -= (2*A/(π*k))*sin(2*π*k*(t - T/2)/T)
    end

    return y
end

function sawtooth_wave_bl(t::Real; A::Real=1.0, T::Real=1.0, band::Real=20.0)::Real
    y::Real = 0
    for k in range(1, floor(band*T))
        y += (2*A/(π*k))*sin(2*π*k*(t - T/2)/T)
    end

    return y
end

function triangular_wave_bl(t::Real; A::Real=1.0, T::Real=1.0, band::Real=20.0)::Real
    y::Real = 0
    for k in range(1, floor(band*T), step = 2)
        y += (4*A*(1 - cos(π*k))/((π*k)^2))*cos(2*pi*k*(t - T/4)/T)
    end

    return y
end

function square_wave_bl(t::Real; A::Real=1.0, T::Real=1.0, band::Real=20.0)
    y::Real = 0
    for k in range(1, floor(band*T), step = 2)
        y += (2*A*(1 - cos(π*k))/(π*k))*sin(2*π*k*t/T)
    end

    return y
end

function pulse_wave_bl(t::Real; D::Real=0.2, A::Real=1.0, T::Real=1.0, band::Real=20.0)
    y::Real = A*D
    for k in range(1, floor(band*T))
        y += (A*sin(2*π*k*D)/(k*π))*cos(2*π*k*t/T)
        y -= (A*(cos(2*π*k*D) - 1)/(k*π))*sin(2*π*k*t/T)
    end

    return y
end

function impulse_repeater_bl(g::Function, t1::Real, t2::Real, band::Real)::Function
    T₀ = t2 - t1
    f₀ = 1/T₀
    ω₀ = 2*π*f₀
    K = convert(UInt64, floor(band/f₀))

    dt = 0.1/band
    t = range(t1, t2, step = dt)

    a₀ = (1/T₀)*sum(g.(t))*dt

    a = zeros(K)
    b = zeros(K)
    for k in range(1, K)
        a[k] = (2/T₀)*sum(g.(t).*cos.(ω₀*k.*t))*dt
        b[k] = (2/T₀)*sum(g.(t).*sin.(ω₀*k.*t))*dt
    end

    function fb(t::Real)::Real
        y = a₀
        for k in range(1, K)
            y += a[k]*cos(ω₀*k*t) + b[k]*sin(ω₀*k*t)
        end

        return y
    end

    return fb
end

function rand_siganl_bl(f1::Real, f2::Real)::Function
    N = 100
    freq = range(f1, f2, N)
    A = rand(N)
    ϕ = 2*π*rand(N)

    f(t::Real) = sum(A.*cos.(2*π.*freq*t .+ ϕ))

    return f
end

# Discrete signals
kronecker(n::Integer)::Real = n == 0 ? 1 : 0
heaviside(n::Integer)::Real = n >= 0 ? 1 : 0

# Windows
rect(N::Integer)::AbstractVector{<:Real} = ones(N)
triang(N::Integer)::AbstractVector{<:Real} = 1 .- 2*abs.(range(0, N-1) .- (N-1)/2)/(N-1)
hanning(N::Integer)::AbstractVector{<:Real} = 0.5*(1 .- cos.(2*π.*range(0, N-1)/(N-1)))
hamming(N::Integer)::AbstractVector{<:Real} = 0.54 .- 0.46*cos.(2*π.*range(0, N-1)/(N-1))
blackman(N::Integer)::AbstractVector{<:Real} = 0.42 .- 0.5*cos.(2*π.*range(0, N-1)/(N-1)) + 0.08*cos.(4*π.*range(0, N-1)/(N-1))

# Signal parameters
mean(x::AbstractVector)::Real = sum(x)/length(x)
peak2peak(x::AbstractVector)::Real = maximum(x) - minimum(x)
energy(x::AbstractVector)::Real = sum(abs2.(x))
power(x::AbstractVector)::Real = energy(x)/length(x)
rms(x::AbstractVector)::Real = sqrt(power(x))

function running_mean(x::AbstractVector, M::Integer)::Vector
    N = length(x)
    y = zeros(N)

    for n in range(1, N)
        n₁ = max(n-M, 1)
        n₂ = min(n+M, N)

        y[n] = mean(x[n₁:n₂])
    end

    return y
end

function running_energy(x::AbstractVector, M::Integer)::Vector
    N = length(x)
    y = zeros(N)

    for n in range(1, N)
        n₁ = max(n-M, 1)
        n₂ = min(n+M, N)

        y[n] = energy(x[n₁:n₂])
    end

    return y
end

function running_power(x::AbstractVector, M::Integer)::Vector
    N = length(x)
    y = zeros(N)

    for n in range(1, N)
        n₁ = max(n-M, 1)
        n₂ = min(n+M, N)

        y[n] = power(x[n₁:n₂])
    end

    return y
end

# Sampling
function interpolate(t::Real; m::AbstractVector, s::AbstractVector, kernel::Function=sinc)::Real
    T = m[2] - m[1]
    y::Real = 0
    for n in range(1, length(m))
        y += s[n]*kernel((t - m[n])/T)
    end

    return y
end

# Quantization
quantize(x::Real; L::AbstractVector)::Real = L[argmin(abs.(x .- L))]
SQNR(N::Integer)::Real = 6.02*N + 1.76
SNR(Psignal, Pnoise)::Real = 10*log10(Psignal/Pnoise)

# DFT
function dtft(f::Real; signal::AbstractVector, fs::Real)
    N = length(signal)
    X = 0
    for n in 0:(N-1)
        X += signal[n+1]*exp(-im*2*π*f*n/fs)
    end

    return X
end

function dft(x::AbstractVector)::Vector
    N = length(x)
    X = zeros(ComplexF64, N)
    for k in range(0, N-1)
        for n in range(0, N-1)
            X[k+1] += x[n+1]*exp(-im*2*π*k*n/N)
        end
    end

    return X
end

function idft(X::AbstractVector)::Vector
    K = length(X)
    x = zeros(ComplexF64, K)
    for n in range(0, K-1)
        for k in range(0, K-1)
            x[n+1] += X[k+1]*exp(im*2*π*k*n/K)
        end
    end

    return x/K
end

function rdft(x::AbstractVector)::Vector
    N = length(x)
    K = convert(UInt64, floor(N/2 + 1))
    X = zeros(ComplexF64, K)
    for k in range(0, K-1)
        for n in range(0, N-1)
            X[k+1] += x[n+1]*exp(-im*2*π*k*n/N)
        end
    end

    return X
end

function irdft(X::AbstractVector, N::Integer)::Vector
    K = length(X)
    XX = vcat(X, conj.(X[(K-(N%2 == 0 ? 1 : 0)):-1:2]))

    return idft(XX)
end

# FFT
function fft_radix2_dit_r(x::AbstractVector)::Vector
    N::UInt64 = length(x)

    if N == 1
        return x
    elseif N == 2
        return [x[1] + x[2], x[1] - x[2]]
    else
        X_even::Vector{ComplexF64} = fft_radix2_dit_r(x[1:2:end])
        X_odd::Vector{ComplexF64} = fft_radix2_dit_r(x[2:2:end])

        X::Vector{ComplexF64} = zeros(ComplexF64, N)
        for i::UInt64 in range(0, N-1)
            m::UInt64 = i%(N/2) + 1
            X[i+1] = X_even[m] + exp(-im*2*π*i/N)*X_odd[m]
        end
    end

    return X
end

function ifft_radix2_dif_r(X::AbstractVector)::Vector
    K::UInt64 = length(X)

    if K == 1
        return X
    elseif K == 2
        return [X[1] + X[2], X[1] - X[2]]/2
    else
        x_even::Vector{ComplexF64} = ifft_radix2_dif_r(X[1:2:end])
        x_odd::Vector{ComplexF64} = ifft_radix2_dif_r(X[2:2:end])

        x::Vector{ComplexF64} = zeros(ComplexF64, K)
        for i::UInt64 in range(0, K-1)
            m::UInt64 = i%(K/2) + 1
            x[i+1] = x_even[m] + exp(im*2*π*i/K)*x_odd[m]
        end
    end

    return x/2
end

function fft(x::AbstractVector)::Vector
    N::UInt64 = length(x)

    if N == 1
        return x
    elseif isinteger(log2(N))
        return fft_radix2_dit_r(x)
    else
        factors::Vector{UInt64} = factor(Vector{UInt64}, N)
        factors_len::UInt64 = length(factors)

        X = convert(Vector{ComplexF64}, x)
        n1::UInt64 = factors[1]

        coprime::Bool = false

        if factors_len != 1
            n2::UInt64 = prod(factors[2:end])
            coprime = gcdx(n1, n2)[1] == 1

            if coprime
                X = zeros(ComplexF64, n1, n2)
                for i in range(0, N-1)
                    X[i%n1 + 1, i%n2 + 1] = x[i+1]
                end
            else
                X = transpose(reshape(X, (n2, n1)))
            end
        end

        w::Vector{ComplexF64} = exp.(-im*2*pi.*range(0, n1-1)/n1)
        W1::Matrix{ComplexF64} = ones(ComplexF64, n1, n1)
        for j in range(1, n1-1), i in range(1, n1-1)
            W1[i+1, j+1] = w[(i*j)%n1 + 1]
        end

        X = W1*X

        if factors_len == 1
            return X
        end

        if !coprime
            w = exp.(-im*2*pi.*range(0, N-1)/N)
            W::Matrix{ComplexF64} = ones(ComplexF64, n1, n2)
            for j in range(1, n2-1), i in range(1, n1-1)
                W[i+1, j+1] = w[(i*j)%N + 1]
            end

            X .*= W
        end

        if factors_len == 2
            if n1 != n2
                w = exp.(-im*2*pi.*range(0, n2-1)/n2)
                W1 = ones(ComplexF64, n2, n2)
                for j in range(1, n2-1), i in range(1, n2-1)
                    W1[i+1, j+1] = w[(i*j)%n2 + 1]
                end
            end

            X *= W1
        else
            for i in range(1, n1)
                X[i, :] = fft(X[i, :])
            end
        end

        if coprime
            XX::Vector{ComplexF64} = zeros(ComplexF64, N)
            for i in range(0, n1-1)
                XX[((i*n2)%n1 + 1):n1:end] = circshift(X[i+1, :], floor(i*n2/n1))
            end

            return XX
        else
            return X[:]
        end
    end
end

function ifft(X::AbstractVector)::Vector
    K = length(X)

    if K == 1
        return X
    else
        fact = factor(Vector{UInt64}, K)

        if length(fact) == 1
            return idft(X)
        else
            k1 = fact[1]
            k2 = prod(fact[2:end])

            w = exp(2*pi*im/K)

            x = transpose(convert(Matrix{ComplexF64}, reshape(copy(X), (k2, k1))))
            W = zeros(ComplexF64, k1, k2)
            for i in range(0, k2-1)
                x[:, i+1] = idft(x[:, i+1])
                for j in range(0, k1-1)
                    W[j+1, i+1] = w^(i*j)
                end
            end

            x .*= W
            for i in range(1, k1)
                x[i, :] = ifft(x[i, :])
            end

            return x[:]
        end
    end
end

# Frequency analysis
fftfreq(N::Integer, fs::Real)::Vector = range(0, (fs - fs/N), N)
rfftfreq(N::Integer, fs::Real)::Vector = range(0, fs/N*(floor(N/2 + 1)), floor(UInt64, N/2 + 1))
amplitude_spectrum(x::AbstractVector, w::AbstractVector=rect(length(x)))::Vector = abs.(fft(x.*w))/sum(w)
power_spectrum(x::AbstractVector, w::AbstractVector=rect(length(x)))::Vector = amplitude_spectrum(x, w).^2
psd(x::AbstractVector, w::AbstractVector=rect(length(x)), fs::Real=1.0)::Vector = power_spectrum(x, w)/(fs*sum(abs2.(w))/abs2(sum(w)))

function periodogram(x::AbstractVector, w::AbstractVector=rect(length(x)), fs::Real=1.0)::Vector
    K::UInt64 = length(w)
    L::UInt64 = K/2
    M::UInt64 = floor((length(x)-K)/(K-L) + 1)

    Pm = zeros(K)

    for i in range(0, M-1)
        Pm += abs2.(fft(x[(i*(K-L) + 1):(K + i*(K-L))].*w))/sum(w.^2)
    end

    Pm /= (M*fs)

    return Pm
end

function stft(x::AbstractVector, w::AbstractVector, L::Integer)::Matrix
    N = length(x)
    K = length(w)
    M = convert(UInt64, floor((N-L)/(K-L)))
    X = zeros(ComplexF64, M, convert(UInt64, floor(K/2 + 1)))

    for m in range(0, M-1)
        n₁ = m*(K-L) + 1
        n₂ = m*(K-L) + K

        X_seg = rdft(x[n₁:n₂].*w)
        X[m+1, :] = X_seg
    end

    return X
end

function istft(X::AbstractMatrix{<:Complex}, w::AbstractVector{<:Real}, L::Integer)::AbstractVector{<:Real}
    K = length(w)
    M = size(X, 1)
    N = M*(K-L) + L
    x = zeros(N)

    for m in range(1, M)
        n₁ = (m-1)*(K-L) + 1
        n₂ = m*(K-L) + L

        x_seg = irdft(X[m, :], K)./w
        x[n₁:n₂] = real.(x_seg)
    end

    return x
end

# LTI systems
function conv(f::Vector, g::Vector)::Vector
    N = length(f)
    M = length(g)
    K = N+M-1

    y = zeros(K)

    for n in range(0, K-1)
        for m in range(0, N-1)
            if n-m >= 0 && n-m < M
                y[n+1] += f[m+1]*g[n-m+1]
            end
        end
    end

    return y
end

function fast_conv(f::Vector, g::Vector)::Vector
    N = length(f)
    M = length(g)
    K = convert(UInt64, 2^ceil(log2(N+M-1)))

    f_zp = vcat(f, zeros(K-N))
    g_zp = vcat(g, zeros(K-M))

    F = fft_radix2_dit_r(f_zp)
    G = fft_radix2_dit_r(g_zp)

    y = ifft_radix2_dif_r(F.*G)

    return y[1:(N+M-1)]
end

function overlap_add(x::Vector, h::Vector, L::Integer)::Vector
    N = length(x)
    M = length(h)
    K = L+M-1

    NL = convert(UInt64, ceil(N/L))

    h_zp = vcat(h, zeros(K-M))
    H = fft_radix2_dit_r(h_zp)

    y = zeros(N+M-1)

    zp = K-L
    for i in range(1, NL)
        n1 = (i-1)*L + 1
        n2 = i*L
        n3 = n2+M-1

        if n2 > N
            n2 = N
            zp = K - (n2 - n1 + 1)
        end

        if n3 > (N+M-1)
            n3 = N+M-1
        end

        x_seg = vcat(x[n1:n2], zeros(zp))
        X = fft_radix2_dit_r(x_seg)
        y_seg = ifft_radix2_dif_r(X.*H)

        y[n1:n3] += real.(y_seg[1:(n3 - n1 + 1)])
    end

    return y
end

function overlap_save(x::Vector, h::Vector, L::Integer)::Vector
    N = length(x)
    M = length(h)
    K = L+M-1

    NL = convert(UInt64, ceil(N/L)) + 1

    h_zp = vcat(h, zeros(K-M))
    H = fft_radix2_dit_r(h_zp)

    y = zeros(N+M-1)

    x_seg = zeros(K)
    zp = 0
    l = 0
    for i in range(1, NL)
        n1 = (i-1)*L + 1
        n2 = n3 = i*L

        if n2 > N
            n2 = N
            zp = n1 > n2 ? L : L-(n2-n1+1)
        end

        if n3 > (N+M-1)
            n3 = N+M-1
            l = L-(n3-n1+1)
        end

        x_seg = vcat(x_seg[(L+1):end], x[n1:n2], zeros(zp))
        X = fft_radix2_dit_r(x_seg)
        y_seg = ifft_radix2_dif_r(X.*H)

        y[n1:n3] = real.(y_seg[(K-L+1):(K-l)])
    end

    return y
end

function lti_filter(b::Vector, a::Vector, x::Vector)::Vector
    N = length(x)
    M = length(b)
    K = length(a)

    y = zeros(N)
    for n in range(0, N-1)
        for m in range(0, M-1)
            if n-m >= 0 && n-m < N
                y[n+1] += b[m+1]*x[n-m+1]
            end
        end

        for k in range(1, K-1)
            if n-k >= 0 && n-k < N
                y[n+1] -= a[k+1]*y[n-k+1]
            end
        end
    end

    return y
end

function filtfilt(b::Vector, a::Vector, x::Vector)::Vector
    y = lti_filter(b, a, x)
    y = y[end:(-1):1]
    y = lti_filter(b, a, y)
    y = y[end:(-1):1]

    return y
end

function lti_amp(f::Real, b::Vector, a::Vector)::Real
    z = exp(-im*2*pi*f)

    Y = evalpoly(z, b)
    X = evalpoly(z, vcat([1], a))

    return abs(Y/X)
end

function lti_phase(f::Real, b::Vector, a::Vector)::Real
    z = exp(-im*2*pi*f)

    Y = evalpoly(z, b)
    X = evalpoly(z, vcat([1], a))

    return angle(Y/X)
end

# FIR filters
function firwin_lp_I(order::Integer, f0::Real)
    m = range(-order/2, order/2)
    h = 2*f0*sinc.(2*f0.*m)

    return h
end

function firwin_hp_I(order::Integer, f0::Real)
    m = range(-order/2, order/2)
    h = -2*f0*sinc.(2*f0.*m)
    h[convert(Int64, order/2) + 1] += 1

    return h
end

function firwin_bp_I(order::Integer, f1::Real, f2::Real)
    m = range(-order/2, order/2)
    h = -2*f1*sinc.(2*f1.*m) .+ 2*f2*sinc.(2*f2.*m)

    return h
end

function firwin_bs_I(order::Integer, f1::Real, f2::Real)
    m = range(-order/2, order/2)
    h = 2*f1*sinc.(2*f1.*m) .- 2*f2*sinc.(2*f2.*m)
    h[convert(Int64, order/2) + 1] += 1

    return h
end

function firwin_lp_II(order::Integer, f0::Real)
    m = range(-order/2, order/2)
    h = 2*f0*sinc.(2*f0.*m)

    return h
end

function firwin_bp_II(order::Integer, f1::Real, f2::Real)
    m = range(-order/2, order/2)
    h = -2*f1*sinc.(2*f1.*m) .+ 2*f2*sinc.(2*f2.*m)

    return h
end

function firwin_diff(order::Integer)
    m = range(-order/2, order/2)
    h = ((-1).^m)./m
    h[convert(Int64, order/2) + 1] = 0

    return h
end

# Resampling
function resample(x::Vector, M::Integer, N::Integer)
    g = zeros(ComplexF64, length(x)*M)
    g[1:M:end] = x
    g = g[1:N:end]

    return g
end

end