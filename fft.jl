using FFTW
using GLMakie
GLMakie.activate!(inline = true)

# Number of points
N = 2^12 - 1
# Sample spacing
Ts = 1 / (1.1 * N)
# Sample rate
fs = 1 / Ts
# Start time
t0 = 0
tmax = t0 + N * Ts

# time coordinate
t = t0:Ts:tmax

# The underlying signal here is the sum of a sine wave at 60 cycles per second
# and its second harmonic (120 cycles per second) at half amplitude. We have
# discrete observations (samples) of this signal at each time `t`, with `fs`
# samples per second.

signal = sin.(2π * 60 * t) + .5 * sin.(2π * 120 * t)

F = fftshift(fft(signal))
freqs =  fftshift(fftfreq(length(t), fs))

fig = Figure(size=(800, 400))
ax1 = Axis(fig[1, 1], title = "Signal",
        limits=(0, 4/60, -1.5, 1.5),
        xlabel="time (s)")
spectrum = abs.(F)
ax2 = Axis(fig[1, 2], title = "Fourier Spectrum",
        limits=(0, 200, minimum(spectrum), maximum(spectrum)),
        xlabel="frequency (Hz)")
lines!(ax1, t, signal)
lines!(ax2, freqs, spectrum)
save("figures/1Dfft.png", fig)
display(fig)