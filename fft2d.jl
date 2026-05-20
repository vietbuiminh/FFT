using FFTW
using GLMakie
using FileIO
using Images

GLMakie.activate!(inline = true)

img = load("figures/image.png")
gray = Gray.(img)
signal = Float32.(gray)

Nx, Ny = size(signal)
dx = 1.0
dy = 1.0

# Nx = 256
# Ny = 256
# dx = 1 / (1.1 * Nx)
# dy = 1 / (1.1 * Ny)

fsx = 1 / dx
fsy = 1 / dy


x = (0:Nx-1) .* dx
y = (0:Ny-1) .* dy



# signal = [sin(2π*(60*xi + 120*yj)) for xi in x, yj in y] #plane wave
# signal = [sin(2π*60*xi) + 0.5*sin(2π*120*yj) for xi in x, yj in y] # separable

# 2D FFT + shift
F = fftshift(fft(signal))
spectrum = abs.(F)

fx = fftshift(fftfreq(Nx, fsx))
fy = fftshift(fftfreq(Ny, fsy))

function inverse_2Dfft(F)
    F_reconstructed = ifftshift(F)
    signal_reconstructed = ifft(F_reconstructed)
    real_reconstructed = real.(signal_reconstructed)
    return real_reconstructed
end

# filter and making mask
FX = repeat(fx, 1, Ny)
FY = repeat(fy', Nx, 1)
radius1 = 0.01
mask = (FX.^2 .+ FY.^2) .>= radius1^2
radius2 = 0.02
mask .&= (FX.^2 .+ FY.^2) .<= radius2^2

filtered_F = F .* mask
filtered_spectrum = abs.(filtered_F)

real_reconstructed = inverse_2Dfft(F)
filtered_reconstructed = inverse_2Dfft(filtered_F)

fig = Figure(size = (1600, 800))
ax1 = GLMakie.Axis(fig[1, 1], yreversed = true,title = "Input Signal (2D)", xlabel = "y", ylabel = "x")
ax2 = GLMakie.Axis(fig[1, 2], title = "Spectrum (2D)", xlabel = "fx (Hz)", ylabel = "fy (Hz)")
ax3 = GLMakie.Axis(fig[1, 3], title = "Log Spectrum (2D)", xlabel = "fx (log Hz)", ylabel = "fy (log Hz)")
ax4 = GLMakie.Axis(fig[1, 4], yreversed = true, title = "Reconstructed Signal (2D)", xlabel = "y", ylabel = "x")

ax6 = GLMakie.Axis(fig[2,2], title = "Filtered Spectrum r1,r2 = ($radius1,$radius2)", xlabel = "fx (Hz)", ylabel = "fy (Hz)")
ax7 = GLMakie.Axis(fig[2,3], title = "Filtered Log Spectrum", xlabel = "fx (Hz)", ylabel = "fy (Hz)")
ax8 = GLMakie.Axis(fig[2,4], yreversed = true, title = "Filtered Reconstructed Signal", xlabel = "y", ylabel = "x")

heatmap!(ax1, x, y, signal')  
heatmap!(ax2, fx, fy, spectrum)  
heatmap!(ax3, fx, fy, log.(1 .+ abs.(F))) 
heatmap!(ax4, fx, fy, real_reconstructed')  
heatmap!(ax6, fx, fy, filtered_spectrum)  
heatmap!(ax7, fx, fy, log.(1 .+ abs.(filtered_F)))  
heatmap!(ax8, fx, fy, filtered_reconstructed') 

display(fig)