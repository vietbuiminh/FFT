using FFTW
using GLMakie
using FileIO
using Images
using Printf
GLMakie.activate!(inline = true)

colormap = :grays

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

function inverse_2Dfft(F)
    F_reconstructed = ifftshift(F)
    signal_reconstructed = ifft(F_reconstructed)
    real_reconstructed = real.(signal_reconstructed)
    return real_reconstructed
end

fx = fftshift(fftfreq(Nx, fsx))
fy = fftshift(fftfreq(Ny, fsy))

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
ax1 = GLMakie.Axis(fig[1, 1], yreversed = true, title = "Input Signal (2D)", xlabel = "y", ylabel = "x")
ax2 = GLMakie.Axis(fig[1, 2], title = "Spectrum (2D)", xlabel = "fx (Hz)", ylabel = "fy (Hz)")
ax3 = GLMakie.Axis(fig[1, 3], title = "Log Spectrum (2D)", xlabel = "fx (log Hz)", ylabel = "fy (log Hz)")
ax4 = GLMakie.Axis(fig[1, 4], yreversed = true, title = "Reconstructed Signal (2D)", xlabel = "y", ylabel = "x")

ax6 = GLMakie.Axis(fig[2,2], title = "Filtered Spectrum r1,r2 = ($radius1,$radius2)", xlabel = "fx (Hz)", ylabel = "fy (Hz)")
ax7 = GLMakie.Axis(fig[2,3], title = "Log Filtered Spectrum", xlabel = "fx (Hz)", ylabel = "fy (Hz)")
ax8 = GLMakie.Axis(fig[2,4], yreversed = true, title = "Filtered Reconstructed Signal", xlabel = "y", ylabel = "x")

heatmap!(ax1, x, y, signal', colormap = colormap)
heatmap!(ax2, fx, fy, spectrum, colormap = colormap)
heatmap!(ax3, fx, fy, log.(1 .+ abs.(F)), colormap = colormap)
heatmap!(ax4, fx, fy, real_reconstructed', colormap = colormap)
heatmap!(ax6, fx, fy, filtered_spectrum, colormap = colormap)
heatmap!(ax7, fx, fy, log.(1 .+ abs.(filtered_F)), colormap = colormap)
heatmap!(ax8, fx, fy, filtered_reconstructed', colormap = colormap)

record(fig, "figures/fft2d_animation.gif", 1:20, framerate = 8) do i
    global ax6, ax7, ax8, radius1, radius2
    GLMakie.empty!(ax6)
    GLMakie.empty!(ax7)
    GLMakie.empty!(ax8)
    # ax6 = GLMakie.Axis(fig[2,2], title = "Filtered Spectrum r1,r2 = ($radius1,$radius2)", xlabel = "fx (Hz)", ylabel = "fy (Hz)")
    # ax7 = GLMakie.Axis(fig[2,3], title = "Filtered Log Spectrum", xlabel = "fx (Hz)", ylabel = "fy (Hz)")
    # ax8 = GLMakie.Axis(fig[2,4], yreversed = true, xreversed = true, title = "Filtered Reconstructed Signal", xlabel = "y", ylabel = "x")

    radius1 = 0.05 - i / 500
    radius2 = 0.22 - i / 100
    mask = (FX.^2 .+ FY.^2) .>= radius1^2
    mask .&= (FX.^2 .+ FY.^2) .<= radius2^2
    filtered_F = F .* mask
    filtered_spectrum = abs.(filtered_F)
    filtered_reconstructed = inverse_2Dfft(filtered_F)
    heatmap!(ax6, fx, fy, filtered_spectrum, colormap = colormap)
    ax6.title = "Filtered Spectrum r1,r2 = ($(@sprintf("%.2f",(round(radius1, digits=3)))),$(@sprintf("%.2f",(round(radius2, digits=3)))))"
    heatmap!(ax7, fx, fy, log.(1 .+ abs.(filtered_F)), colormap = colormap)
    heatmap!(ax8, fx, fy, filtered_reconstructed', colormap = colormap)
end

display(fig)