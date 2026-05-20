# make sure to run either fft2d.jl or fft2d_animate.jl before running this file, as it relies on the variables defined there
# I will use spectrum and F for this 
# include("fft2d.jl")

function generate_2d_cosine(freq_x, freq_y, x, y)
    image =  [cos(2π * (freq_x * xi + freq_y * yi)) for xi in x, yi in y]
    return image
end

mag = spectrum # because spectrum is abs.(F)
mag_vec = vec(mag)
idx = sortperm(vec(mag); rev = true)
mid = div(size(F, 1), 2) + 1
center_idx = LinearIndices(F)[mid, mid]
idx = filter(i -> i != center_idx, idx)
max_coeff_mag = maximum(mag_vec)
F_sorted = zeros(ComplexF64, size(F))
recons = Vector{Array{Float64, 2}}()

# Only process every other index (symmetric frequencies share the same magnitude)
idx_subset = idx[1:2:1000]

for i in idx_subset
    F_sorted[i] = F[i]
    push!(recons, inverse_2Dfft(F_sorted))
end

cart = CartesianIndices(size(F))
first_ix, first_iy = Tuple(cart[idx_subset[1]])
freq_x = fx[first_ix]
freq_y = fy[first_iy]
cosine_image = generate_2d_cosine(freq_x, freq_y, x, y; amp = max_coeff_mag > 0 ? 1.0 : 0.0)

ani = Figure(size = (800, 400))
ax_recon = GLMakie.Axis(ani[1, 2], yreversed = true, title = "Reconstruction FFT Coefficients", xlabel = "y", ylabel = "x")
ax_wave = GLMakie.Axis(ani[1, 1], yreversed = true, title = "Cosine Wave", xlabel = "y", ylabel = "x")
heatmap!(ax_wave, x, y, cosine_image', colormap = :grays, colorrange = (-1, 1))

record(ani, "figures/fft2d_reconstruction.gif", 1:length(idx_subset), framerate = 30) do i
    GLMakie.empty!(ax_recon)
    GLMakie.empty!(ax_wave)
    ix, iy = Tuple(cart[idx_subset[i]])
    freq_x = fx[ix]
    freq_y = fy[iy]
    coeff_mag = mag_vec[idx_subset[i]]
    cosine_image = generate_2d_cosine(freq_x, freq_y, x, y)
    ax_recon.title = "Reconstruction (|F|=$(round(coeff_mag, digits=1)))"
    heatmap!(ax_recon, x, y, recons[i]', colormap = :grays)
    heatmap!(ax_wave, x, y, cosine_image', colormap = :grays, colorrange = (-1, 1))
end
display(ani)