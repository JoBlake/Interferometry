# Quick analysis of the reconstructed image
using FITSIO
using Plots
using Statistics

println("="^70)
println("RECONSTRUCTED IMAGE ANALYSIS")
println("="^70)

# Prompt user for FITS file
println("\nEnter the path to your FITS file:")
println("(Press Enter to use: reconstructed_image.fits)")
print("> ")

user_input = readline()

if isempty(strip(user_input))
    # Use default
    filename = "reconstructed_image.fits"
    println("Using default file: $filename")
else
    filename = strip(user_input)
end

# Check if file exists
if !isfile(filename)
    println("\n⚠ Error: File not found: $filename")
    println("\nAvailable .fits files in current directory:")
    fits_files = filter(f -> endswith(f, ".fits"), readdir("."))
    if isempty(fits_files)
        println("  (none found)")
    else
        for f in fits_files
            println("  - $f")
        end
    end
    println("\nPress Enter to exit...")
    readline()
    exit(1)
end

println("\nLoading $filename...")

# Read the reconstructed FITS file
if isfile(filename)
    f = FITS(filename)
    img = read(f[1])
    close(f)

    println("\nImage statistics:")
    println("  Size: $(size(img))")
    println("  Min value: $(minimum(img))")
    println("  Max value: $(maximum(img))")
    println("  Mean value: $(mean(img))")
    println("  Std dev: $(std(img))")
    println("  Total flux: $(sum(img))")

    # Find peak location
    max_idx = argmax(img)
    println("\n  Peak location: $max_idx")
    println("  Peak value: $(img[max_idx])")

    # Check concentration
    flux_in_peak_pixel = img[max_idx] / sum(img) * 100
    println("  Flux in peak pixel: $(round(flux_in_peak_pixel, digits=2))%")

    # Count significant pixels (> 1% of peak)
    threshold = 0.01 * maximum(img)
    n_significant = sum(img .> threshold)
    println("  Pixels > 1% of peak: $n_significant ($(round(n_significant/length(img)*100, digits=2))%)")

    # Check if image is centered
    center = size(img) .÷ 2
    center_region = img[center[1]-10:center[1]+10, center[2]-10:center[2]+10]
    center_flux = sum(center_region) / sum(img) * 100
    println("\n  Flux in central 21×21 pixels: $(round(center_flux, digits=1))%")

    # Check for multiple peaks
    println("\n  Looking for multiple bright regions...")
    sorted_vals = sort(vec(img), rev=true)
    top_10_flux = sum(sorted_vals[1:10]) / sum(img) * 100
    top_100_flux = sum(sorted_vals[1:100]) / sum(img) * 100
    println("  Top 10 pixels contain: $(round(top_10_flux, digits=1))% of flux")
    println("  Top 100 pixels contain: $(round(top_100_flux, digits=1))% of flux")

    if top_10_flux > 90
        println("\n  ⚠ VERY CONCENTRATED - Nearly point-like or has strong peaks")
    elseif top_100_flux < 30
        println("\n  ⚠ VERY DIFFUSE - Flux spread across many pixels")
    else
        println("\n  ✓ Moderate concentration - Extended structure")
    end

    # Create contour plot
    println("\nCreating contour plot...")
    levels = [0.01, 0.05, 0.1, 0.2, 0.4, 0.6, 0.8, 0.95] .* maximum(img)
    p = contour(img',
                levels=levels,
                xlabel="Pixel X",
                ylabel="Pixel Y",
                title="Brightness Contours - $filename",
                aspect_ratio=:equal,
                fill=true,
                c=:hot)

    # Generate output filename based on input filename
    base_name = replace(filename, ".fits" => "")
    contour_file = base_name * "_contours.png"
    savefig(p, contour_file)
    println("Saved: $contour_file")

    # Create radial profile
    println("\nCalculating radial profile...")
    center = size(img) .÷ 2
    x = repeat(1:size(img,1), 1, size(img,2)) .- center[1]
    y = repeat((1:size(img,2))', size(img,1), 1) .- center[2]
    r = sqrt.(x.^2 .+ y.^2)

    max_r = minimum(size(img)) ÷ 2
    r_bins = 0:max_r
    profile = zeros(length(r_bins))

    for i in 1:length(r_bins)
        if i < length(r_bins)
            mask = (r .>= r_bins[i]) .& (r .< r_bins[i+1])
            if sum(mask) > 0
                profile[i] = mean(img[mask])
            end
        end
    end

    p2 = plot(r_bins[1:end-1], profile[1:end-1],
              xlabel="Radius (pixels)",
              ylabel="Mean brightness",
              title="Radial Profile - $filename",
              legend=false,
              lw=2)

    profile_file = base_name * "_radial_profile.png"
    savefig(p2, profile_file)
    println("Saved: $profile_file")
end

println("\n" * "="^70)
println("Press Enter to exit...")
readline()
