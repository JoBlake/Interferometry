using FITSIO
using Statistics

# Check both initial and reconstructed images
println("="^60)
println("Checking initial_model.png (via FITS if available)...")
println("="^60)

# First check reconstructed
println("\nChecking reconstructed_image_1.fits...")
f = FITS("reconstructed_image_1.fits")
img = read(f[1])
close(f)

println("Image shape: ", size(img))
println("Min value: ", minimum(img))
println("Max value: ", maximum(img))
println("Mean value: ", mean(img))
println("Sum: ", sum(img))
println("Number of non-zero pixels: ", count(x -> x != 0, img))
println("Number of NaN pixels: ", count(isnan, img))
println("Number of Inf pixels: ", count(isinf, img))

# Show histogram of values
println("\nValue distribution:")
println("  Pixels > 0: ", count(x -> x > 0, img))
println("  Pixels == 0: ", count(x -> x == 0, img))
println("  Pixels < 0: ", count(x -> x < 0, img))

# Sample some pixel values
println("\nSample pixel values (center region):")
center = size(img, 1) ÷ 2
for i in -2:2
    for j in -2:2
        val = img[center + i, center + j]
        print("$(round(val, sigdigits=3)) ")
    end
    println()
end
