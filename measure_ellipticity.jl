# Measure ellipticity and size of reconstructed object
using FITSIO
using Statistics
using LinearAlgebra

println("="^70)
println("OBJECT SIZE AND ELLIPTICITY MEASUREMENT")
println("="^70)

# Read the reconstructed FITS file
if !isfile("reconstructed_image.fits")
    println("\n⚠ Error: reconstructed_image.fits not found!")
    println("Press Enter to exit...")
    readline()
    exit()
end

f = FITS("reconstructed_image.fits")
img = read(f[1])
header = read_header(f[1])
close(f)

# Get pixel scale from FITS header (in degrees)
pixscale_deg = abs(get(header, "CDELT1", 0.0))  # degrees/pixel
pixscale_mas = pixscale_deg * 3600000  # milliarcseconds/pixel
pixscale_arcsec = pixscale_deg * 3600  # arcseconds/pixel

println("\nImage parameters:")
println("  Size: $(size(img))")
println("  Pixel scale: $(round(pixscale_mas, digits=4)) mas/pixel")
println("  Pixel scale: $(round(pixscale_arcsec, digits=6)) arcsec/pixel")

# Calculate moments to find centroid and axes
npix = size(img, 1)
x_coords = repeat(1:npix, 1, npix)
y_coords = repeat((1:npix)', npix, 1)

# Centroid
total_flux = sum(img)
x_center = sum(x_coords .* img) / total_flux
y_center = sum(y_coords .* img) / total_flux

println("\nCentroid:")
println("  X: $(round(x_center, digits=2)) pixels")
println("  Y: $(round(y_center, digits=2)) pixels")

# Second moments (for ellipse fitting)
x_shifted = x_coords .- x_center
y_shifted = y_coords .- y_center

Mxx = sum(x_shifted.^2 .* img) / total_flux
Myy = sum(y_shifted.^2 .* img) / total_flux
Mxy = sum(x_shifted .* y_shifted .* img) / total_flux

# Eigenvalues of moment matrix give major/minor axes
moment_matrix = [Mxx Mxy; Mxy Myy]
eigenvals = eigvals(moment_matrix)
sort!(eigenvals, rev=true)  # Largest first

# Convert eigenvalues to FWHM (assuming Gaussian-like)
# For 2D Gaussian: sigma^2 = second moment, FWHM = 2.355*sigma
major_axis_sigma = sqrt(eigenvals[1])
minor_axis_sigma = sqrt(eigenvals[2])

major_axis_fwhm_pix = 2.355 * major_axis_sigma
minor_axis_fwhm_pix = 2.355 * minor_axis_sigma

# Convert to angular size
major_axis_mas = major_axis_fwhm_pix * pixscale_mas
minor_axis_mas = minor_axis_fwhm_pix * pixscale_mas
major_axis_arcsec = major_axis_fwhm_pix * pixscale_arcsec
minor_axis_arcsec = minor_axis_fwhm_pix * pixscale_arcsec

# Ellipticity
ellipticity = major_axis_fwhm_pix / minor_axis_fwhm_pix

# Position angle (from eigenvectors)
eigenvecs = eigvecs(moment_matrix)
major_axis_vec = eigenvecs[:, 1]  # Eigenvector of largest eigenvalue
pa_rad = atan(major_axis_vec[2], major_axis_vec[1])
pa_deg = rad2deg(pa_rad)

println("\n" * "="^70)
println("MEASURED OBJECT PROPERTIES")
println("="^70)

println("\nSize (FWHM):")
println("  Major axis: $(round(major_axis_mas, digits=3)) mas = $(round(major_axis_arcsec, digits=6)) arcsec")
println("  Minor axis: $(round(minor_axis_mas, digits=3)) mas = $(round(minor_axis_arcsec, digits=6)) arcsec")
println("  Mean diameter: $(round((major_axis_mas + minor_axis_mas)/2, digits=3)) mas")

println("\nEllipticity:")
println("  Major/Minor ratio: $(round(ellipticity, digits=3))")
println("  Ellipticity (1 - b/a): $(round(1 - 1/ellipticity, digits=3))")

println("\nPosition angle:")
println("  PA: $(round(pa_deg, digits=1))° (East of North)")

# Alternative measurement: Half-light radius
sorted_pixels = sort(vec(img), rev=true)
cumsum_flux = cumsum(sorted_pixels)
half_light_idx = findfirst(cumsum_flux .>= total_flux / 2)
half_light_pixels = half_light_idx

println("\nHalf-light properties:")
println("  Half-light contained in: $half_light_pixels pixels")
println("  Effective radius: $(round(sqrt(half_light_pixels/π), digits=1)) pixels")
println("  Effective diameter: $(round(2*sqrt(half_light_pixels/π)*pixscale_mas, digits=3)) mas")

# Contour-based measurement (at 50% of peak)
threshold = 0.5 * maximum(img)
mask = img .> threshold
n_pixels_above = sum(mask)

if n_pixels_above > 0
    # Find extent in X and Y
    rows_with_flux = any(mask, dims=2)[:]
    cols_with_flux = any(mask, dims=1)[:]

    y_extent = sum(rows_with_flux)
    x_extent = sum(cols_with_flux)

    x_size_mas = x_extent * pixscale_mas
    y_size_mas = y_extent * pixscale_mas
    x_size_arcsec = x_extent * pixscale_arcsec
    y_size_arcsec = y_extent * pixscale_arcsec

    println("\nSize at 50% contour:")
    println("  X extent: $(round(x_size_mas, digits=3)) mas = $(round(x_size_arcsec, digits=6)) arcsec")
    println("  Y extent: $(round(y_size_mas, digits=3)) mas = $(round(y_size_arcsec, digits=6)) arcsec")
    println("  Aspect ratio (X/Y): $(round(x_extent/y_extent, digits=3))")
end

# Summary
println("\n" * "="^70)
println("SUMMARY")
println("="^70)
println("\nMost reliable measurement (moment-based FWHM):")
println("  Major axis: $(round(major_axis_mas, digits=2)) mas ($(round(major_axis_arcsec, digits=5)) arcsec)")
println("  Minor axis: $(round(minor_axis_mas, digits=2)) mas ($(round(minor_axis_arcsec, digits=5)) arcsec)")
println("  Ellipticity (major/minor): $(round(ellipticity, digits=2))")
println("  Position angle: $(round(pa_deg, digits=0))°")

if ellipticity < 1.1
    println("\n  → Nearly circular object (ellipticity < 1.1)")
elseif ellipticity < 1.3
    println("\n  → Slightly elongated (ellipticity 1.1-1.3)")
elseif ellipticity < 2.0
    println("\n  → Moderately elliptical (ellipticity 1.3-2.0)")
else
    println("\n  → Highly elongated (ellipticity > 2.0)")
end

println("\n" * "="^70)
println("Press Enter to exit...")
readline()
