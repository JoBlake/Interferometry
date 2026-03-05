# Image reconstruction
using OITOOLS
using FITSIO
using OIFITS
using Plots

# Workaround: Define missing OIFITS functions that OITOOLS expects
Core.eval(OIFITS, quote
    function load(filename::String)
        # Try reading with hack_revn option for older OIFITS files
        return read(OIDataSet, filename; hack_revn=1)
    end

    function load(file)
        return read(OIDataSet, file; hack_revn=1)
    end

    function select(ds::OIDataSet, extname::String)
        # Select specific table types from the dataset
        if extname == "OI_WAVELENGTH"
            return ds.instr
        elseif extname == "OI_TARGET"
            return [ds.target]
        elseif extname == "OI_VIS"
            return ds.vis
        elseif extname == "OI_VIS2"
            return ds.vis2
        elseif extname == "OI_T3"
            return ds.t3
        elseif extname == "OI_FLUX"
            return ds.flux
        elseif extname == "OI_ARRAY"
            return ds.array
        elseif extname == "OI_CORREL"
            return ds.correl
        elseif extname == "OI_INSPOL"
            return ds.inspol
        else
            return []
        end
    end
end)

# Add compatibility properties for OI_TARGET
# In OIFITS 2.0, OI_TARGET.list contains OITargetEntry objects
# OITOOLS expects to access properties like target_id directly on OI_TARGET

function Base.getproperty(target::OIFITS.OI_TARGET, sym::Symbol)
    # First check if it's a standard field
    if sym in (:revn, :list, :extname)
        return getfield(target, sym)
    end

    # For other properties, extract from the list field
    list = getfield(target, :list)
    if sym == :target_id && isa(list, AbstractVector)
        return [item.target_id for item in list]
    elseif hasfield(OIFITS.OITargetEntry, sym) && isa(list, AbstractVector)
        # Generic extraction for any OITargetEntry field
        return [getfield(item, sym) for item in list]
    else
        # Fallback to error
        error("type OI_TARGET has no property $sym")
    end
end

# Helper function to save FITS - defined early so it can be used later
function save_fits(image, filename, pixsize)
    f = FITS(filename, "w")

    # Create header
    header = FITSHeader(["SIMPLE", "BITPIX", "NAXIS", "NAXIS1", "NAXIS2",
                        "CDELT1", "CDELT2", "CRPIX1", "CRPIX2",
                        "CTYPE1", "CTYPE2", "CUNIT1", "CUNIT2"],
                       [true, -32, 2, size(image, 1), size(image, 2),
                        pixsize/3600000, pixsize/3600000,
                        size(image, 1)/2, size(image, 2)/2,
                        "RA---TAN", "DEC--TAN", "deg", "deg"],
                       ["", "32-bit floating point", "2D image",
                        "Width", "Height", "Pixel scale RA", "Pixel scale Dec",
                        "Reference pixel X", "Reference pixel Y",
                        "Coordinate type", "Coordinate type",
                        "Units", "Units"])

    write(f, Float32.(image), header=header)
    close(f)

    println("FITS file saved: $filename")
end

# File path
filename = "C:/Users/johnb/OneDrive/Documents/Projects/Astronomy/MIRCX_L2.2025May20.V_CVn.MIRCX_IDL.GHS.AVG6m.oifits"

println("Loading OIFITS data...")

# Now use OITOOLS readoifits function
data = readoifits(filename)

println("Data loaded successfully")
println("Data type: ", typeof(data))
println("Data size: ", size(data))

# readoifits returns an array of OIdata objects
# For single-wavelength data, it's typically a 1-element array
if isa(data, AbstractArray)
    if length(data) == 1
        data = data[1]  # Extract the single OIdata object
    else
        # For multi-wavelength, we might need the first channel or all of them
        println("Multi-wavelength data detected: $(length(data)) channels")
        data = data[1]  # Use first wavelength channel for now
    end
end

println("Number of V² measurements: ", data.nv2)
println("Number of T3 measurements: ", data.nt3amp)

# Set up image reconstruction parameters
println("\nSetting up image reconstruction...")

# Image size and pixel scale
npix = 256  # Number of pixels (power of 2 works best)
fov = 12.0  # Field of view in milliarcseconds (mas) - reduced to avoid periodic artifacts
pixsize = fov / npix  # Pixel size in mas (~0.047 mas/pixel)

println("Image size: $npix x $npix pixels")
println("Pixel scale: $pixsize mas/pixel")
println("Field of view: $fov mas")

# Create initial image as a uniform disk model
# Disk diameter in mas
disk_diameter = 4.0  # mas
disk_radius_pixels = (disk_diameter / 2.0) / pixsize  # Convert to pixels

println("Creating initial uniform disk model...")
println("  Disk diameter: $disk_diameter mas")
println("  Disk radius: $disk_radius_pixels pixels")

# Create coordinate grids
center = npix / 2 + 0.5
x = repeat(1:npix, 1, npix) .- center
y = repeat((1:npix)', npix, 1) .- center
r = sqrt.(x.^2 .+ y.^2)  # Radial distance from center in pixels

# Create uniform disk: 1 inside radius, 0 outside
initial_image = zeros(Float64, npix, npix)
initial_image[r .<= disk_radius_pixels] .= 1.0

# Normalize to unit flux
initial_image = initial_image ./ sum(initial_image)

# Optional: Save initial model for comparison
println("Saving initial model image...")
extent = fov / 2
coords = range(-extent, extent, length=npix)
p_initial = heatmap(coords, coords,
            initial_image',
            aspect_ratio=:equal,
            xlabel="RA offset (mas)",
            ylabel="Dec offset (mas)",
            title="Initial Model - Uniform Disk ($disk_diameter mas)",
            color=:hot,
            clims=(0, maximum(initial_image)))
savefig(p_initial, "initial_model.png")
println("Saved: initial_model.png")

# Regularization parameters
# Available regularizers: "centering", "tv", "tvsq", "l1l2", "l1l2w", "l1hyp",
#                         "l2sq", "compactness", "radialvar", "entropy", "support"

# Strategy 1: Current (balanced smoothness + edges)
regularizers = [
    ("centering", 1e-2),   # Strong centering to avoid edge artifacts
    ("entropy", 1e-3),     # Maximum entropy regularization
    ("tv", 1e-4),          # Total variation (edge preservation)
    ("l2sq", 1e-5)         # L2 squared smoothness
]

# Strategy 2: Edge-preserving (uncomment to try)
# regularizers = [
#     ("centering", 1e-2),
#     ("tv", 1e-3),        # Dominant: preserves sharp features
#     ("l2sq", 1e-5)
# ]

# Strategy 3: Maximum smoothness (uncomment to try)
# regularizers = [
#     ("centering", 1e-2),
#     ("entropy", 1e-2),   # Dominant: very smooth result
# ]

# Strategy 4: Compactness (for point-like sources, uncomment to try)
# regularizers = [
#     ("compactness", 1e-2),  # Keeps object compact and centered
#     ("entropy", 1e-3),
# ]

# Run image reconstruction
println("\nRunning image reconstruction...")
println("This may take several minutes...")

try
    # Basic reconstruction using OITOOLS
    println("Starting image reconstruction with maximum entropy regularization...")

    # Setup NFFT (Non-uniform FFT) for the reconstruction
    # This creates the mapping between image and visibilities
    println("Setting up NFFT...")
    nfft_plan = setup_nfft(data, npix, pixsize)

    # reconstruct expects: (initial_image, data, nfft_plan; keyword_args)
    println("Running reconstruction algorithm...")
    println("Note: Look for colored output showing V2, T3A, T3P values on each iteration")
    println("      These are REDUCED chi-squared values (chi²/N_measurements)")
    println("      Target: values close to 1.0 indicate good fit\n")

    reconstructed_image = reconstruct(
        initial_image,
        data,
        nfft_plan,
        maxiter=2000,        # Increased, but should converge before this
        verb=true,
        regularizers=regularizers,
        ftol=(1e-4, 1e-6),   # Relaxed relative tolerance: stop when chi² changes < 0.01%
        xtol=(1e-4, 1e-6),   # Relaxed relative tolerance for image changes
        gtol=(1e-4, 1e-6)    # Relaxed relative tolerance for gradient
    )
    
    println("\\nReconstruction complete!")

    # Display final image statistics
    println("\\nFinal image statistics:")
    println("  Min pixel value: $(minimum(reconstructed_image))")
    println("  Max pixel value: $(maximum(reconstructed_image))")
    println("  Total flux: $(sum(reconstructed_image))")
    println("  Image dimensions: $(size(reconstructed_image))")

    # Display and save results
    println("\\nCreating visualization...")
    
    # Create coordinate arrays for plotting
    extent = npix * pixsize / 2
    coords = range(-extent, extent, length=npix)
    
    # Plot the reconstructed image
    p = heatmap(coords, coords,
                reconstructed_image',
                aspect_ratio=:equal,
                xlabel="RA offset (mas)",
                ylabel="Dec offset (mas)",
                title="Reconstructed Image - V CVn",
                color=:hot,
                clims=(0, maximum(reconstructed_image)))
    
    savefig(p, "reconstructed_image.png")
    println("Saved: reconstructed_image.png")
    
    # Also save as FITS file
    println("Saving FITS file...")
    save_fits(reconstructed_image, "reconstructed_image.fits", pixsize)
    
catch e
    println("\\nError during reconstruction:")
    println(e)
    println("\\nStacktrace:")
    for (exc, bt) in Base.catch_stack()
        showerror(stdout, exc, bt)
        println()
    end
end

println("\nImage reconstruction complete!")
println("Press Enter to exit...")
readline()
