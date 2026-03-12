# Image reconstruction
using OITOOLS
using FITSIO
using OIFITS
using Plots
using Statistics

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

# Prompt user for OIFITS file
println("="^70)
println("OIFITS IMAGE RECONSTRUCTION")
println("="^70)
println("\nEnter the path to your OIFITS file:")
println("(Press Enter to use: MIRCX_L2.2025May20.V_CVn.MIRCX_IDL.GHS.AVG6m.oifits)")
print("> ")

user_input = readline()

if isempty(strip(user_input))
    # Use default
    filename = "MIRCX_L2.2025May20.V_CVn.MIRCX_IDL.GHS.AVG6m.oifits"
    println("Using default file: $filename")
else
    filename = strip(user_input)
end

# Check if file exists
if !isfile(filename)
    println("\n⚠ Error: File not found: $filename")
    println("\nAvailable .oifits files in current directory:")
    oifits_files = filter(f -> endswith(f, ".oifits"), readdir("."))
    if isempty(oifits_files)
        println("  (none found)")
    else
        for f in oifits_files
            println("  - $f")
        end
    end
    println("\nPress Enter to exit...")
    readline()
    exit(1)
end

println("\n" * "="^70)
println("Loading OIFITS data from: $filename")
println("="^70)

# Now use OITOOLS readoifits function
data = nothing
try
    global data = readoifits(filename)
    println("Data loaded successfully")
    println("Data type: ", typeof(data))
    println("Data size: ", size(data))
catch e
    if occursin("OI_ARRAY", string(e))
        println("\n⚠ ERROR: This OIFITS file is missing the OI_ARRAY table.")
        println("\nThe OI_ARRAY table contains telescope array information that OITOOLS")
        println("requires for image reconstruction. Without it, reconstruction cannot proceed.")
        println("\nPossible solutions:")
        println("  1. Find a version of this data file that includes OI_ARRAY")
        println("  2. Use quick_diagnostic.jl to analyze the data quality")
        println("  3. Contact the data provider about the missing table")
        println("\nNote: Some older OIFITS v1 files may not have this table.")
        println("\nPress Enter to exit...")
        readline()
        exit(1)
    else
        println("\n⚠ ERROR: Failed to load OIFITS file:")
        println(e)
        println("\nPress Enter to exit...")
        readline()
        exit(1)
    end
end

# readoifits returns an array of OIdata objects
# For single-wavelength data, it's typically a 1-element array
if isa(data, AbstractArray)
    if length(data) == 1
        global data = data[1]  # Extract the single OIdata object
    else
        # For multi-wavelength, we might need the first channel or all of them
        println("Multi-wavelength data detected: $(length(data)) channels")
        global data = data[1]  # Use first wavelength channel for now
    end
end

println("Number of V² measurements: ", data.nv2)
println("Number of T3 measurements: ", data.nt3amp)

# Set up image reconstruction parameters
println("\nSetting up image reconstruction...")

# Image size and pixel scale
npix = 256  # Number of pixels (increased for better sampling)
fov = 100.0  # Large FOV for extended sources (~5-6x typical source size)
pixsize = fov / npix  # Pixel size in mas (~0.39 mas/pixel)

println("Image size: $npix x $npix pixels")
println("Pixel scale: $pixsize mas/pixel")
println("Field of view: $fov mas")

# Create initial image as a uniform disk model
# Disk diameter in mas
disk_diameter = 4.0  # mas
disk_radius_pixels = (disk_diameter / 2.0) / pixsize  # Convert to pixels

# Choose initial model type
initial_model = "gaussian"  # Use Gaussian - smooth starting point

# Create coordinate grids (used by all models)
center = npix / 2 + 0.5
x = repeat(reshape(1:npix, 1, npix), npix, 1) .- center  # Row vector repeated down
y = repeat(reshape(1:npix, npix, 1), 1, npix) .- center  # Column vector repeated across
r = sqrt.(x.^2 .+ y.^2)  # Radial distance from center in pixels

if initial_model == "disk"
    println("Creating initial uniform disk model...")

    # AUTO-ADJUST disk size based on mean V²
    # These thresholds match quick_diagnostic.jl recommendations
    mean_v2 = mean(data.v2)

    if mean_v2 > 0.7
        disk_diameter = 2.0  # Very compact, ~1-3 mas
    elseif mean_v2 > 0.4
        disk_diameter = 4.0  # Compact, ~3-6 mas
    elseif mean_v2 > 0.2
        disk_diameter = 8.0  # Moderately resolved, ~6-12 mas
    elseif mean_v2 > 0.1
        disk_diameter = 18.0  # Well resolved, ~12-25 mas
    elseif mean_v2 > 0.05
        disk_diameter = 25.0  # Very extended
    else
        disk_diameter = 30.0  # Extremely extended, >25 mas
    end

    println("  Data mean V²: $(round(mean_v2, digits=4))")
    println("  Auto-selected disk diameter: $disk_diameter mas")

    disk_radius_pixels = (disk_diameter / 2.0) / pixsize
    println("  Disk radius: $(round(disk_radius_pixels, digits=1)) pixels")

    # Create uniform disk: 1 inside radius, 0 outside
    initial_image = zeros(Float64, npix, npix)
    initial_image[r .<= disk_radius_pixels] .= 1.0

elseif initial_model == "gaussian"
    println("Creating initial Gaussian model...")

    # AUTO-ADJUST Gaussian size based on mean V²
    # These thresholds match quick_diagnostic.jl recommendations
    mean_v2 = mean(data.v2)

    if mean_v2 > 0.7
        fwhm_mas = 2.0  # Very compact, ~1-3 mas
    elseif mean_v2 > 0.4
        fwhm_mas = 4.0  # Compact, ~3-6 mas
    elseif mean_v2 > 0.2
        fwhm_mas = 8.0  # Moderately resolved, ~6-12 mas
    elseif mean_v2 > 0.1
        fwhm_mas = 18.0  # Well resolved, ~12-25 mas
    elseif mean_v2 > 0.05
        fwhm_mas = 25.0  # Very extended
    else
        fwhm_mas = 30.0  # Extremely extended, >25 mas
    end

    println("  Data mean V²: $(round(mean_v2, digits=4))")
    println("  Auto-selected FWHM: $fwhm_mas mas")

    fwhm_pixels = fwhm_mas / pixsize
    sigma_pixels = fwhm_pixels / 2.355  # Convert FWHM to sigma
    println("  Sigma: $(round(sigma_pixels, digits=1)) pixels")

    # Create 2D Gaussian
    initial_image = exp.(-(r.^2) ./ (2 * sigma_pixels^2))

else  # "flat"
    println("Creating flat uniform initial model...")
    # Completely uniform - no structure
    initial_image = ones(Float64, npix, npix)
end

# Normalize to unit flux
initial_image = initial_image ./ sum(initial_image)
println("  Peak pixel value: $(round(maximum(initial_image), sigdigits=4))")

# Save initial model and check if it matches data
println("\nSaving initial model image...")
extent = fov / 2
coords = range(-extent, extent, length=npix)
p_initial = heatmap(coords, coords,
            initial_image',
            aspect_ratio=:equal,
            xlabel="RA offset (mas)",
            ylabel="Dec offset (mas)",
            title="Initial Model - $initial_model",
            color=:hot,
            clims=(0, maximum(initial_image)))
savefig(p_initial, "initial_model.png")
println("Saved: initial_model.png")

# Quick check: how well does initial model match data?
println("\nChecking initial model visibility amplitudes...")
nfft_test = setup_nfft(data, npix, pixsize)
cvis_initial = image_to_vis(initial_image, nfft_test[1])
v2_initial = vis_to_v2(cvis_initial, data.indx_v2)
println("  Data V² range: $(round(minimum(data.v2), digits=3)) to $(round(maximum(data.v2), digits=3))")
println("  Model V² range: $(round(minimum(v2_initial), digits=3)) to $(round(maximum(v2_initial), digits=3))")
println("  Data V² mean: $(round(mean(data.v2), digits=3))")
println("  Model V² mean: $(round(mean(v2_initial), digits=3))")

if abs(mean(v2_initial) - mean(data.v2)) < 0.05
    println("  ✓ Initial model V² amplitude matches data well")
elseif mean(v2_initial) > mean(data.v2) + 0.1
    println("  ⚠ Initial model too compact - try larger size or 'flat' initial model")
elseif mean(v2_initial) < mean(data.v2) - 0.1
    println("  ⚠ Initial model too extended - try smaller size")
else
    println("  → Reasonable starting point")
end

# Regularization parameters
# Available regularizers: "centering", "tv", "tvsq", "l1l2", "l1l2w", "l1hyp",
#                         "l2sq", "compactness", "radialvar", "entropy", "support"

# Moderate regularization - balance between data fidelity and stability
# Reduced from original but not too weak
regularizers = [
    ("centering", 1e-4),    # Keep centered
    ("entropy", 5e-5),      # Light smoothing for stability
    ("tv", 1e-5),           # Very light total variation for edges
]

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
        maxiter=5000,  # More iterations for better convergence
        verb=true,
        regularizers=regularizers,
        weights=[1.0, 1.0, 1.0],  # Equal weight to all data - V2, T3 amplitude, T3 phase
        ftol=(1e-5, 1e-7),  # Tighter convergence tolerances
        xtol=(1e-5, 1e-7),
        gtol=(1e-5, 1e-7)
    )
    
    println("\nReconstruction complete!")

    # Display final image statistics
    println("\nFinal image statistics:")
    println("  Min pixel value: $(minimum(reconstructed_image))")
    println("  Max pixel value: $(maximum(reconstructed_image))")
    println("  Total flux: $(sum(reconstructed_image))")
    println("  Image dimensions: $(size(reconstructed_image))")

    # Display and save results
    println("\nCreating visualization...")
    
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
    println("\nError during reconstruction:")
    println(e)
    println("\nStacktrace:")
    for (exc, bt) in Base.catch_stack()
        showerror(stdout, exc, bt)
        println()
    end
end

println("\nImage reconstruction complete!")
println("Press Enter to exit...")
readline()
