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
println("(Press Enter to use: Pi_GRU.oifits)")
print("> ")

user_input = readline()

if isempty(strip(user_input))
    # Use default
    filename = "Pi_GRU.oifits"
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
fov = 60.0  # FOV optimized for well-resolved ~15-20 mas stellar disks
pixsize = fov / npix  # Pixel size in mas (~0.23 mas/pixel)

println("Image size: $npix x $npix pixels")
println("Pixel scale: $pixsize mas/pixel")
println("Field of view: $fov mas")

# Prompt user for initial model configuration
println("\n" * "="^70)
println("INITIAL MODEL CONFIGURATION")
println("="^70)

# Prompt for initial model type
println("\nSelect initial model type:")
println("  1. Uniform disk (recommended for Pi_GRU)")
println("  2. Gaussian")
println("  3. Flattened Gaussian (disk-like)")
print("Enter choice [1-3, default: 1]: ")
model_choice = strip(readline())

if isempty(model_choice)
    initial_model = "disk"
    println("  Using default: Uniform disk")
elseif model_choice == "1"
    initial_model = "disk"
elseif model_choice == "2"
    initial_model = "gaussian"
elseif model_choice == "3"
    initial_model = "flattened_gaussian"
else
    println("⚠ Invalid choice. Using default: disk")
    initial_model = "disk"
end

# Auto-estimate size based on data
mean_v2 = mean(data.v2)
if mean_v2 > 0.7
    suggested_diameter = 2.0  # Very compact, ~1-3 mas
elseif mean_v2 > 0.4
    suggested_diameter = 4.0  # Compact, ~3-6 mas
elseif mean_v2 > 0.2
    suggested_diameter = 8.0  # Moderately resolved, ~6-12 mas
elseif mean_v2 > 0.1
    suggested_diameter = 18.0  # Well resolved, ~12-25 mas
elseif mean_v2 > 0.05
    suggested_diameter = 25.0  # Very extended
else
    suggested_diameter = 30.0  # Extremely extended, >25 mas
end

# Prompt for initial disk/source size
println("\nData analysis: Mean V² = $(round(mean_v2, digits=4))")
println("Suggested initial diameter: $suggested_diameter mas")
print("Enter initial diameter in mas [default: $suggested_diameter]: ")
diameter_input = strip(readline())

if isempty(diameter_input)
    disk_diameter = suggested_diameter
    println("  Using suggested: $disk_diameter mas")
else
    disk_diameter = tryparse(Float64, diameter_input)
    if disk_diameter === nothing || disk_diameter <= 0
        println("⚠ Invalid input. Using suggested: $suggested_diameter mas")
        disk_diameter = suggested_diameter
    end
end

disk_radius_pixels = (disk_diameter / 2.0) / pixsize  # Convert to pixels
println("  Diameter: $disk_diameter mas")
println("  Radius: $(round(disk_radius_pixels, digits=1)) pixels")

# Create coordinate grids (used by all models)
center = npix / 2 + 0.5
x = repeat(reshape(1:npix, 1, npix), npix, 1) .- center  # Row vector repeated down
y = repeat(reshape(1:npix, npix, 1), 1, npix) .- center  # Column vector repeated across
r = sqrt.(x.^2 .+ y.^2)  # Radial distance from center in pixels

if initial_model == "disk"
    println("\nCreating initial uniform disk model...")
    println("  Using user-specified diameter: $disk_diameter mas")

    # Create uniform disk: 1 inside radius, 0 outside
    initial_image = zeros(Float64, npix, npix)
    initial_image[r .<= disk_radius_pixels] .= 1.0

elseif initial_model == "gaussian"
    println("\nCreating initial Gaussian model...")
    println("  Using user-specified diameter: $disk_diameter mas")

    # Use disk_diameter as FWHM
    fwhm_mas = disk_diameter
    fwhm_pixels = fwhm_mas / pixsize
    sigma_pixels = fwhm_pixels / 2.355  # Convert FWHM to sigma
    println("  FWHM: $fwhm_mas mas")
    println("  Sigma: $(round(sigma_pixels, digits=1)) pixels")

    # Create pure 2D Gaussian
    initial_image = exp.(-(r.^2) ./ (2 * sigma_pixels^2))

elseif initial_model == "flattened_gaussian"
    println("\nCreating initial flattened Gaussian model...")
    println("  Using user-specified diameter: $disk_diameter mas")

    # Use disk_diameter as FWHM
    fwhm_mas = disk_diameter
    fwhm_pixels = fwhm_mas / pixsize
    sigma_pixels = fwhm_pixels / 2.355  # Convert FWHM to sigma
    println("  FWHM: $fwhm_mas mas")
    println("  Sigma: $(round(sigma_pixels, digits=1)) pixels")

    # Create 2D Gaussian
    initial_image = exp.(-(r.^2) ./ (2 * sigma_pixels^2))

    # Flatten the peak - clip values above threshold
    # Use 70% of peak to create a flat top in the center
    flatten_threshold = 0.7
    peak_value = maximum(initial_image)
    flatten_level = flatten_threshold * peak_value
    initial_image[initial_image .> flatten_level] .= flatten_level

    println("  Peak flattening: Values above $(round(flatten_threshold*100, digits=0))% of peak set to plateau")
    println("  This creates a more disk-like profile")

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
println("Initial model statistics:")
println("  Min: $(minimum(initial_image))")
println("  Max: $(maximum(initial_image))")
println("  Mean: $(mean(initial_image))")
println("  Sum: $(sum(initial_image))")

extent = fov / 2
coords = range(-extent, extent, length=npix)

# For initial model, use appropriate color scale
init_min = minimum(initial_image)
init_max = maximum(initial_image)
if init_max - init_min < 0.1 * init_max
    init_clims = (init_min, init_max)
else
    init_clims = (0, init_max)
end

p_initial = heatmap(coords, coords,
            initial_image',
            aspect_ratio=:equal,
            xlabel="RA offset (mas)",
            ylabel="Dec offset (mas)",
            title="Initial Model - $initial_model",
            color=:hot,
            clims=init_clims)
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

# Setup NFFT once (reused for all reconstructions)
println("\nSetting up NFFT...")
nfft_plan = setup_nfft(data, npix, pixsize)

# Interactive loop for regularization parameters
reconstruction_count = 0
while true
    global reconstruction_count
    reconstruction_count += 1

    println("\n" * "="^70)
    println("REGULARIZATION PARAMETER INPUT (Reconstruction #$reconstruction_count)")
    println("="^70)
    println("\nAvailable regularizers: centering, entropy, tv (Total Variation)")
    println("Typical ranges: 1e-10 (very weak) to 1e-5 (strong)")
    println("\nOptimized defaults for Pi_GRU (well-resolved, asymmetric source):")
    println("  - Centering: 1e-6 (minimal image drift)")
    println("  - Entropy:   1e-8 (very weak smoothing - let data dominate)")
    println("  - TV:        0 (disabled - can cause ring artifacts)")
    println("\nPress Enter for defaults, or enter custom value")
    println("Type 'exit' at any prompt to quit\n")

    # Input for centering parameter
    print("Enter CENTERING parameter [1e-6]: ")
    centering_input = strip(readline())
    if lowercase(centering_input) == "exit"
        println("\nExiting reconstruction loop...")
        break
    end
    if isempty(centering_input)
        centering_value = 1e-6
        println("  Using default: 1e-6")
    else
        centering_value = tryparse(Float64, centering_input)
        if centering_value === nothing
            println("⚠ Invalid input. Using default: 1e-6")
            centering_value = 1e-6
        end
    end

    # Input for entropy parameter
    print("Enter ENTROPY parameter [1e-8]: ")
    entropy_input = strip(readline())
    if lowercase(entropy_input) == "exit"
        println("\nExiting reconstruction loop...")
        break
    end
    if isempty(entropy_input)
        entropy_value = 1e-8
        println("  Using default: 1e-8")
    else
        entropy_value = tryparse(Float64, entropy_input)
        if entropy_value === nothing
            println("⚠ Invalid input. Using default: 1e-8")
            entropy_value = 1e-8
        end
    end

    # Input for TV parameter
    print("Enter TV (Total Variation) parameter [0]: ")
    tv_input = strip(readline())
    if lowercase(tv_input) == "exit"
        println("\nExiting reconstruction loop...")
        break
    end
    if isempty(tv_input)
        tv_value = 0.0
        println("  Using default: 0 (TV disabled)")
    else
        tv_value = tryparse(Float64, tv_input)
        if tv_value === nothing
            println("⚠ Invalid input. Using default: 0")
            tv_value = 0.0
        end
    end

    # Build regularizers list (only include non-zero values)
    regularizers = []
    if centering_value > 0
        push!(regularizers, ("centering", centering_value))
    end
    if entropy_value > 0
        push!(regularizers, ("entropy", entropy_value))
    end
    if tv_value > 0
        push!(regularizers, ("tv", tv_value))
    end

    println("\n" * "-"^70)
    println("Regularization settings:")
    println("  Centering: $centering_value")
    println("  Entropy:   $entropy_value")
    tv_status = tv_value == 0 ? " (disabled)" : ""
    println("  TV:        $tv_value$tv_status")
    println("  Active regularizers: $(length(regularizers))")
    println("-"^70)

    # Run image reconstruction
    println("\nRunning image reconstruction...")
    println("Configuration: DATA-DRIVEN mode with custom regularization")
    println("  - Uniform disk initial model")
    println("  - Minimal regularization (data dominates)")
    println("  - Extended iterations for full convergence")
    println("This may take several minutes...")

    try
        # Basic reconstruction using OITOOLS
        println("Starting image reconstruction...")

        # reconstruct expects: (initial_image, data, nfft_plan; keyword_args)
        println("Running reconstruction algorithm...")
        println("Note: Look for colored output showing V2, T3A, T3P values on each iteration")
        println("      These are REDUCED chi-squared values (chi²/N_measurements)")
        println("      Target: values close to 1.0 indicate good fit\n")

        reconstructed_image = reconstruct(
            initial_image,
            data,
            nfft_plan,
            maxiter=5000,  # Sufficient iterations for well-sampled data
            verb=true,
            regularizers=regularizers,
            weights=[1.0, 1.0, 1.0],  # Equal weight to all data - V2, T3 amplitude, T3 phase
            ftol=(1e-5, 1e-7),  # Balanced convergence criteria
            xtol=(1e-5, 1e-7),
            gtol=(1e-5, 1e-7)
        )

        println("\nReconstruction complete!")

        # Display final image statistics
        println("\nFinal image statistics:")
        img_min = minimum(reconstructed_image)
        img_max = maximum(reconstructed_image)
        img_mean = mean(reconstructed_image)
        img_std = std(reconstructed_image)
        dynamic_range = img_max - img_min

        println("  Min pixel value: $img_min")
        println("  Max pixel value: $img_max")
        println("  Mean pixel value: $img_mean")
        println("  Std dev: $img_std")
        println("  Dynamic range: $dynamic_range")
        println("  Total flux: $(sum(reconstructed_image))")
        println("  Image dimensions: $(size(reconstructed_image))")

        # Check if image is too uniform
        relative_range = dynamic_range / img_max
        if relative_range < 0.01
            println("\n⚠ WARNING: Image is very uniform (dynamic range < 1%)")
            println("  This suggests regularization is too strong.")
            println("  Try weaker parameters (e.g., 1e-8 to 1e-9)")
        elseif relative_range < 0.1
            println("\n⚠ Note: Image has low contrast (dynamic range < 10%)")
            println("  Consider reducing regularization for more structure.")
        end

        # Display and save results
        println("\nCreating visualization...")

        # Create coordinate arrays for plotting
        extent = npix * pixsize / 2
        coords = range(-extent, extent, length=npix)

        # Choose appropriate color limits
        # For very uniform images, use min-max range to show any variation
        # For normal images, use 0-max to show structure clearly
        if relative_range < 0.1
            color_limits = (img_min, img_max)
            println("  Using full dynamic range for color scale (min to max)")
        else
            color_limits = (0, img_max)
            println("  Using standard color scale (0 to max)")
        end

        # Plot the reconstructed image with parameter info in title
        param_str = "C=$(centering_value) E=$(entropy_value) TV=$(tv_value)"
        p = heatmap(coords, coords,
                    reconstructed_image',
                    aspect_ratio=:equal,
                    xlabel="RA offset (mas)",
                    ylabel="Dec offset (mas)",
                    title="Reconstructed Image #$reconstruction_count\n$param_str",
                    color=:hot,
                    clims=color_limits)

        # Save with unique filenames
        png_filename = "reconstructed_image_$(reconstruction_count).png"
        fits_filename = "reconstructed_image_$(reconstruction_count).fits"

        savefig(p, png_filename)
        println("Saved: $png_filename")

        # Also save as FITS file
        println("Saving FITS file...")
        save_fits(reconstructed_image, fits_filename, pixsize)

        println("\n✓ Reconstruction #$reconstruction_count complete!")
        println("  Files saved: $png_filename, $fits_filename")

    catch e
        println("\n⚠ Error during reconstruction:")
        println(e)
        println("\nStacktrace:")
        for (exc, bt) in Base.catch_stack()
            showerror(stdout, exc, bt)
            println()
        end
        println("\nPress Enter to continue or type 'exit' to quit...")
        continue_input = strip(readline())
        if lowercase(continue_input) == "exit"
            break
        end
    end
end

println("\n" * "="^70)
println("Image reconstruction session complete!")
println("Total reconstructions performed: $reconstruction_count")
println("="^70)
