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
println("\\nSetting up image reconstruction...")

# Image size and pixel scale
npix = 128  # Number of pixels (power of 2 works best)
pixsize = 0.1  # Pixel size in milliarcseconds (mas)

# Create initial image (uniform disk as starting point)
initial_image = ones(npix, npix)
initial_image = initial_image ./ sum(initial_image)  # Normalize

# Regularization parameters
regularizers = [
    ("entropy", 1e-3),     # Maximum entropy regularization
    ("tv", 1e-4),          # Total variation (edge preservation)
    ("l2", 1e-5)           # L2 smoothness
]

println("Image size: $npix x $npix pixels")
println("Pixel scale: $pixsize mas/pixel")
println("Field of view: $(npix * pixsize) mas")

# Run image reconstruction
println("\\nRunning image reconstruction...")
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
    reconstructed_image = reconstruct(
        initial_image,
        data,
        nfft_plan,
        maxiter=400,
        verb=true,
        regularizers=[("centering", 1e-3)]
    )
    
    println("\\nReconstruction complete!")
    
    # Display and save results
    println("Creating visualization...")
    
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
