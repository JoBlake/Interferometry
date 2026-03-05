# Check UV coverage and determine appropriate FOV
using OITOOLS
using FITSIO
using OIFITS
using Plots
using Statistics  # For mean function

# Add compatibility workarounds
Core.eval(OIFITS, quote
    function load(filename::String)
        return read(OIDataSet, filename; hack_revn=1)
    end
    function load(file)
        return read(OIDataSet, file; hack_revn=1)
    end
    function select(ds::OIDataSet, extname::String)
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

function Base.getproperty(target::OIFITS.OI_TARGET, sym::Symbol)
    if sym in (:revn, :list, :extname)
        return getfield(target, sym)
    end
    list = getfield(target, :list)
    if sym == :target_id && isa(list, AbstractVector)
        return [item.target_id for item in list]
    elseif hasfield(OIFITS.OITargetEntry, sym) && isa(list, AbstractVector)
        return [getfield(item, sym) for item in list]
    else
        error("type OI_TARGET has no property $sym")
    end
end

# File path
filename = "C:/Users/johnb/OneDrive/Documents/Projects/Astronomy/MIRCX_L2.2025May20.V_CVn.MIRCX_IDL.GHS.AVG6m.oifits"

println("Loading OIFITS data...")
data = readoifits(filename)

if isa(data, AbstractArray)
    data = data[1]
end

println("\n" * "="^60)
println("UV Coverage Analysis")
println("="^60)

# Get UV coordinates (stored as spatial frequencies in OITOOLS: baseline/wavelength)
u = data.uv[:,1]  # u spatial frequency (cycles/radian)
v = data.uv[:,2]  # v spatial frequency (cycles/radian)
spatial_freq = sqrt.(u.^2 .+ v.^2)  # Spatial frequency magnitude
wavelength = mean(data.uv_lam)  # Mean wavelength in meters

# Convert spatial frequencies back to baseline lengths
baseline = spatial_freq .* wavelength  # Baseline in meters

println("\nWavelength: $(wavelength*1e6) μm")
println("\nBaseline statistics:")
println("  Minimum: $(minimum(baseline)) m")
println("  Maximum: $(maximum(baseline)) m")
println("  Mean: $(mean(baseline)) m")

println("\nSpatial frequency coverage:")
println("  Minimum: $(minimum(spatial_freq)) cycles/rad")
println("  Maximum: $(maximum(spatial_freq)) cycles/rad")

# Maximum resolution (finest scale you can see)
max_resolution_rad = 1.0 / maximum(spatial_freq)
max_resolution_mas = max_resolution_rad * 206265000  # Convert to mas

println("\n" * "="^60)
println("RECOMMENDED IMAGE PARAMETERS")
println("="^60)

println("\nMaximum angular resolution:")
println("  θ_max = λ / B_max = $(max_resolution_mas) mas")
println("  (This is the finest detail you can resolve)")

# Nyquist sampling: need 2 pixels per resolution element
recommended_pixsize = max_resolution_mas / 2.0
println("\nRecommended pixel size (Nyquist sampling):")
println("  pixsize = θ_max / 2 = $(round(recommended_pixsize, digits=4)) mas/pixel")

# Field of view from shortest baseline (largest scale)
largest_scale_rad = 1.0 / minimum(spatial_freq)
largest_scale_mas = largest_scale_rad * 206265000
println("\nLargest angular scale (from shortest baseline):")
println("  θ_large = λ / B_min = $(round(largest_scale_mas, digits=3)) mas")
println("  (FOV much larger than this is not well-constrained by the data)")

# Recommended FOV: 3-5 times the object size
disk_diameter = 4.0  # mas
println("\nFor a $(disk_diameter) mas object:")
println("  Recommended FOV: $(round(3*disk_diameter, digits=1))-$(round(5*disk_diameter, digits=1)) mas (3-5× object size)")

# Recommended image size
println("\n" * "="^60)
println("IMAGE SIZE RECOMMENDATIONS:")
println("="^60)
for npix in [128, 256, 512]
    fov_for_this_npix = npix * recommended_pixsize
    println("\nOption: $npix × $npix pixels")
    println("  Pixel size: $(round(recommended_pixsize, digits=4)) mas/pixel")
    println("  FOV: $(round(fov_for_this_npix, digits=2)) mas")

    # Evaluate
    if fov_for_this_npix < 2 * disk_diameter
        println("  ⚠ Too small - won't contain the object")
    elseif fov_for_this_npix < 3 * disk_diameter
        println("  ⚠ Marginal - very tight around object")
    elseif fov_for_this_npix >= 3*disk_diameter && fov_for_this_npix <= 6*disk_diameter
        println("  ✓✓ EXCELLENT - Good padding around object")
    elseif fov_for_this_npix > 6*disk_diameter && fov_for_this_npix <= largest_scale_mas
        println("  ✓ Good - May have some empty space")
    else
        println("  ⚠ Too large - Risk of periodic artifacts")
    end
end

println("\n" * "="^60)
println("UV COVERAGE PLOT")
println("="^60)

# Plot UV coverage
p1 = scatter(u, v,
            label="UV points",
            xlabel="u (m)",
            ylabel="v (m)",
            aspect_ratio=:equal,
            title="UV Coverage",
            markersize=2)
scatter!(p1, -u, -v, label="", markersize=2)  # Add symmetric points

savefig(p1, "uv_coverage.png")
println("Saved: uv_coverage.png")

# Plot baseline histogram
p2 = histogram(baseline,
              xlabel="Baseline length (m)",
              ylabel="Count",
              title="Baseline Distribution",
              legend=false,
              bins=30)
savefig(p2, "baseline_distribution.png")
println("Saved: baseline_distribution.png")

println("\nDone! Check the plots to understand your UV coverage.")
println("Press Enter to exit...")
readline()
