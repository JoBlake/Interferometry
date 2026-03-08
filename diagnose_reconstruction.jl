# Diagnostic script to check NFFT setup and data quality
using OITOOLS
using FITSIO
using OIFITS
using Plots
using Statistics

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

println("="^70)
println("INTERFEROMETRIC DATA RECONSTRUCTION DIAGNOSTICS")
println("="^70)

# Prompt for filename
println("\nEnter the path to your OIFITS file:")
println("(Press Enter to use: MIRCX_L2.2025May20.V_CVn.MIRCX_IDL.GHS.AVG6m.oifits)")
print("> ")

user_input = readline()

if isempty(strip(user_input))
    filename = "MIRCX_L2.2025May20.V_CVn.MIRCX_IDL.GHS.AVG6m.oifits"
    println("Using default file: $filename")
else
    filename = strip(user_input)
end

if !isfile(filename)
    println("\n⚠ Error: File not found: $filename")
    println("\nPress Enter to exit...")
    readline()
    exit(1)
end

# Load data
println("\n1. Loading OIFITS data from: $filename...")
data = nothing
try
    data = readoifits(filename)
    if isa(data, AbstractArray)
        data = data[1]
    end
    println("   ✓ Data loaded successfully")
    println("   - V² measurements: $(data.nv2)")
    println("   - T3 amplitude measurements: $(data.nt3amp)")
    println("   - T3 phase measurements: $(data.nt3phi)")
    println("   - Total: $(data.nv2 + data.nt3amp + data.nt3phi) measurements")
catch e
    if occursin("OI_ARRAY", string(e))
        println("\n⚠ ERROR: This OIFITS file is missing the OI_ARRAY table.")
        println("\nThe OI_ARRAY table contains telescope array information that OITOOLS")
        println("requires for diagnostics. Without it, this script cannot proceed.")
        println("\nPossible solutions:")
        println("  1. Use quick_diagnostic.jl instead (supports files without OI_ARRAY)")
        println("  2. Find a version of this data file that includes OI_ARRAY")
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

# Check data quality
println("\n2. Checking data quality...")

v2_snr = data.v2 ./ data.v2_err
t3amp_snr = data.t3amp ./ data.t3amp_err
t3phi_snr = abs.(data.t3phi) ./ data.t3phi_err

println("   Signal-to-Noise Ratio statistics:")
println("   - V² SNR: min=$(round(minimum(v2_snr), digits=1)), max=$(round(maximum(v2_snr), digits=1)), median=$(round(median(v2_snr), digits=1))")
println("   - T3amp SNR: min=$(round(minimum(t3amp_snr), digits=1)), max=$(round(maximum(t3amp_snr), digits=1)), median=$(round(median(t3amp_snr), digits=1))")
println("   - T3phi SNR: min=$(round(minimum(t3phi_snr), digits=1)), max=$(round(maximum(t3phi_snr), digits=1)), median=$(round(median(t3phi_snr), digits=1))")

if median(v2_snr) < 3 || median(t3amp_snr) < 3
    println("   ⚠ WARNING: Low SNR data - reconstruction may be difficult")
else
    println("   ✓ SNR looks good")
end

# Check V2 values
println("\n   V² value statistics:")
println("   - Min: $(round(minimum(data.v2), digits=4))")
println("   - Max: $(round(maximum(data.v2), digits=4))")
println("   - Mean: $(round(mean(data.v2), digits=4))")

if maximum(data.v2) > 0.5
    println("   ✓ Good visibility amplitude - source is partially resolved")
elseif maximum(data.v2) < 0.1
    println("   ⚠ WARNING: Very low V² - source may be over-resolved")
else
    println("   ✓ Source is well-resolved")
end

# Check closure phases
println("\n   Closure phase statistics:")
println("   - RMS T3phi: $(round(std(data.t3phi), digits=2))°")
println("   - |T3phi| > 5°: $(sum(abs.(data.t3phi) .> 5)) measurements")

if std(data.t3phi) < 2
    println("   ⚠ NOTE: Small closure phases - source may be nearly centrosymmetric")
else
    println("   ✓ Significant closure phases detected - source has asymmetries")
end

# Test NFFT setup
println("\n3. Testing NFFT setup...")
npix = 256
fov = 12.0
pixsize = fov / npix

nfft_plan = setup_nfft(data, npix, pixsize)
println("   ✓ NFFT plan created successfully")
println("   - Image size: $npix × $npix")
println("   - Pixel scale: $(round(pixsize, digits=4)) mas/pixel")
println("   - FOV: $fov mas")

# Test forward model with uniform disk
println("\n4. Testing forward model (uniform disk → visibilities)...")
disk_diameter = 4.0
disk_radius_pixels = (disk_diameter / 2.0) / pixsize

center = npix / 2 + 0.5
x = repeat(1:npix, 1, npix) .- center
y = repeat((1:npix)', npix, 1) .- center
r = sqrt.(x.^2 .+ y.^2)

test_image = zeros(Float64, npix, npix)
test_image[r .<= disk_radius_pixels] .= 1.0
test_image = test_image ./ sum(test_image)

# Compute model visibilities
cvis_model = image_to_vis(test_image, nfft_plan[1])
v2_model = vis_to_v2(cvis_model, data.indx_v2)
t3_model, t3amp_model, t3phi_model = vis_to_t3(cvis_model, data.indx_t3_1, data.indx_t3_2, data.indx_t3_3)

println("   Model V² range: $(round(minimum(v2_model), digits=4)) to $(round(maximum(v2_model), digits=4))")
println("   Data V² range: $(round(minimum(data.v2), digits=4)) to $(round(maximum(data.v2), digits=4))")

# Compute chi-squared for uniform disk model
chi2_v2_disk = sum(((v2_model - data.v2) ./ data.v2_err).^2) / data.nv2
chi2_t3amp_disk = sum(((t3amp_model - data.t3amp) ./ data.t3amp_err).^2) / data.nt3amp
chi2_t3phi_disk = sum(((t3phi_model - data.t3phi) ./ data.t3phi_err).^2) / data.nt3phi

println("\n   Uniform disk ($disk_diameter mas) chi-squared:")
println("   - V²: $(round(chi2_v2_disk, digits=2))")
println("   - T3amp: $(round(chi2_t3amp_disk, digits=2))")
println("   - T3phi: $(round(chi2_t3phi_disk, digits=2))")

if chi2_v2_disk < 2 && chi2_t3amp_disk < 2 && chi2_t3phi_disk < 2
    println("   ⚠ WARNING: Uniform disk fits data very well!")
    println("      → Data may truly support a simple uniform disk")
    println("      → Or: Source size matches your model closely")
elseif chi2_v2_disk > 100
    println("   ⚠ WARNING: Uniform disk fits data very poorly!")
    println("      → Check if disk size ($disk_diameter mas) is appropriate")
else
    println("   ✓ Uniform disk is a poor fit - reconstruction should improve this")
end

# Test with point source
println("\n5. Testing with point source model...")
point_image = zeros(Float64, npix, npix)
point_image[div(npix,2)+1, div(npix,2)+1] = 1.0

cvis_point = image_to_vis(point_image, nfft_plan[1])
v2_point = vis_to_v2(cvis_point, data.indx_v2)

chi2_v2_point = sum(((v2_point - data.v2) ./ data.v2_err).^2) / data.nv2

println("   Point source chi-squared V²: $(round(chi2_v2_point, digits=2))")

if chi2_v2_point < chi2_v2_disk
    println("   ⚠ WARNING: Point source fits better than disk!")
    println("      → Source may be much smaller than $disk_diameter mas")
else
    println("   ✓ Extended source confirmed")
end

# Check UV coverage uniformity
println("\n6. Checking UV coverage uniformity...")
u = data.uv[:,1]
v = data.uv[:,2]
spatial_freq = sqrt.(u.^2 .+ v.^2)

# Divide UV plane into annuli
freq_bins = range(minimum(spatial_freq), maximum(spatial_freq), length=10)
bin_counts = [sum((spatial_freq .>= freq_bins[i]) .& (spatial_freq .< freq_bins[i+1])) for i in 1:length(freq_bins)-1]

println("   UV coverage distribution:")
for i in 1:length(bin_counts)
    freq_range = "$(round(freq_bins[i]/1e13, digits=1))-$(round(freq_bins[i+1]/1e13, digits=1))e13"
    bar = "▮"^(div(bin_counts[i], maximum(bin_counts)) * 30)
    println("   $freq_range: $bar ($(bin_counts[i]))")
end

if minimum(bin_counts) < 0.1 * maximum(bin_counts)
    println("   ⚠ WARNING: Uneven UV coverage - some spatial frequencies poorly sampled")
else
    println("   ✓ UV coverage is reasonably uniform")
end

# Summary and recommendations
println("\n" * "="^70)
println("SUMMARY AND RECOMMENDATIONS")
println("="^70)

issues = []
recommendations = []

if median(v2_snr) < 3
    push!(issues, "Low SNR data")
    push!(recommendations, "Consider increasing regularization weights")
end

if chi2_v2_disk < 2 && chi2_t3phi_disk < 2
    push!(issues, "Data fits simple uniform disk very well")
    push!(recommendations, "This may be the correct model - complex features may not be justified")
    push!(recommendations, "Try larger or smaller disk diameters (2-8 mas)")
end

if chi2_v2_disk > 100
    push!(issues, "Uniform disk is a very poor fit")
    push!(recommendations, "Try different disk sizes: $(round(disk_diameter*0.5, digits=1)), $(round(disk_diameter*2, digits=1)) mas")
    push!(recommendations, "Source may be significantly different from circular")
end

if std(data.t3phi) < 2
    push!(issues, "Very small closure phases")
    push!(recommendations, "Source may be nearly centrosymmetric - asymmetric features may not be present")
end

if length(issues) > 0
    println("\n⚠ Potential Issues Detected:")
    for (i, issue) in enumerate(issues)
        println("   $i. $issue")
    end
end

if length(recommendations) > 0
    println("\n💡 Recommendations:")
    for (i, rec) in enumerate(recommendations)
        println("   $i. $rec")
    end
else
    println("\n✓ No major issues detected - NFFT and data appear healthy")
    println("\nIf reconstruction still looks like initial model:")
    println("  1. Reduce regularization weights by another 10×")
    println("  2. Try flat uniform initial image instead of disk")
    println("  3. Increase maxiter to 3000-5000")
    println("  4. Check final iteration output for regularization vs chi² balance")
end

println("\n" * "="^70)
println("Press Enter to exit...")
readline()
