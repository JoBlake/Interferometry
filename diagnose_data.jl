using OITOOLS
using OIFITS
using Statistics

# Prompt user for OIFITS file
println("="^70)
println("OIFITS DATA DIAGNOSTICS")
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

println("\nLoading $filename...")

# Use the same workaround as reconstruct.jl
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

data = readoifits(filename)
if isa(data, AbstractArray) && length(data) == 1
    data = data[1]
elseif isa(data, AbstractArray)
    println("Multi-wavelength data: $(length(data)) channels, using first")
    data = data[1]
end

println("\n" * "="^70)
println("DATA DIAGNOSTICS - $filename")
println("="^70)

# V² statistics
println("\nVisibility² (V²) Data:")
println("  Number of measurements: $(data.nv2)")
println("  Mean V²: $(round(mean(data.v2), digits=4))")
println("  Min V²: $(round(minimum(data.v2), digits=4))")
println("  Max V²: $(round(maximum(data.v2), digits=4))")
println("  Std V²: $(round(std(data.v2), digits=4))")

# T3 amplitude statistics
if data.nt3amp > 0
    println("\nClosure Amplitude (T3) Data:")
    println("  Number of measurements: $(data.nt3amp)")
    println("  Mean T3 amp: $(round(mean(data.t3amp), digits=4))")
    println("  Min T3 amp: $(round(minimum(data.t3amp), digits=4))")
    println("  Max T3 amp: $(round(maximum(data.t3amp), digits=4))")
end

# T3 phase statistics
if data.nt3phi > 0
    println("\nClosure Phase (T3φ) Data:")
    println("  Number of measurements: $(data.nt3phi)")
    println("  Mean |T3φ|: $(round(mean(abs.(data.t3phi)), digits=2))°")
    println("  Max |T3φ|: $(round(maximum(abs.(data.t3phi)), digits=2))°")
    println("  Std T3φ: $(round(std(data.t3phi), digits=2))°")
end

# UV coverage
println("\nUV Coverage:")
println("  U range: $(round(minimum(data.uv[1,:]), sigdigits=3)) to $(round(maximum(data.uv[1,:]), sigdigits=3)) cycles/rad")
println("  V range: $(round(minimum(data.uv[2,:]), sigdigits=3)) to $(round(maximum(data.uv[2,:]), sigdigits=3)) cycles/rad")

# Estimate source size from V²
mean_v2 = mean(data.v2)
println("\n" * "="^70)
println("SOURCE SIZE ESTIMATE")
println("="^70)
if mean_v2 > 0.7
    println("  Mean V² = $mean_v2 → Very compact source (~1-3 mas)")
    println("  Recommended FOV: 10-20 mas")
elseif mean_v2 > 0.4
    println("  Mean V² = $mean_v2 → Compact source (~3-6 mas)")
    println("  Recommended FOV: 20-30 mas")
elseif mean_v2 > 0.2
    println("  Mean V² = $mean_v2 → Moderately resolved (~6-12 mas)")
    println("  Recommended FOV: 30-50 mas")
else
    println("  Mean V² = $mean_v2 → Well resolved (>12 mas)")
    println("  Recommended FOV: 50-100 mas")
end

println("\n" * "="^70)
