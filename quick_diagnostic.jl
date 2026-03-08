# Quick diagnostic to understand the data
using OITOOLS
using FITSIO
using OIFITS
using Statistics

Core.eval(OIFITS, quote
    function load(filename::String)
        return read(OIDataSet, filename; hack_revn=1)
    end
    function select(ds::OIDataSet, extname::String)
        if extname == "OI_WAVELENGTH"
            return ds.instr
        elseif extname == "OI_TARGET"
            return [ds.target]
        elseif extname == "OI_VIS2"
            return ds.vis2
        elseif extname == "OI_T3"
            return ds.t3
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

filename = "2009-11-eps_Aur-avg5.oifits"
println("Quick Data Check for: $filename\n")

# Try reading with OITOOLS, fallback to direct OIFITS reading
data = nothing
try
    data = readoifits(filename)
    if isa(data, AbstractArray)
        data = data[1]
    end
catch e
    println("⚠ OITOOLS readoifits failed (missing OI_ARRAY table)")
    println("  Using direct OIFITS read instead...\n")

    # Read with OIFITS package directly
    oifits_data = OIFITS.load(filename)

    # Extract V² and T3 data manually
    vis2_tables = OIFITS.select(oifits_data, "OI_VIS2")
    t3_tables = OIFITS.select(oifits_data, "OI_T3")

    if length(vis2_tables) > 0
        # Get V² data from first table
        v2_data = vis2_tables[1].vis2data
        v2_err = vis2_tables[1].vis2err

        println("V² Statistics:")
        println("  Min: $(round(minimum(v2_data), digits=4))")
        println("  Max: $(round(maximum(v2_data), digits=4))")
        println("  Mean: $(round(mean(v2_data), digits=4))")
        println("  Median: $(round(median(v2_data), digits=4))")

        mean_v2 = mean(v2_data)
    else
        println("⚠ No V² data found in file")
        println("Press Enter to exit...")
        readline()
        exit(1)
    end

    if length(t3_tables) > 0
        t3phi_data = t3_tables[1].t3phi
        println("\nClosure phase RMS: $(round(std(t3phi_data), digits=1))°")
    end

    # Provide recommendation
    println("\nWhat this means:")
    if mean_v2 > 0.7
        println("  → Very compact source, diameter ~1-3 mas")
        suggested_size = 2.0
    elseif mean_v2 > 0.4
        println("  → Compact source, diameter ~3-6 mas")
        suggested_size = 4.0
    elseif mean_v2 > 0.2
        println("  → Moderately resolved, diameter ~6-12 mas")
        suggested_size = 8.0
    elseif mean_v2 > 0.1
        println("  → Well resolved, diameter ~12-25 mas")
        suggested_size = 18.0
    else
        println("  → Very extended/resolved, diameter >25 mas")
        suggested_size = 30.0
    end

    println("\nRECOMMENDED: Use initial disk diameter = $(suggested_size) mas")
    println("\nNOTE: This file is missing OI_ARRAY table.")
    println("      Image reconstruction with OITOOLS may fail.")
    println("      Consider using parametric model fitting instead.")

    println("\nPress Enter to exit...")
    readline()
    exit(0)
end

println("V² Statistics:")
println("  Min: $(round(minimum(data.v2), digits=4))")
println("  Max: $(round(maximum(data.v2), digits=4))")
println("  Mean: $(round(mean(data.v2), digits=4))")
println("  Median: $(round(median(data.v2), digits=4))")

println("\nWhat this means:")
if mean(data.v2) > 0.7
    println("  → Very compact source, diameter ~1-3 mas")
    suggested_size = 2.0
elseif mean(data.v2) > 0.4
    println("  → Compact source, diameter ~3-6 mas")
    suggested_size = 4.0
elseif mean(data.v2) > 0.2
    println("  → Moderately resolved, diameter ~6-12 mas")
    suggested_size = 8.0
elseif mean(data.v2) > 0.1
    println("  → Well resolved, diameter ~12-25 mas")
    suggested_size = 18.0
else
    println("  → Very extended/resolved, diameter >25 mas")
    suggested_size = 30.0
end

println("\nRECOMMENDED: Use initial disk diameter = $(suggested_size) mas")

println("\nClosure phase RMS: $(round(std(data.t3phi), digits=1))°")
if std(data.t3phi) < 5
    println("  → Nearly symmetric source")
elseif std(data.t3phi) < 20
    println("  → Some asymmetry")
else
    println("  → Highly asymmetric")
end

println("\nPress Enter to exit...")
readline()
