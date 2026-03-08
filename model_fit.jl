# Alternative to image reconstruction: Parametric model fitting
# This fits simple geometric models (disk, binary, etc.) to avoid reconstruction artifacts

using OITOOLS
using FITSIO
using OIFITS
using Plots

# Add compatibility workarounds (same as reconstruct.jl)
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
println("PARAMETRIC MODEL FITTING (Alternative to Image Reconstruction)")
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

# Load data
println("\n" * "="^70)
println("Loading OIFITS data from: $filename")
println("="^70)

data = nothing
try
    data = readoifits(filename)
    if isa(data, AbstractArray)
        data = data[1]
    end
    println("  ✓ Loaded $(data.nv2) V² + $(data.nt3amp) T3amp + $(data.nt3phi) T3phi measurements")
catch e
    if occursin("OI_ARRAY", string(e))
        println("\n⚠ ERROR: This OIFITS file is missing the OI_ARRAY table.")
        println("\nThe OI_ARRAY table contains telescope array information that OITOOLS")
        println("requires for model fitting. Without it, parametric fitting cannot proceed.")
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

# Try fitting a simple uniform disk model
println("\n" * "="^70)
println("MODEL 1: UNIFORM DISK")
println("="^70)

println("\nFitting uniform disk model to V² data...")
println("(This avoids closure phase artifacts)\n")

# Create uniform disk model
model = create_model()
component = create_component(type="ud")  # Uniform disk

# Set initial parameters
set_param!(component, "diameter", 6.0)  # Start with 6 mas
set_param!(component, "f", 1.0)         # Flux = 1

add_component!(model, component)

println("Initial parameters:")
println("  Diameter: 6.0 mas")
println("  Flux: 1.0")

# Fit to V² data only
println("\nFitting to V² data...")
result = fit_model_nlopt(model, data, weights=[1.0, 0.0, 0.0])  # V² only

println("\n" * "-"^70)
println("BEST-FIT UNIFORM DISK:")
println("-"^70)
println("  Diameter: $(round(result.params[1], digits=3)) mas")
println("  Flux: $(round(result.params[2], digits=3))")
println("  Chi²: $(round(result.chi2, digits=2))")
println("  Reduced chi²: $(round(result.chi2/data.nv2, digits=2))")

# Try fitting an elliptical disk
println("\n" * "="^70)
println("MODEL 2: ELLIPTICAL DISK")
println("="^70)

model2 = create_model()
component2 = create_component(type="ellipse_uniform")

# Set initial parameters [major_axis, axis_ratio, PA, flux]
set_param!(component2, "diameter", 7.0)
set_param!(component2, "ratio", 0.7)    # b/a = 0.7 (30% elongated)
set_param!(component2, "pa", 45.0)      # Position angle
set_param!(component2, "f", 1.0)

add_component!(model2, component2)

println("\nFitting elliptical disk to V² data...")
result2 = fit_model_nlopt(model2, data, weights=[1.0, 0.0, 0.0])

println("\n" * "-"^70)
println("BEST-FIT ELLIPTICAL DISK:")
println("-"^70)
println("  Major axis: $(round(result2.params[1], digits=3)) mas")
println("  Axis ratio (b/a): $(round(result2.params[2], digits=3))")
println("  Position angle: $(round(result2.params[3], digits=1))°")
println("  Flux: $(round(result2.params[4], digits=3))")
println("  Chi²: $(round(result2.chi2, digits=2))")
println("  Reduced chi²: $(round(result2.chi2/data.nv2, digits=2))")

println("\n" * "="^70)
println("SUMMARY")
println("="^70)

if result2.chi2 < result.chi2 * 0.95
    println("\n✓ Elliptical disk fits significantly better than circular disk")
    println("  → Source is elongated")
    ellipticity = 1.0 / result2.params[2]
    println("  → Ellipticity (major/minor): $(round(ellipticity, digits=2))")
else
    println("\n✓ Circular disk fits data adequately")
    println("  → Source is nearly circular")
end

println("\nNOTE: These fits use V² data only to avoid closure phase artifacts.")
println("      For more complex models, try adding binary components or asymmetries.")

println("\n" * "="^70)
println("Press Enter to exit...")
readline()
