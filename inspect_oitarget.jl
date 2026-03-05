using OIFITS
using AstroFITS

# Open the OIFITS file
filename = "2009-11-eps_Aur-avg5.oifits"

if isfile(filename)
    println("Reading $filename\n")

    # Read the dataset
    ds = OIFITS.read(OIFITS.OIDataSet, filename; hack_revn=1)

    println("Dataset loaded successfully!")
    println("\nOI_TARGET structure:")
    println("  Type: ", typeof(ds.target))
    println("  Fields: ", fieldnames(typeof(ds.target)))

    println("\nInspecting ds.target:")
    for field in fieldnames(typeof(ds.target))
        println("  $field = ", getfield(ds.target, field))
    end

    # Check if target has a list
    if hasfield(typeof(ds.target), :list)
        list = getfield(ds.target, :list)
        println("\n  list type: ", typeof(list))
        if length(list) > 0
            println("  list[1] type: ", typeof(list[1]))
            println("  list[1] fields: ", fieldnames(typeof(list[1])))
            println("  list[1] contents: ", list[1])
        end
    end
else
    println("File not found: $filename")
    println("Available files:")
    for f in readdir(".")
        if endswith(f, ".oifits")
            println("  - $f")
        end
    end
end
