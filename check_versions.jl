using Pkg

println("Checking installed package versions:\n")

# Get package info
deps = Pkg.dependencies()

packages = ["OITOOLS", "OIFITS", "OIFITS", "AstroFITS", "FITSIO"]

for pkg_name in packages
    for (uuid, dep) in deps
        if dep.name == pkg_name
            println("$pkg_name:")
            println("  Version: $(dep.version)")
            println("  Source: $(dep.source)")
            println()
        end
    end
end

# Check what functions OIFITS actually exports
println("\n" * "="^60)
println("Functions/types exported by OIFITS:")
println("="^60)
using OIFITS
for name in names(OIFITS, all=false)
    println("  - $name")
end

println("\n" * "="^60)
println("Functions/types exported by OITOOLS:")
println("="^60)
using OITOOLS
for name in names(OITOOLS, all=false)
    println("  - $name")
end
