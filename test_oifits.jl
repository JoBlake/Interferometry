# Test what's available in OIFITS
using OIFITS

println("OIFITS package loaded")
println("Available names in OIFITS:")
for name in names(OIFITS)
    println("  - ", name)
end
