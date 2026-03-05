# Nico interferometry
using FITSIO
using Printf

# Open the FITS file
filename = "C:/Users/johnb/OneDrive/Documents/Projects/Astronomy/MIRCX_L2.2025May20.V_CVn.MIRCX_IDL.GHS.AVG6m.oifits"
f = FITS(filename)

println("="^80)
println("OIFITS File Analysis")
println("File: ", basename(filename))
println("="^80)
println()

# Display all HDUs in the file
println("File Structure:")
println("-"^80)
for i in 1:length(f)
    hdu = f[i]
    hdu_type = typeof(hdu)
    
    # Get extension name if it exists
    extname = ""
    try
        extname = read_key(hdu, "EXTNAME")[1]
    catch
        extname = "PRIMARY"
    end
    
    println("HDU $i: $extname ($(hdu_type))")
end
println()

# Function to display header keywords
function display_header(hdu, title)
    println("="^80)
    println(title)
    println("="^80)
    
    header = read_header(hdu)
    
    # Get all keywords
    for key in keys(header)
        # Skip some verbose keywords
        if key ∉ ["COMMENT", "HISTORY", ""]
            value = header[key]
            # Get comment if available
            try
                comment = get_comment(header, key)
                println("  $key = $value")
                if !isempty(comment)
                    println("      # $comment")
                end
            catch
                println("  $key = $value")
            end
        end
    end
    println()
end

# 1. PRIMARY HDU
display_header(f[1], "PRIMARY HDU")

# Function to find HDU by extension name
function find_hdu(fits_file, extname)
    for i in 1:length(fits_file)
        try
            if read_key(fits_file[i], "EXTNAME")[1] == extname
                return i
            end
        catch
            continue
        end
    end
    return nothing
end

# 2. OI_TARGET
target_idx = find_hdu(f, "OI_TARGET")
if target_idx !== nothing
    println("="^80)
    println("OI_TARGET - Observed Targets")
    println("="^80)
    
    hdu = f[target_idx]
    
    # Read the table columns
    target_id = read(hdu, "TARGET_ID")
    target = read(hdu, "TARGET")
    raep0 = read(hdu, "RAEP0")
    decep0 = read(hdu, "DECEP0")
    equinox = read(hdu, "EQUINOX")
    
    println("Number of targets: ", length(target))
    println()
    
    for i in 1:length(target)
        println("Target $i:")
        println("  ID: ", target_id[i])
        println("  Name: ", strip(target[i]))
        println("  RA (J2000): ", raep0[i], " degrees")
        println("  Dec (J2000): ", decep0[i], " degrees")
        println("  Equinox: ", equinox[i])
        println()
    end
else
    println("OI_TARGET not found in file")
    println()
end

# 3. OI_ARRAY
array_idx = find_hdu(f, "OI_ARRAY")
if array_idx !== nothing
    println("="^80)
    println("OI_ARRAY - Telescope Array Configuration")
    println("="^80)
    
    hdu = f[array_idx]
    
    # Display header info
    try
        arrname = read_key(hdu, "ARRNAME")[1]
        frame = read_key(hdu, "FRAME")[1]
        arrayx = read_key(hdu, "ARRAYX")[1]
        arrayy = read_key(hdu, "ARRAYY")[1]
        arrayz = read_key(hdu, "ARRAYZ")[1]
        
        println("Array Name: ", arrname)
        println("Reference Frame: ", frame)
        println("Array Center (X, Y, Z): ($arrayx, $arrayy, $arrayz) meters")
        println()
    catch
    end
    
    # Read station data
    tel_name = read(hdu, "TEL_NAME")
    sta_name = read(hdu, "STA_NAME")
    staxyz = read(hdu, "STAXYZ")
    diameter = read(hdu, "DIAMETER")
    
    println("Number of stations: ", length(tel_name))
    println()
    
    for i in 1:length(tel_name)
        println("Station $i:")
        println("  Telescope: ", strip(tel_name[i]))
        println("  Station: ", strip(sta_name[i]))
        println("  Position (X,Y,Z): (", staxyz[1,i], ", ", staxyz[2,i], ", ", staxyz[3,i], ") meters")
        println("  Diameter: ", diameter[i], " meters")
        println()
    end
else
    println("OI_ARRAY not found in file")
    println()
end

# 4. OI_WAVELENGTH
wave_idx = find_hdu(f, "OI_WAVELENGTH")
if wave_idx !== nothing
    println("="^80)
    println("OI_WAVELENGTH - Spectral Channels")
    println("="^80)
    
    hdu = f[wave_idx]
    
    # Display header info
    try
        insname = read_key(hdu, "INSNAME")[1]
        println("Instrument: ", insname)
        println()
    catch
    end
    
    # Read wavelength data
    eff_wave = read(hdu, "EFF_WAVE")
    eff_band = read(hdu, "EFF_BAND")
    
    println("Number of spectral channels: ", length(eff_wave))
    println()
    
    # Convert to microns for readability
    println("Channel | Wavelength (μm) | Bandwidth (μm)")
    println("-"^50)
    for i in 1:length(eff_wave)
        wave_micron = eff_wave[i] * 1e6
        band_micron = eff_band[i] * 1e6
        println(@sprintf("%7d | %15.6f | %14.6f", i, wave_micron, band_micron))
    end
    println()
    
    println("Wavelength range: ", minimum(eff_wave)*1e6, " - ", maximum(eff_wave)*1e6, " μm")
    println()
else
    println("OI_WAVELENGTH not found in file")
    println()
end

# 5. OI_VIS2 - Squared Visibility Data
vis2_idx = find_hdu(f, "OI_VIS2")
if vis2_idx !== nothing
    println("="^80)
    println("OI_VIS2 - Squared Visibility Data")
    println("="^80)
    
    hdu = f[vis2_idx]
    
    # Read data
    vis2data = read(hdu, "VIS2DATA")
    vis2err = read(hdu, "VIS2ERR")
    ucoord = read(hdu, "UCOORD")
    vcoord = read(hdu, "VCOORD")
    
    # Calculate spatial frequencies (baseline lengths)
    baseline = sqrt.(ucoord.^2 .+ vcoord.^2)
    
    println("Number of measurements: ", length(ucoord))
    println("Number of spectral channels: ", size(vis2data, 1))
    println()
    
    # Statistics
    println("Baseline range: ", minimum(baseline), " - ", maximum(baseline), " meters")
    println("V² range: ", minimum(vis2data), " - ", maximum(vis2data))
    println()
else
    println("OI_VIS2 not found in file")
    println()
end

# 6. OI_T3 - Closure Phase Data
t3_idx = find_hdu(f, "OI_T3")
if t3_idx !== nothing
    println("="^80)
    println("OI_T3 - Closure Phase Data")
    println("="^80)
    
    hdu = f[t3_idx]
    
    # Read data
    t3amp = read(hdu, "T3AMP")
    t3phi = read(hdu, "T3PHI")
    u1coord = read(hdu, "U1COORD")
    v1coord = read(hdu, "V1COORD")
    u2coord = read(hdu, "U2COORD")
    v2coord = read(hdu, "V2COORD")
    
    # Calculate baseline lengths for the triangle
    baseline1 = sqrt.(u1coord.^2 .+ v1coord.^2)
    baseline2 = sqrt.(u2coord.^2 .+ v2coord.^2)
    
    println("Number of closure phase measurements: ", length(u1coord))
    println("Number of spectral channels: ", size(t3phi, 1))
    println()
    
    # Statistics
    println("Baseline 1 range: ", minimum(baseline1), " - ", maximum(baseline1), " meters")
    println("Baseline 2 range: ", minimum(baseline2), " - ", maximum(baseline2), " meters")
    println("Closure phase range: ", minimum(t3phi), " - ", maximum(t3phi), " degrees")
    println("T3 amplitude range: ", minimum(t3amp), " - ", maximum(t3amp))
    println()
else
    println("OI_T3 not found in file")
    println()
end

# Close the file
close(f)

println("="^80)
println("Creating visualizations...")
println("="^80)

# Now create plots using Plots.jl
using Plots

# Reopen file to get data for plotting
f = FITS(filename)

# Plot 1: UV Coverage
println("Generating UV coverage plot...")
vis2_idx = find_hdu(f, "OI_VIS2")
if vis2_idx !== nothing
    hdu = f[vis2_idx]
    ucoord = read(hdu, "UCOORD")
    vcoord = read(hdu, "VCOORD")
    
    p1 = scatter(ucoord, vcoord, 
                 label="Observed baselines",
                 marker=:circle,
                 markersize=3,
                 xlabel="U (meters)",
                 ylabel="V (meters)",
                 title="UV Coverage",
                 aspect_ratio=:equal,
                 grid=true)
    
    # Add mirror points (interferometry symmetry)
    scatter!(p1, -ucoord, -vcoord, 
            label="Mirror baselines",
            marker=:circle,
            markersize=3,
            alpha=0.5)
    
    savefig(p1, "uv_coverage.png")
    println("Saved: uv_coverage.png")
end

# Plot 2: Squared Visibility vs Baseline
println("Generating V² vs baseline plot...")
if vis2_idx !== nothing
    hdu = f[vis2_idx]
    vis2data = read(hdu, "VIS2DATA")
    vis2err = read(hdu, "VIS2ERR")
    ucoord = read(hdu, "UCOORD")
    vcoord = read(hdu, "VCOORD")
    
    baseline = sqrt.(ucoord.^2 .+ vcoord.^2)
    
    # Convert to spatial frequency (cycles/radian)
    wave_idx = find_hdu(f, "OI_WAVELENGTH")
    if wave_idx !== nothing
        wave_hdu = f[wave_idx]
        eff_wave = read(wave_hdu, "EFF_WAVE")
        
        # Use first wavelength channel for plotting
        spatial_freq = baseline ./ eff_wave[1]
        
        # Plot first spectral channel
        p2 = scatter(spatial_freq, vis2data[1, :],
                    yerr=vis2err[1, :],
                    label="Channel 1",
                    marker=:circle,
                    markersize=4,
                    xlabel="Spatial Frequency (cycles/rad)",
                    ylabel="V²",
                    title="Squared Visibility",
                    legend=:topright,
                    grid=true,
                    ylims=(0, 1.1))
        
        savefig(p2, "visibility_squared.png")
        println("Saved: visibility_squared.png")
    end
end

# Plot 3: Closure Phase vs Time or Baseline
println("Generating closure phase plot...")
t3_idx = find_hdu(f, "OI_T3")
if t3_idx !== nothing
    hdu = f[t3_idx]
    t3phi = read(hdu, "T3PHI")
    t3phierr = read(hdu, "T3PHIERR")
    u1coord = read(hdu, "U1COORD")
    v1coord = read(hdu, "V1COORD")
    
    # Maximum baseline for each triangle
    baseline_max = sqrt.(u1coord.^2 .+ v1coord.^2)
    
    # Plot first spectral channel
    p3 = scatter(baseline_max, t3phi[1, :],
                yerr=t3phierr[1, :],
                label="Channel 1",
                marker=:circle,
                markersize=4,
                xlabel="Maximum Baseline (meters)",
                ylabel="Closure Phase (degrees)",
                title="Closure Phase",
                legend=:topright,
                grid=true)
    
    # Add zero line
    hline!(p3, [0], linestyle=:dash, color=:black, label="Zero")
    
    savefig(p3, "closure_phase.png")
    println("Saved: closure_phase.png")
end

# Plot 4: Multi-wavelength V² plot
println("Generating multi-wavelength V² plot...")
if vis2_idx !== nothing && wave_idx !== nothing
    vis2_hdu = f[vis2_idx]
    wave_hdu = f[wave_idx]
    
    vis2data = read(vis2_hdu, "VIS2DATA")
    ucoord = read(vis2_hdu, "UCOORD")
    vcoord = read(vis2_hdu, "VCOORD")
    eff_wave = read(wave_hdu, "EFF_WAVE")
    
    baseline = sqrt.(ucoord.^2 .+ vcoord.^2)
    
    # Create multi-channel plot
    p4 = plot(xlabel="Baseline (meters)",
              ylabel="V²",
              title="Squared Visibility - All Channels",
              legend=:topright,
              grid=true)
    
    # Plot up to 8 channels to avoid clutter
    n_channels = min(8, size(vis2data, 1))
    colors = palette(:rainbow, n_channels)
    
    for i in 1:n_channels
        scatter!(p4, baseline, vis2data[i, :],
                marker=:circle,
                markersize=2,
                alpha=0.6,
                color=colors[i],
                label=@sprintf("%.3f μm", eff_wave[i]*1e6))
    end
    
    savefig(p4, "visibility_multiwav.png")
    println("Saved: visibility_multiwav.png")
end

close(f)

println()
println("Generated plots:")
println("  - uv_coverage.png")
println("  - visibility_squared.png")
println("  - closure_phase.png")
println("  - visibility_multiwav.png")
println()

# Display plots in interactive windows
println("Displaying plots in interactive windows...")
println("Close the plot windows to exit.")
println()

# Show all plots together
plot(p1, p2, p3, p4, layout=(2,2), size=(1200, 1000))
gui()  # Display the plots and keep window open

println("Press Enter to exit...")
readline()