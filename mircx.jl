"""
# MIRCX.jl

A Julia module for astronomy and astrophysics calculations, providing computational
tools for astronomical data processing and analysis.

Version: 0.1.0
"""
module MIRCX

using Dates
using OIFITS
using Plots

# Exports - Coordinate Transformations
export ra_dec_to_galactic, galactic_to_ra_dec

# Exports - Distance and Magnitude
export distance_modulus, apparent_to_absolute_mag, absolute_to_apparent_mag
export mag_to_luminosity, luminosity_to_mag

# Exports - Celestial Mechanics
export angular_separation, position_angle
export proper_motion_total, radial_velocity_correction

# Exports - Time Conversions
export jd_to_calendar, calendar_to_jd, mjd_to_jd, jd_to_mjd
export calculate_gst, calculate_lst

# Exports - Observational Astronomy
export calculate_airmass, atmospheric_extinction
export altitude_to_airmass, hour_angle

# Exports - Utilities
export deg_to_rad, rad_to_deg, hours_to_deg, deg_to_hours
export initialize, generate_oifits_plots

# Physical and Astronomical Constants
const C_LIGHT = 299792458.0              # Speed of light in m/s
const PC_TO_M = 3.0856775814913673e16    # Parsec to meters
const AU_TO_M = 1.495978707e11           # Astronomical Unit to meters
const MJD_OFFSET = 2400000.5             # Offset between JD and MJD

# Coordinate transformation constants
const GALACTIC_POLE_RA = 192.85948      # RA of galactic north pole (degrees, J2000)
const GALACTIC_POLE_DEC = 27.12825      # Dec of galactic north pole (degrees, J2000)
const GALACTIC_CENTER_LON = 122.932     # Galactic longitude of celestial pole (degrees)

# Global variable for OIFITS filename
OIFITS_FILENAME = ""

#==============================================================================#
# INITIALIZATION AND SETUP
#==============================================================================#

"""
    generate_oifits_plots(filename::String)

Read OIFITS file and generate diagnostic plots.

# Arguments
- `filename::String`: Path to OIFITS file

# Returns
- `Bool`: true if plots were generated successfully, false otherwise

# Generates
- uv_coverage.png: UV plane coverage plot
- visibility_squared.png: Squared visibility vs spatial frequency
- closure_phase.png: Closure phase measurements
- visibility_multiwav.png: Visibility at multiple wavelengths
"""
function generate_oifits_plots(filename::String)
    try
        # Read OIFITS file
        oifits_data = OIFITS.load(filename)

        # Extract data tables
        vis2_tables = [table for table in oifits_data if isa(table, OIFITS.OIVis2)]
        vis_tables = [table for table in oifits_data if isa(table, OIFITS.OIVis)]
        t3_tables = [table for table in oifits_data if isa(table, OIFITS.OIT3)]

        println("\nGenerating plots from OIFITS data...")

        # Plot 1: UV Coverage
        if !isempty(vis2_tables) || !isempty(vis_tables)
            p1 = plot(title="UV Coverage", xlabel="U (m)", ylabel="V (m)",
                     aspect_ratio=:equal, legend=false)

            for table in vis2_tables
                scatter!(p1, table.ucoord, table.vcoord, markersize=2, alpha=0.5)
                scatter!(p1, -table.ucoord, -table.vcoord, markersize=2, alpha=0.5)
            end

            for table in vis_tables
                scatter!(p1, table.ucoord, table.vcoord, markersize=2, alpha=0.5)
                scatter!(p1, -table.ucoord, -table.vcoord, markersize=2, alpha=0.5)
            end

            savefig(p1, "uv_coverage.png")
        end

        # Plot 2: Visibility Squared
        if !isempty(vis2_tables)
            p2 = plot(title="Squared Visibility", xlabel="Spatial Frequency (cycles/rad)",
                     ylabel="V²", legend=false)

            for table in vis2_tables
                spatial_freq = sqrt.(table.ucoord.^2 .+ table.vcoord.^2) ./
                              mean(table.eff_wave)
                scatter!(p2, spatial_freq, table.vis2data,
                        yerror=table.vis2err, markersize=3, alpha=0.6)
            end

            savefig(p2, "visibility_squared.png")
        end

        # Plot 3: Closure Phase
        if !isempty(t3_tables)
            p3 = plot(title="Closure Phase", xlabel="Measurement Index",
                     ylabel="Closure Phase (degrees)", legend=false)

            for table in t3_tables
                errorbar_data = table.t3phi
                scatter!(p3, 1:length(errorbar_data), errorbar_data,
                        yerror=table.t3phierr, markersize=3, alpha=0.6)
            end

            savefig(p3, "closure_phase.png")
        end

        # Plot 4: Visibility at Multiple Wavelengths
        if !isempty(vis_tables)
            p4 = plot(title="Visibility (Multi-wavelength)",
                     xlabel="Spatial Frequency (cycles/rad)",
                     ylabel="|V|", legend=:topright)

            for table in vis_tables
                for i in 1:size(table.visamp, 2)
                    spatial_freq = sqrt.(table.ucoord.^2 .+ table.vcoord.^2) ./
                                  table.eff_wave[i]
                    scatter!(p4, spatial_freq, table.visamp[:, i],
                            label="λ = $(round(table.eff_wave[i]*1e6, digits=2)) μm",
                            markersize=3, alpha=0.6)
                end
            end

            savefig(p4, "visibility_multiwav.png")
        end

        println("\nGenerated plots:")
        println("  - uv_coverage.png")
        println("  - visibility_squared.png")
        println("  - closure_phase.png")
        println("  - visibility_multiwav.png")

        return true

    catch e
        @error "Failed to generate plots: $e"
        return false
    end
end

"""
    initialize()

Initialize the MIRCX module by prompting user for OIFITS file name.
This should be called when starting a new session.

# Example
```julia
using MIRCX
initialize()
```
"""
function initialize()
    println("MIRCX.jl - Astronomy and Astrophysics Calculation Module")
    println("=" ^ 60)
    print("Enter OIFITS file name (or press Enter to skip): ")
    global OIFITS_FILENAME = readline()

    if isempty(OIFITS_FILENAME)
        println("No OIFITS file specified. You can set it later.")
    else
        if isfile(OIFITS_FILENAME)
            println("OIFITS file set to: $OIFITS_FILENAME")
            println("\nReading OIFITS file and generating plots...")
            generate_oifits_plots(OIFITS_FILENAME)
        else
            @warn "File '$OIFITS_FILENAME' not found. Please verify the path."
        end
    end
    println("=" ^ 60)
    return OIFITS_FILENAME
end

#==============================================================================#
# UTILITY FUNCTIONS
#==============================================================================#

"""
    deg_to_rad(degrees::Real)

Convert degrees to radians.

# Arguments
- `degrees::Real`: Angle in degrees

# Returns
- `Float64`: Angle in radians

# Example
```julia
rad = deg_to_rad(180.0)  # Returns π
```
"""
deg_to_rad(degrees::Real) = degrees * π / 180.0

"""
    rad_to_deg(radians::Real)

Convert radians to degrees.

# Arguments
- `radians::Real`: Angle in radians

# Returns
- `Float64`: Angle in degrees

# Example
```julia
deg = rad_to_deg(π)  # Returns 180.0
```
"""
rad_to_deg(radians::Real) = radians * 180.0 / π

"""
    hours_to_deg(hours::Real)

Convert hours to degrees (1 hour = 15 degrees).

# Arguments
- `hours::Real`: Angle in hours

# Returns
- `Float64`: Angle in degrees
"""
hours_to_deg(hours::Real) = hours * 15.0

"""
    deg_to_hours(degrees::Real)

Convert degrees to hours (15 degrees = 1 hour).

# Arguments
- `degrees::Real`: Angle in degrees

# Returns
- `Float64`: Angle in hours
"""
deg_to_hours(degrees::Real) = degrees / 15.0

#==============================================================================#
# COORDINATE TRANSFORMATIONS
#==============================================================================#

"""
    ra_dec_to_galactic(ra::Real, dec::Real; degrees=true)

Convert equatorial coordinates (RA, Dec) to galactic coordinates (l, b).

# Arguments
- `ra::Real`: Right Ascension
- `dec::Real`: Declination
- `degrees::Bool`: If true (default), inputs/outputs are in degrees; otherwise radians

# Returns
- `Tuple{Float64, Float64}`: Galactic longitude (l) and latitude (b)

# Example
```julia
l, b = ra_dec_to_galactic(266.4, -29.0)  # Galactic center
```

# Reference
Based on the transformation matrix from J2000 equatorial to galactic coordinates.
"""
function ra_dec_to_galactic(ra::Real, dec::Real; degrees=true)
    # Convert to radians if needed
    ra_rad = degrees ? deg_to_rad(ra) : ra
    dec_rad = degrees ? deg_to_rad(dec) : dec

    # Validate declination range
    if degrees && (dec < -90.0 || dec > 90.0)
        throw(ArgumentError("Declination must be between -90 and 90 degrees"))
    end

    # Galactic pole position
    pole_ra = deg_to_rad(GALACTIC_POLE_RA)
    pole_dec = deg_to_rad(GALACTIC_POLE_DEC)
    lon_cp = deg_to_rad(GALACTIC_CENTER_LON)

    # Calculate galactic latitude
    sin_b = sin(dec_rad) * sin(pole_dec) +
            cos(dec_rad) * cos(pole_dec) * cos(ra_rad - pole_ra)
    b = asin(sin_b)

    # Calculate galactic longitude
    y = cos(dec_rad) * sin(ra_rad - pole_ra)
    x = sin(dec_rad) * cos(pole_dec) -
        cos(dec_rad) * sin(pole_dec) * cos(ra_rad - pole_ra)
    l = lon_cp - atan(y, x)

    # Normalize longitude to [0, 2π)
    l = mod(l, 2π)

    return degrees ? (rad_to_deg(l), rad_to_deg(b)) : (l, b)
end

"""
    galactic_to_ra_dec(l::Real, b::Real; degrees=true)

Convert galactic coordinates (l, b) to equatorial coordinates (RA, Dec).

# Arguments
- `l::Real`: Galactic longitude
- `b::Real`: Galactic latitude
- `degrees::Bool`: If true (default), inputs/outputs are in degrees; otherwise radians

# Returns
- `Tuple{Float64, Float64}`: Right Ascension (RA) and Declination (Dec)

# Example
```julia
ra, dec = galactic_to_ra_dec(0.0, 0.0)  # Galactic center direction
```
"""
function galactic_to_ra_dec(l::Real, b::Real; degrees=true)
    # Convert to radians if needed
    l_rad = degrees ? deg_to_rad(l) : l
    b_rad = degrees ? deg_to_rad(b) : b

    # Galactic pole position
    pole_ra = deg_to_rad(GALACTIC_POLE_RA)
    pole_dec = deg_to_rad(GALACTIC_POLE_DEC)
    lon_cp = deg_to_rad(GALACTIC_CENTER_LON)

    # Calculate declination
    sin_dec = sin(b_rad) * sin(pole_dec) +
              cos(b_rad) * cos(pole_dec) * sin(lon_cp - l_rad)
    dec = asin(sin_dec)

    # Calculate right ascension
    y = cos(b_rad) * cos(lon_cp - l_rad)
    x = sin(b_rad) * cos(pole_dec) -
        cos(b_rad) * sin(pole_dec) * sin(lon_cp - l_rad)
    ra = pole_ra + atan(y, x)

    # Normalize RA to [0, 2π)
    ra = mod(ra, 2π)

    return degrees ? (rad_to_deg(ra), rad_to_deg(dec)) : (ra, dec)
end

#==============================================================================#
# DISTANCE AND MAGNITUDE CALCULATIONS
#==============================================================================#

"""
    distance_modulus(distance_pc::Real)

Calculate the distance modulus for a given distance.

# Arguments
- `distance_pc::Real`: Distance in parsecs

# Returns
- `Float64`: Distance modulus (m - M)

# Example
```julia
μ = distance_modulus(10.0)  # Returns 0.0 for 10 pc
```

# Formula
μ = 5 * log₁₀(d) - 5, where d is distance in parsecs
"""
function distance_modulus(distance_pc::Real)
    if distance_pc <= 0
        throw(ArgumentError("Distance must be positive"))
    end
    return 5.0 * log10(distance_pc) - 5.0
end

"""
    apparent_to_absolute_mag(apparent_mag::Real, distance_pc::Real; extinction=0.0)

Convert apparent magnitude to absolute magnitude.

# Arguments
- `apparent_mag::Real`: Apparent magnitude
- `distance_pc::Real`: Distance in parsecs
- `extinction::Real`: Extinction correction (default: 0.0)

# Returns
- `Float64`: Absolute magnitude

# Example
```julia
M = apparent_to_absolute_mag(5.0, 100.0)  # m=5, d=100pc
```

# Formula
M = m - μ - A, where μ is distance modulus and A is extinction
"""
function apparent_to_absolute_mag(apparent_mag::Real, distance_pc::Real; extinction=0.0)
    μ = distance_modulus(distance_pc)
    return apparent_mag - μ - extinction
end

"""
    absolute_to_apparent_mag(absolute_mag::Real, distance_pc::Real; extinction=0.0)

Convert absolute magnitude to apparent magnitude.

# Arguments
- `absolute_mag::Real`: Absolute magnitude
- `distance_pc::Real`: Distance in parsecs
- `extinction::Real`: Extinction correction (default: 0.0)

# Returns
- `Float64`: Apparent magnitude

# Example
```julia
m = absolute_to_apparent_mag(0.0, 100.0)  # M=0, d=100pc
```
"""
function absolute_to_apparent_mag(absolute_mag::Real, distance_pc::Real; extinction=0.0)
    μ = distance_modulus(distance_pc)
    return absolute_mag + μ + extinction
end

"""
    mag_to_luminosity(mag::Real, mag_ref::Real=0.0, L_ref::Real=1.0)

Convert magnitude to luminosity relative to a reference.

# Arguments
- `mag::Real`: Magnitude to convert
- `mag_ref::Real`: Reference magnitude (default: 0.0)
- `L_ref::Real`: Reference luminosity (default: 1.0)

# Returns
- `Float64`: Luminosity in units of L_ref

# Formula
L/L_ref = 10^(-0.4 * (M - M_ref))
"""
function mag_to_luminosity(mag::Real, mag_ref::Real=0.0, L_ref::Real=1.0)
    return L_ref * 10.0^(-0.4 * (mag - mag_ref))
end

"""
    luminosity_to_mag(luminosity::Real, mag_ref::Real=0.0, L_ref::Real=1.0)

Convert luminosity to magnitude.

# Arguments
- `luminosity::Real`: Luminosity in units of L_ref
- `mag_ref::Real`: Reference magnitude (default: 0.0)
- `L_ref::Real`: Reference luminosity (default: 1.0)

# Returns
- `Float64`: Magnitude
"""
function luminosity_to_mag(luminosity::Real, mag_ref::Real=0.0, L_ref::Real=1.0)
    if luminosity <= 0
        throw(ArgumentError("Luminosity must be positive"))
    end
    return mag_ref - 2.5 * log10(luminosity / L_ref)
end

#==============================================================================#
# CELESTIAL MECHANICS
#==============================================================================#

"""
    angular_separation(ra1::Real, dec1::Real, ra2::Real, dec2::Real; degrees=true)

Calculate the angular separation between two celestial positions.

# Arguments
- `ra1::Real, dec1::Real`: First position (RA, Dec)
- `ra2::Real, dec2::Real`: Second position (RA, Dec)
- `degrees::Bool`: If true (default), inputs/outputs are in degrees; otherwise radians

# Returns
- `Float64`: Angular separation

# Example
```julia
sep = angular_separation(10.0, 20.0, 15.0, 25.0)
```

# Reference
Uses the Vincenty formula for better numerical stability.
"""
function angular_separation(ra1::Real, dec1::Real, ra2::Real, dec2::Real; degrees=true)
    # Convert to radians if needed
    α1 = degrees ? deg_to_rad(ra1) : ra1
    δ1 = degrees ? deg_to_rad(dec1) : dec1
    α2 = degrees ? deg_to_rad(ra2) : ra2
    δ2 = degrees ? deg_to_rad(dec2) : dec2

    # Vincenty formula for better numerical stability
    Δα = α2 - α1

    num = sqrt((cos(δ2) * sin(Δα))^2 +
               (cos(δ1) * sin(δ2) - sin(δ1) * cos(δ2) * cos(Δα))^2)
    den = sin(δ1) * sin(δ2) + cos(δ1) * cos(δ2) * cos(Δα)

    separation = atan(num, den)

    return degrees ? rad_to_deg(separation) : separation
end

"""
    position_angle(ra1::Real, dec1::Real, ra2::Real, dec2::Real; degrees=true)

Calculate the position angle from position 1 to position 2.

# Arguments
- `ra1::Real, dec1::Real`: First position (RA, Dec)
- `ra2::Real, dec2::Real`: Second position (RA, Dec)
- `degrees::Bool`: If true (default), inputs/outputs are in degrees; otherwise radians

# Returns
- `Float64`: Position angle measured East of North (0-360 degrees or 0-2π radians)

# Example
```julia
pa = position_angle(10.0, 20.0, 15.0, 25.0)
```
"""
function position_angle(ra1::Real, dec1::Real, ra2::Real, dec2::Real; degrees=true)
    # Convert to radians if needed
    α1 = degrees ? deg_to_rad(ra1) : ra1
    δ1 = degrees ? deg_to_rad(dec1) : dec1
    α2 = degrees ? deg_to_rad(ra2) : ra2
    δ2 = degrees ? deg_to_rad(dec2) : dec2

    Δα = α2 - α1

    y = sin(Δα)
    x = cos(δ1) * tan(δ2) - sin(δ1) * cos(Δα)

    pa = atan(y, x)

    # Normalize to [0, 2π)
    pa = mod(pa, 2π)

    return degrees ? rad_to_deg(pa) : pa
end

"""
    proper_motion_total(pm_ra::Real, pm_dec::Real)

Calculate total proper motion from RA and Dec components.

# Arguments
- `pm_ra::Real`: Proper motion in RA (mas/yr, includes cos(dec) factor)
- `pm_dec::Real`: Proper motion in Dec (mas/yr)

# Returns
- `Float64`: Total proper motion (mas/yr)

# Example
```julia
pm_total = proper_motion_total(10.0, 15.0)
```
"""
function proper_motion_total(pm_ra::Real, pm_dec::Real)
    return sqrt(pm_ra^2 + pm_dec^2)
end

"""
    radial_velocity_correction(rv_helio::Real, v_bary::Real=0.0)

Apply barycentric correction to heliocentric radial velocity.

# Arguments
- `rv_helio::Real`: Heliocentric radial velocity (km/s)
- `v_bary::Real`: Barycentric correction velocity (km/s, default: 0.0)

# Returns
- `Float64`: Barycentric radial velocity (km/s)

# Example
```julia
rv_bary = radial_velocity_correction(50.0, 2.5)
```
"""
function radial_velocity_correction(rv_helio::Real, v_bary::Real=0.0)
    return rv_helio + v_bary
end

#==============================================================================#
# TIME CONVERSIONS
#==============================================================================#

"""
    calendar_to_jd(year::Int, month::Int, day::Int, hour::Int=0, minute::Int=0, second::Real=0.0)

Convert calendar date to Julian Date.

# Arguments
- `year::Int`: Year
- `month::Int`: Month (1-12)
- `day::Int`: Day of month
- `hour::Int`: Hour (0-23, default: 0)
- `minute::Int`: Minute (0-59, default: 0)
- `second::Real`: Second (0-60, default: 0.0)

# Returns
- `Float64`: Julian Date

# Example
```julia
jd = calendar_to_jd(2000, 1, 1, 12, 0, 0)  # J2000.0
```

# Reference
Based on Meeus, "Astronomical Algorithms" (1998)
"""
function calendar_to_jd(year::Int, month::Int, day::Int,
                        hour::Int=0, minute::Int=0, second::Real=0.0)
    # Adjust for January and February being months 13 and 14 of the previous year
    a = div(14 - month, 12)
    y = year + 4800 - a
    m = month + 12 * a - 3

    # Calculate JD at noon
    jd = day + div(153 * m + 2, 5) + 365 * y + div(y, 4) - div(y, 100) + div(y, 400) - 32045

    # Add fractional day
    fraction = (hour - 12) / 24.0 + minute / 1440.0 + second / 86400.0

    return Float64(jd) + fraction
end

"""
    jd_to_calendar(jd::Real)

Convert Julian Date to calendar date.

# Arguments
- `jd::Real`: Julian Date

# Returns
- `DateTime`: Calendar date and time

# Example
```julia
dt = jd_to_calendar(2451545.0)  # Returns DateTime for J2000.0
```
"""
function jd_to_calendar(jd::Real)
    # Algorithm from Meeus
    jd_adjusted = jd + 0.5
    z = floor(Int, jd_adjusted)
    f = jd_adjusted - z

    if z >= 2299161
        α = floor(Int, (z - 1867216.25) / 36524.25)
        a = z + 1 + α - div(α, 4)
    else
        a = z
    end

    b = a + 1524
    c = floor(Int, (b - 122.1) / 365.25)
    d = floor(Int, 365.25 * c)
    e = floor(Int, (b - d) / 30.6001)

    day = b - d - floor(Int, 30.6001 * e)
    month = e < 14 ? e - 1 : e - 13
    year = month > 2 ? c - 4716 : c - 4715

    # Calculate time
    day_fraction = f
    hours = floor(Int, day_fraction * 24)
    minutes = floor(Int, (day_fraction * 24 - hours) * 60)
    seconds = ((day_fraction * 24 - hours) * 60 - minutes) * 60

    return DateTime(year, month, day, hours, minutes, floor(Int, seconds))
end

"""
    jd_to_mjd(jd::Real)

Convert Julian Date to Modified Julian Date.

# Arguments
- `jd::Real`: Julian Date

# Returns
- `Float64`: Modified Julian Date

# Example
```julia
mjd = jd_to_mjd(2451545.0)
```
"""
jd_to_mjd(jd::Real) = jd - MJD_OFFSET

"""
    mjd_to_jd(mjd::Real)

Convert Modified Julian Date to Julian Date.

# Arguments
- `mjd::Real`: Modified Julian Date

# Returns
- `Float64`: Julian Date

# Example
```julia
jd = mjd_to_jd(51544.5)
```
"""
mjd_to_jd(mjd::Real) = mjd + MJD_OFFSET

"""
    calculate_gst(jd::Real; degrees=true)

Calculate Greenwich Sidereal Time.

# Arguments
- `jd::Real`: Julian Date
- `degrees::Bool`: If true (default), output is in degrees; otherwise hours

# Returns
- `Float64`: Greenwich Sidereal Time

# Example
```julia
gst = calculate_gst(2451545.0)
```

# Reference
Simplified formula accurate to ~0.1 seconds
"""
function calculate_gst(jd::Real; degrees=true)
    # Days since J2000.0
    T = (jd - 2451545.0) / 36525.0

    # GST in hours at 0h UT
    gst = 280.46061837 + 360.98564736629 * (jd - 2451545.0) +
          0.000387933 * T^2 - T^3 / 38710000.0

    # Normalize to [0, 360)
    gst = mod(gst, 360.0)

    return degrees ? gst : deg_to_hours(gst)
end

"""
    calculate_lst(jd::Real, longitude::Real; degrees=true)

Calculate Local Sidereal Time.

# Arguments
- `jd::Real`: Julian Date
- `longitude::Real`: Observer's longitude (East positive)
- `degrees::Bool`: If true (default), inputs/outputs are in degrees; otherwise hours

# Returns
- `Float64`: Local Sidereal Time

# Example
```julia
lst = calculate_lst(2451545.0, -77.0)  # Longitude 77°W
```
"""
function calculate_lst(jd::Real, longitude::Real; degrees=true)
    gst = calculate_gst(jd; degrees=true)
    lon = degrees ? longitude : hours_to_deg(longitude)

    lst = gst + lon
    lst = mod(lst, 360.0)

    return degrees ? lst : deg_to_hours(lst)
end

#==============================================================================#
# OBSERVATIONAL ASTRONOMY
#==============================================================================#

"""
    altitude_to_airmass(altitude::Real; degrees=true)

Calculate airmass from altitude angle.

# Arguments
- `altitude::Real`: Altitude angle above horizon
- `degrees::Bool`: If true (default), input is in degrees; otherwise radians

# Returns
- `Float64`: Airmass (dimensionless)

# Example
```julia
X = altitude_to_airmass(30.0)  # 30 degrees altitude
```

# Reference
Uses Hardie (1962) formula for altitudes > 15°, plane-parallel for lower altitudes.
"""
function altitude_to_airmass(altitude::Real; degrees=true)
    alt = degrees ? altitude : rad_to_deg(altitude)

    if alt <= 0
        return Inf
    elseif alt >= 90
        return 1.0
    end

    z = 90.0 - alt  # Zenith angle
    z_rad = deg_to_rad(z)

    # Hardie (1962) formula
    sec_z = 1.0 / cos(z_rad)

    if alt >= 15.0
        # More accurate formula for higher altitudes
        airmass = sec_z - 0.0018167 * (sec_z - 1) -
                  0.002875 * (sec_z - 1)^2 -
                  0.0008083 * (sec_z - 1)^3
    else
        # Simple plane-parallel approximation for low altitudes
        airmass = sec_z
    end

    return airmass
end

"""
    calculate_airmass(altitude::Real; degrees=true)

Alias for altitude_to_airmass for consistency with spec.
"""
calculate_airmass(altitude::Real; degrees=true) = altitude_to_airmass(altitude; degrees=degrees)

"""
    atmospheric_extinction(airmass::Real, k::Real=0.2)

Calculate atmospheric extinction magnitude.

# Arguments
- `airmass::Real`: Airmass value
- `k::Real`: Extinction coefficient in magnitudes (default: 0.2 for V-band)

# Returns
- `Float64`: Extinction in magnitudes

# Example
```julia
A = atmospheric_extinction(2.0, 0.15)  # 2 airmasses, k=0.15
```

# Formula
A = k * X, where X is airmass
"""
function atmospheric_extinction(airmass::Real, k::Real=0.2)
    if airmass < 0
        throw(ArgumentError("Airmass must be non-negative"))
    end
    return k * airmass
end

"""
    hour_angle(lst::Real, ra::Real; degrees=true)

Calculate hour angle from local sidereal time and right ascension.

# Arguments
- `lst::Real`: Local Sidereal Time
- `ra::Real`: Right Ascension
- `degrees::Bool`: If true (default), inputs/outputs are in degrees; otherwise hours

# Returns
- `Float64`: Hour angle

# Example
```julia
ha = hour_angle(120.0, 100.0)  # LST=120°, RA=100°
```
"""
function hour_angle(lst::Real, ra::Real; degrees=true)
    ha = lst - ra

    if degrees
        # Normalize to [-180, 180]
        ha = mod(ha + 180.0, 360.0) - 180.0
    else
        # Normalize to [-12, 12] hours
        ha = mod(ha + 12.0, 24.0) - 12.0
    end

    return ha
end

#==============================================================================#
# MODULE INITIALIZATION MESSAGE
#==============================================================================#

function __init__()
    initialize()
end

end # module MIRCX
