# MIRCX.jl Specification

## Overview
MIRCX.jl is a Julia module for astronomy and astrophysics calculations, providing computational tools for astronomical data processing and analysis.

## Purpose
Provide a robust, efficient set of functions for common astronomical computations used in observational astronomy and astrophysics research.

## Requirements

### Functional Requirements

#### 1. Coordinate Transformations
- Convert between equatorial (RA/Dec) and galactic coordinates
- Convert between different coordinate systems (ICRS, FK5, etc.)
- Handle proper motion corrections
- Support epoch transformations (J2000, B1950, etc.)

#### 2. Distance and Magnitude Calculations
- Calculate distance modulus
- Convert between apparent and absolute magnitudes
- Compute luminosity from magnitude and distance
- Handle extinction corrections

#### 3. Celestial Mechanics
- Calculate angular separation between celestial objects
- Compute position angles
- Calculate proper motion in RA and Dec
- Determine radial velocities and corrections

#### 4. Time Conversions
- Convert between Julian Date (JD) and calendar dates
- Support Modified Julian Date (MJD)
- Handle different time scales (UTC, UT1, TT, TAI)
- Calculate sidereal time (LST, GST)

#### 5. Observational Astronomy
- Calculate airmass for given altitude
- Compute atmospheric extinction
- Determine rise, transit, and set times
- Calculate hour angle and altitude/azimuth

#### Plot results
- Generater the plots for uv_coverage, visibility_squared, closure_phase and visibility_multiwav

### Non-Functional Requirements

#### Performance
- Functions should execute efficiently for both scalar and array inputs
- Vectorized operations where possible
- Optimized for Julia's type system

#### Accuracy
- Numerical precision appropriate for professional astronomy (typically 1e-6 or better)
- Use established algorithms from astronomical references
- Include proper handling of edge cases

#### Usability
- At initiation prompt user for OIFITS file name
- Clear, descriptive function names
- Comprehensive documentation strings
- Consistent parameter ordering (e.g., RA before Dec)
- Support for common astronomical units (degrees, hours, radians)

## Module Structure

```julia
module MIRCX

# Exports
export ra_dec_to_galactic, galactic_to_ra_dec
export distance_modulus, apparent_to_absolute_mag
export angular_separation, position_angle
export jd_to_calendar, calendar_to_jd
export calculate_airmass, atmospheric_extinction

# Constants
const C_LIGHT = 299792458.0  # Speed of light in m/s
const PC_TO_M = 3.0856775814913673e16  # Parsec to meters
const AU_TO_M = 1.495978707e11  # AU to meters

# Implementation sections...
end
```

## Input/Output Specifications

### Coordinate Functions
- **Input**: Angles in degrees or radians (clearly specified)
- **Output**: Angles in same units as input or as specified
- **Convention**: RA in hours or degrees, Dec in degrees

### Time Functions
- **Input**: Calendar dates as DateTime objects or numeric JD
- **Output**: Numeric values for JD/MJD
- **Precision**: Double precision floating point

### Distance/Magnitude Functions
- **Input**: Magnitudes (dimensionless), distances in parsecs
- **Output**: Same units as appropriate for calculation
- **Range**: Handle extreme values gracefully

## Error Handling
- Validate input ranges (e.g., Dec must be -90 to +90 degrees)
- Throw descriptive errors for invalid inputs
- Handle NaN and Inf appropriately
- Document assumptions and limitations

## Testing Requirements
- Unit tests for each major function
- Test edge cases (poles, zero crossing, etc.)
- Validate against known astronomical values
- Performance benchmarks for critical functions

## Dependencies
- Base Julia (v1.6+)
- Dates module (standard library)
- Optional: AstroLib.jl for validation/comparison

## Documentation
- Module-level documentation describing purpose
- Function docstrings with:
  - Description of what the function does
  - Parameter descriptions with types and units
  - Return value description
  - Example usage
  - References to algorithms/papers if applicable

## Future Enhancements
- Support for additional coordinate systems
- Barycentric corrections
- Precession and nutation calculations
- Integration with FITS file handling
- Parallax and distance calculations
- Kepler orbit calculations

## References
- Meeus, J. "Astronomical Algorithms" (1998)
- The Astronomical Almanac
- USNO Circular 179
- IAU SOFA library documentation

## Version History
- v0.1.0 (planned): Initial implementation with core coordinate and time functions
