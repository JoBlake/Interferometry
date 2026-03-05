# OIFITS Image Reconstruction

This Julia script performs image reconstruction from optical interferometry data in OIFITS format using the OITOOLS package.

## Overview

The script reads OIFITS files containing visibility and closure phase measurements from optical interferometers and reconstructs an image of the observed astronomical object using iterative optimization techniques.

## Requirements

### Julia Packages

- **OITOOLS** (v0.8.0) - Optical interferometry data processing and image reconstruction
- **OIFITS** (v2.0.0) - OIFITS file format reading/writing
- **FITSIO** (v0.17.5) - FITS file I/O
- **AstroFITS** (v1.1.0) - Alternative FITS library used by OIFITS
- **Plots** - Visualization

### Installation

```julia
using Pkg
Pkg.add("OITOOLS")
Pkg.add("OIFITS")
Pkg.add("FITSIO")
Pkg.add("AstroFITS")
Pkg.add("Plots")
```

## Usage

### Basic Usage

1. Edit the `filename` variable in `reconstruct.jl` to point to your OIFITS file:
```julia
filename = "path/to/your/file.oifits"
```

2. Run the script:
```bash
julia reconstruct.jl
```

### Output

The script generates two output files:
- **reconstructed_image.png** - Visualization of the reconstructed image
- **reconstructed_image.fits** - FITS file containing the image data with WCS headers

## Configuration

### Image Parameters

Adjust these parameters in the script to customize the reconstruction:

```julia
npix = 128        # Image size in pixels (power of 2 recommended)
pixsize = 0.1     # Pixel scale in milliarcseconds (mas)
maxiter = 400     # Maximum number of iterations
```

### Regularization

The reconstruction uses regularization to constrain the solution. Current setting:

```julia
regularizers=[("centering", 1e-3)]
```

Available regularization options:
- `"centering"` - Keeps the image centered
- `"tv"` - Total variation (preserves sharp edges)
- `"entropy"` - Maximum entropy (produces smoother images)
- `"l2"` - L2 smoothness penalty

Example with multiple regularizers:
```julia
regularizers=[
    ("entropy", 1e-3),
    ("tv", 1e-4),
    ("centering", 1e-5)
]
```

## Compatibility Fixes

This script includes compatibility patches for known issues between OITOOLS v0.8.0 and OIFITS v2.0.0:

1. **Missing `OIFITS.load()` function** - Implemented as wrapper around `OIFITS.read()`
2. **Missing `OIFITS.select()` function** - Implemented to filter dataset tables
3. **OI_TARGET structure changes** - Added property accessor for `target_id` field
4. **Missing FOV column** - Uses `hack_revn=1` option for older OIFITS files

These patches are automatically applied when the script runs.

## Data Structure

The script expects OIFITS files containing:
- **OI_WAVELENGTH** table(s) - Spectral channel information
- **OI_TARGET** table - Target information
- **OI_VIS2** table(s) - Squared visibility measurements
- **OI_T3** table(s) - Closure phase/amplitude measurements

## Troubleshooting

### "Too many iterations" warning

The reconstruction didn't fully converge. Solutions:
- Increase `maxiter` (try 500-1000)
- Adjust regularization weights
- Try different regularizers

### Memory issues

Reduce the image size:
```julia
npix = 64  # Smaller image
```

### Poor image quality

Try adjusting:
- **Pixel scale** - Should match your target's expected size
- **Field of view** - `npix * pixsize` should encompass the object
- **Regularization** - Balance between smoothness and fitting the data

### File not found errors

Ensure the OIFITS file path is correct and uses forward slashes or escaped backslashes:
```julia
# Good
filename = "C:/Users/yourname/data/file.oifits"
# Also good
filename = "C:\\Users\\yourname\\data\\file.oifits"
```

## Technical Details

### Reconstruction Algorithm

The script uses:
1. **NFFT (Non-uniform FFT)** - Efficiently maps between image and UV plane
2. **Iterative optimization** - Minimizes chi-squared with regularization
3. **Regularization** - Constrains the solution to be physical

### Image Output

The FITS file includes:
- **BITPIX**: -32 (32-bit floating point)
- **WCS headers**: RA/Dec coordinates with TAN projection
- **Pixel scale**: Converted to degrees for standard compliance

## Example

```julia
# Reconstruct eps Aurigae from OIFITS data
filename = "2009-11-eps_Aur-avg5.oifits"

# Settings for this target
npix = 128
pixsize = 0.1  # mas
maxiter = 400

# The script will output:
# Number of V² measurements: 240
# Number of T3 measurements: 240
# Image size: 128 x 128 pixels
# Field of view: 12.8 mas
```

## References

- **OIFITS Standard**: https://www.aanda.org/articles/aa/pdf/2005/33/aa2811.pdf
- **OITOOLS.jl**: https://github.com/fabienbaron/OITOOLS.jl
- **OIFITS.jl**: https://github.com/emmt/OIFITS.jl

## Known Issues

- OITOOLS v0.8.0 and OIFITS v2.0.0 have compatibility issues that are patched in this script
- Some older OIFITS files may be missing required columns (handled via `hack_revn=1`)
- Convergence warnings are common and don't necessarily indicate bad results

## License

This script is provided as-is for astronomical data analysis.

## Author

Created for processing MIRCX and optical interferometry data.
