# Optical Interferometry Data Analysis Toolkit

A comprehensive suite of Julia tools for analyzing optical/infrared interferometry data in OIFITS format. Includes image reconstruction, diagnostics, parametric modeling, and analysis utilities for MIRCX and other interferometric data.

## Overview

This toolkit provides end-to-end analysis capabilities for optical interferometry observations:
- **Image reconstruction** from visibility and closure phase data
- **Data quality diagnostics** and coverage analysis
- **Parametric model fitting** as an alternative to reconstruction
- **Post-reconstruction analysis** tools for measuring source properties

## Available Tools

### Main Programs

#### `reconstruct.jl` - Image Reconstruction
Performs model-independent image reconstruction from OIFITS data using regularized optimization. Produces reconstructed images in both PNG and FITS formats.

**Features:**
- Auto-adjusts initial model size based on visibility amplitudes
- Supports multiple regularization methods (entropy, TV, centering)
- Interactive file selection
- Detailed convergence monitoring

**Usage:** `julia reconstruct.jl`

**Output:** `reconstructed_image.png`, `reconstructed_image.fits`, `initial_model.png`

#### `model_fit.jl` - Parametric Model Fitting
Alternative to image reconstruction that fits simple geometric models (uniform disk, elliptical disk) to interferometric data. Useful when you want to avoid reconstruction artifacts or when the source is expected to be simple.

**Features:**
- Fits uniform circular disk
- Fits elliptical disk with position angle
- Uses only V² data to avoid closure phase artifacts
- Provides chi-squared statistics for model comparison

**Usage:** `julia model_fit.jl`

### Diagnostic Tools

#### `quick_diagnostic.jl` - Fast Data Quality Check
Quick analysis of OIFITS data to understand source characteristics before reconstruction. Works even with files missing the OI_ARRAY table.

**Features:**
- Visibility amplitude statistics
- Source size estimation from mean V²
- Closure phase RMS analysis
- Automatic recommendations for initial model parameters

**Usage:** `julia quick_diagnostic.jl`

**Best for:** First look at new data, troubleshooting files

#### `check_uv_coverage.jl` - UV Coverage Analysis
Analyzes the spatial frequency coverage of your interferometric observations and recommends optimal image parameters.

**Features:**
- Calculates maximum angular resolution
- Recommends pixel size (Nyquist sampling)
- Determines appropriate field of view
- Generates UV coverage plots and baseline histograms

**Usage:** `julia check_uv_coverage.jl`

**Output:** `uv_coverage.png`, `baseline_distribution.png`

**Best for:** Planning reconstruction parameters

#### `diagnose_reconstruction.jl` - Comprehensive Reconstruction Diagnostics
In-depth analysis of data quality and forward modeling tests to diagnose reconstruction issues.

**Features:**
- Signal-to-noise ratio statistics
- Tests NFFT setup with forward models
- Evaluates uniform disk and point source models
- UV coverage uniformity analysis
- Provides specific troubleshooting recommendations

**Usage:** `julia diagnose_reconstruction.jl`

**Best for:** Understanding why reconstruction isn't working as expected

### Analysis Tools

#### `analyze_result.jl` - Reconstructed Image Analysis
Analyzes the properties of a reconstructed image, providing statistics and visualizations.

**Features:**
- Image statistics (min, max, mean, flux distribution)
- Flux concentration analysis
- Radial profile computation
- Contour plot generation

**Usage:** `julia analyze_result.jl`

**Output:** `reconstructed_contours.png`, `radial_profile.png`

**Best for:** Understanding the structure of your reconstructed image

#### `measure_ellipticity.jl` - Source Property Measurement
Measures physical properties of reconstructed sources using moment analysis.

**Features:**
- Centroid location
- Major/minor axis FWHM measurements
- Ellipticity calculation
- Position angle determination
- Half-light radius
- Multiple measurement methods (moments, contours)

**Usage:** `julia measure_ellipticity.jl`

**Best for:** Getting quantitative measurements from reconstructions

### Utility Tools

#### `check_versions.jl` - Package Version Checker
Displays installed versions of all required packages and their exported functions. Useful for debugging compatibility issues.

**Usage:** `julia check_versions.jl`

#### `inspect_oitarget.jl` - OIFITS Structure Inspector
Low-level inspection tool for examining the OI_TARGET table structure in OIFITS files. Useful for debugging OIFITS compatibility issues.

**Usage:** `julia inspect_oitarget.jl`

## Recommended Workflow

1. **First time with new data:**
   ```bash
   julia quick_diagnostic.jl          # Understand your data
   julia check_uv_coverage.jl         # Determine image parameters
   ```

2. **Image reconstruction:**
   ```bash
   julia reconstruct.jl               # Run reconstruction
   ```

3. **If reconstruction looks wrong:**
   ```bash
   julia diagnose_reconstruction.jl   # Diagnose issues
   # Adjust parameters in reconstruct.jl based on recommendations
   ```

4. **Analyze results:**
   ```bash
   julia analyze_result.jl            # View image properties
   julia measure_ellipticity.jl       # Get quantitative measurements
   ```

5. **Alternative approach (if reconstruction has artifacts):**
   ```bash
   julia model_fit.jl                 # Try parametric fitting instead
   ```

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

### Quick Start

Most tools are interactive and will prompt for the OIFITS file:

```bash
julia reconstruct.jl
# Press Enter to use default file or type path to your .oifits file
```

### Common Outputs

Different tools generate various outputs:
- **Images (PNG):** `reconstructed_image.png`, `initial_model.png`, `uv_coverage.png`, `baseline_distribution.png`, `reconstructed_contours.png`, `radial_profile.png`
- **FITS files:** `reconstructed_image.fits` (with WCS headers)
- **Console output:** Statistics, recommendations, and diagnostic information

## Configuration

### Image Reconstruction Parameters (`reconstruct.jl`)

Adjust these parameters to customize the reconstruction:

```julia
# Image grid
npix = 256        # Image size in pixels (power of 2 recommended)
fov = 100.0       # Field of view in mas
pixsize = fov / npix  # Pixel scale (auto-calculated)

# Initial model
initial_model = "gaussian"  # Options: "gaussian", "disk", "flat"

# Optimization
maxiter = 5000    # Maximum number of iterations
ftol = (1e-5, 1e-7)  # Function tolerance
xtol = (1e-5, 1e-7)  # Parameter tolerance
gtol = (1e-5, 1e-7)  # Gradient tolerance

# Data weighting
weights = [1.0, 1.0, 1.0]  # [V², T3amp, T3phi]
```

### Regularization

The reconstruction uses multiple regularizers to constrain the solution:

```julia
regularizers = [
    ("centering", 1e-4),   # Keep image centered
    ("entropy", 5e-5),     # Light smoothing for stability
    ("tv", 1e-5),          # Very light total variation for edges
]
```

**Available regularizers:**
- `"centering"` - Keeps the image centered (always recommended)
- `"entropy"` - Maximum entropy (smoothness)
- `"tv"` - Total variation (preserves sharp edges)
- `"tvsq"` - Squared total variation (smoother edges)
- `"l1l2"`, `"l1l2w"`, `"l1hyp"` - Sparsity penalties
- `"l2sq"` - L2 smoothness
- `"compactness"` - Favors compact sources
- `"radialvar"` - Minimizes radial variations

**Tuning tips:**
- **Too much regularization:** Reconstruction looks like initial model → reduce weights by 10×
- **Too little regularization:** Noisy, unstable result → increase weights by 2-5×
- **Balance:** Final chi² should be ~1.0 for good data fit

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
