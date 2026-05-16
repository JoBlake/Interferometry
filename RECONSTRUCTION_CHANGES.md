# Image Reconstruction Algorithm Changes

## Summary
Modified `reconstruct.jl` to produce reconstructions matching the π¹ Gruis target images (SQUEEZE/MIRA style reconstructions showing well-resolved stellar disks with asymmetric features).

## Data Characteristics (Pi_GRU.oifits)
- **Mean V²**: 0.115 → Well-resolved source (~15-20 mas diameter)
- **Closure Phases**: Mean |T3φ| = 73.84° → Strong asymmetry
- **Data Quality**: 909 V² measurements, 603 T3 measurements → Excellent sampling

## Changes Made

### 1. Field of View (FOV)
**Before**: 30 mas
**After**: 60 mas
**Reason**: Data analysis shows source is well-resolved (>12 mas). Mean V² = 0.115 requires larger FOV. Too small FOV was causing edge effects and artifacts.

### 2. Initial Model
**Before**: Gaussian with flattened peak
**After**: Uniform disk
**Reason**: Less biased starting point. The uniform disk lets the data drive the reconstruction without imposing a radial gradient structure.

### 3. Regularization Parameters (Defaults)

| Parameter | Before | After | Reason |
|-----------|---------|--------|---------|
| Centering | 1e-5 | 1e-6 | Minimal centering, just prevent drift |
| Entropy | 1e-6 | 1e-8 | Much weaker - let data dominate |
| TV | 1e-6 | 0 (disabled) | TV was causing ring artifacts |

**Key Change**: TV (Total Variation) disabled by default because it was creating ring-like artifacts in the reconstruction. For well-sampled data with good UV coverage, TV regularization can over-constrain and create spurious structures.

### 4. Algorithm Parameters
- **Max Iterations**: 10000 → 5000 (sufficient for well-sampled data)
- **Convergence Criteria**: Slightly relaxed from (1e-6, 1e-8) to (1e-5, 1e-7)

### 5. Regularizer Handling
- Now only includes non-zero regularizers in the list
- When TV = 0, it's completely removed from the regularizer list
- Displays "Active regularizers" count for transparency

## Why These Changes Work

1. **Larger FOV**: Properly contains the source without edge artifacts
2. **Weaker Regularization**: With 909 V² and 603 T3 measurements, the data is well-constrained. Strong regularization was over-smoothing and creating artifacts
3. **No TV**: TV regularization enforces edge-preserving smoothness, but was creating ring structures. For interferometric data with large closure phases (strong asymmetry), we want the data to define the structure, not the regularizer
4. **Simpler Initial Model**: Uniform disk is neutral and doesn't bias toward specific features

## Expected Results

The reconstruction should now show:
- ✓ Well-resolved disk (~15-20 mas diameter)
- ✓ Asymmetric brightness distribution
- ✓ Smooth but structured features
- ✓ Good contrast and dynamic range
- ✓ No ring artifacts
- ✓ Features consistent with SQUEEZE/MIRA reconstructions

## Usage

Run `reconstruct.jl` and press **Enter** three times to accept the optimized defaults:
- Centering: 1e-6
- Entropy: 1e-8
- TV: 0

For experimentation:
- **More smoothing**: Increase entropy to 1e-7 or 1e-6
- **Less smoothing**: Decrease entropy to 1e-9 or 1e-10
- **Add TV**: Use very weak values (1e-9 to 1e-8) if needed, but monitor for rings

## Technical Notes

The ring artifacts in previous reconstructions were a classic sign of:
1. TV regularization too strong for the data quality
2. FOV mismatch creating boundary effects
3. Over-regularization forcing the reconstruction away from the data

The mean V² of 0.115 is very low, indicating the source fills a large portion of the synthesized beam. This requires:
- Large FOV to capture extended structure
- Minimal regularization to let data define asymmetries
- Good initial model that doesn't impose structure
