# Figure Reproducibility Tests

## Overview

The `test_figures.py` script automatically verifies that the Jupyter notebooks in this directory continue to produce the correct figure outputs as the code evolves.

## What It Does

1. **Executes each notebook**
   - Runs `fig_11.ipynb` → generates `fig_11.jpg`
   - Runs `fig_12.ipynb` → generates `fig_12.jpg`
   - Runs `fig_13_14.ipynb` → generates `fig_13.jpg` and `fig_14.jpg`

2. **Compares generated figures**
   - Compares newly generated JPG files against baseline figures in `submitted_figures/`
   - Uses multiple comparison metrics (pixel-level and perceptual)

3. **Reports differences**
   - Provides detailed metrics on any differences found
   - Clear pass/fail indicators with visual status symbols

## Installation

The required dependencies are included in `environment.yml`:
- numpy
- pillow
- scikit-image
- nbconvert

If you need to install them separately:
```bash
conda install pillow scikit-image nbconvert
```

## Usage

### Basic Usage

Run all tests from the `fig_11-14` directory:
```bash
python test_figures.py
```

### Advanced Options

```bash
# Specify custom directories
python test_figures.py --notebooks-dir . --submitted-dir submitted_figures

# Save detailed results to JSON
python test_figures.py --json-output test_results.json

# View help
python test_figures.py --help
```

## Understanding the Output

### Comparison Metrics

- **Identical**: Images match perfectly (pixel-for-pixel)
- **Mean Absolute Error (MAE)**: Average pixel difference across all pixels
- **RMSE**: Root mean squared error of pixel differences
- **Max Pixel Difference**: Maximum difference found in any single pixel
- **Pixels Different**: Percentage of pixels that differ between images
- **SSIM (Structural Similarity Index)**: Perceptual similarity metric (0-1 scale)

### SSIM Interpretation

- **>0.99**: Excellent - negligible differences
- **>0.95**: Good - minor differences
- **>0.90**: Fair - noticeable differences
- **<0.90**: Poor - significant differences

### Exit Codes

- **0**: All tests passed
- **1**: One or more tests failed or no tests were run

### Example Output

```
======================================================================
Figure Reproducibility Test
======================================================================
Notebooks directory: /path/to/fig_11-14
Baseline figures: /path/to/submitted_figures

──────────────────────────────────────────────────────────────────────
Testing: fig_11.ipynb
──────────────────────────────────────────────────────────────────────
  Executing notebook: fig_11.ipynb
  ✓ Notebook executed successfully

  fig_11.jpg:
    ✓ IDENTICAL - Images match perfectly!

──────────────────────────────────────────────────────────────────────
Testing: fig_12.ipynb
──────────────────────────────────────────────────────────────────────
  Executing notebook: fig_12.ipynb
  ✓ Notebook executed successfully

  fig_12.jpg:
    Shape: (2400, 3200, 3)
    Mean Absolute Error: 0.0234
    RMSE: 0.1523
    Max Pixel Difference: 12.0
    Pixels Different: 2.34%
    SSIM (Structural Similarity): 0.999234
    ✓ EXCELLENT - Negligible differences

======================================================================
Test Summary
======================================================================
Notebooks tested: 3
  Passed: 3
  Failed: 0

Figures compared: 4
  Identical: 3
  Different: 1
======================================================================

✓ ALL TESTS PASSED
```

## When to Be Concerned

- **SSIM < 0.95**: Investigate the differences - may indicate code changes affecting results
- **Execution failures**: Check that the notebooks can run in the current environment
- **Missing figures**: Verify that notebooks are saving figures with the expected filenames

## Expected Minor Differences

Minor differences (SSIM > 0.99) may occur due to:
- Floating-point precision variations across systems
- Random number generation (if not seeded)
- System font rendering differences
- Timestamp-based data if using current time

## Use Cases

### Development Workflow
Run tests before committing changes:
```bash
python test_figures.py
# If passed, safe to commit
git add .
git commit -m "Update analysis code"
```

### Continuous Integration
Add to GitHub Actions, GitLab CI, etc.:
```yaml
- name: Test figure reproducibility
  run: python test_figures.py --json-output test_results.json
```

### Regression Testing
Detect unintended changes to analysis output after:
- Updating dependencies
- Refactoring code
- Fixing bugs
- Adding new features

## JSON Output

Save detailed metrics for programmatic analysis:
```bash
python test_figures.py --json-output results.json
```

Example JSON structure:
```json
{
  "notebooks_tested": 3,
  "notebooks_passed": 3,
  "figures_compared": 4,
  "details": {
    "fig_11.ipynb": {
      "execution": {"success": true},
      "comparisons": {
        "fig_11.jpg": {
          "identical": true,
          "mae": 0.0,
          "ssim": 1.0
        }
      }
    }
  }
}
```

## Troubleshooting

### Import Errors
Ensure all required packages are installed:
```bash
conda install pillow scikit-image nbconvert
```

### Jupyter Not Found
Install Jupyter if needed:
```bash
conda install jupyter nbconvert
```

### Notebook Execution Timeout
Long-running notebooks may need increased timeout. Edit the `ExecutePreprocessor.timeout` value in `test_figures.py` (default: 600 seconds).

### Memory Issues
If notebooks consume too much memory, run tests individually by modifying the `notebook_outputs` dictionary in `test_figures.py`:
```python
self.notebook_outputs = {
    'fig_11.ipynb': ['fig_11.jpg'],
    # 'fig_12.ipynb': ['fig_12.jpg'],  # Commented out
    # 'fig_13_14.ipynb': ['fig_13.jpg', 'fig_14.jpg']
}
```

## Technical Details

- **Execution timeout**: Default 10 minutes per notebook (configurable in script)
- **Execution mode**: Notebooks run sequentially to avoid memory issues
- **In-place execution**: Notebooks are executed in-place (preserves outputs)
- **Image formats**: Currently supports JPG (easily extensible to PNG, etc.)
