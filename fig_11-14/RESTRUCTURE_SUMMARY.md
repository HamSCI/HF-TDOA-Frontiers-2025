# Project Restructuring Summary

## Changes Made

The project has been reorganized to improve modularity and testing capabilities.

## New Directory Structure

```
.
├── README.md                    # Main project documentation
├── environment.yml              # Conda environment specification
├── setup.py                     # Package installation script
├── MANIFEST.in                  # Package manifest for distribution
├── test_figures.py              # Figure reproducibility test script
├── TEST_README.md               # Testing documentation
├── fig_11/                      # Figure 11 analysis
│   └── fig_11.ipynb             # Figure 11 notebook (generates fig_11.jpg here)
├── fig_12/                      # Figure 12 analysis
│   └── fig_12.ipynb             # Figure 12 notebook (generates fig_12.jpg here)
├── fig_13_14/                   # Figures 13-14 analysis
│   └── fig_13_14.ipynb          # Figures 13-14 notebook (generates fig_13.jpg and fig_14.jpg here)
├── data/                        # WAV recordings and CSV data (shared)
├── templates/                   # Reference chirp templates (shared)
├── submitted_figures/           # Baseline figures for testing (shared)
└── hf_tdoa/                     # HF TDOA analysis package (shared)
    ├── __init__.py
    ├── hf_tdoa_lib.py
    └── ...
```

## Key Changes

### 1. Notebook Organization

**Before:** All notebooks in root directory
```
fig_11-14/
├── fig_11.ipynb
├── fig_12.ipynb
└── fig_13_14.ipynb
```

**After:** Each notebook in its own subdirectory
```
fig_11-14/
├── fig_11/
│   └── fig_11.ipynb
├── fig_12/
│   └── fig_12.ipynb
└── fig_13_14/
    └── fig_13_14.ipynb
```

**Benefits:**
- Better organization
- Figures are generated in the same directory as their notebook
- Easier to test individual notebooks
- Clear separation of outputs

### 2. Package Installation

**New:** `hf_tdoa` is now installed as an editable package

**Installation:**
```bash
conda env create -f environment.yml
conda activate hf-tdoa
pip install -e .
```

**Benefits:**
- Library can be imported from any location: `import hf_tdoa as tdoa`
- Notebooks in subdirectories can easily access the library
- Changes to library code are immediately available (no reinstall needed)
- Follows Python best practices for development

### 3. Test Script Updates

The `test_figures.py` script has been updated to:
- Look for notebooks in subdirectories (`fig_11/`, `fig_12/`, `fig_13_14/`)
- Find generated figures in those subdirectories
- Compare them against baselines in `submitted_figures/`

**Usage remains the same:**
```bash
python test_figures.py
```

## Migration Guide

If you have existing work, follow these steps:

### 1. Update Your Environment

```bash
# Pull the latest changes
git pull

# Update conda environment
conda env update -f environment.yml

# Install hf_tdoa package
pip install -e .
```

### 2. No Notebook Changes Needed

The notebooks already import `hf_tdoa` correctly and will work with the new structure once the package is installed.

### 3. Run Tests

Verify everything works:
```bash
python test_figures.py
```

## Shared Resources

The following directories are shared across all notebooks:

- **`data/`** - All notebooks access WAV files and CSV data from here
- **`templates/`** - Shared chirp templates
- **`hf_tdoa/`** - Core analysis library (installed as package)
- **`submitted_figures/`** - Baseline figures for comparison testing

## File Outputs

Each notebook generates its figures in its own subdirectory:

- `fig_11/fig_11.ipynb` → `fig_11/fig_11.jpg`
- `fig_12/fig_12.ipynb` → `fig_12/fig_12.jpg`
- `fig_13_14/fig_13_14.ipynb` → `fig_13_14/fig_13.jpg` and `fig_13_14/fig_14.jpg`

## Advantages of New Structure

1. **Modularity**: Each analysis is self-contained in its own directory
2. **Scalability**: Easy to add new figure notebooks without cluttering root
3. **Testing**: Test framework can validate each notebook independently
4. **Package Management**: Library follows Python best practices
5. **Collaboration**: Clearer organization for team development
6. **Reproducibility**: Better isolation of outputs for each analysis

## Technical Details

### setup.py

The `setup.py` file enables installation of `hf_tdoa` as a proper Python package:

- Defines package metadata (name, version, author, etc.)
- Lists all dependencies
- Specifies Python version requirement (3.12+)
- Enables editable installation with `pip install -e .`

### Editable Installation

The `-e` flag in `pip install -e .` creates an "editable" installation:

- Creates a link to the source code instead of copying it
- Changes to source files are immediately available
- No need to reinstall after modifying library code
- Perfect for development workflows

### Import Resolution

When you run `import hf_tdoa`, Python now finds it because:

1. `pip install -e .` registered the package with Python
2. Python's import system knows where to find it
3. Works from any directory, not just the project root

## Questions?

See the updated documentation:
- [README.md](README.md) - Main project documentation
- [TEST_README.md](TEST_README.md) - Testing documentation
