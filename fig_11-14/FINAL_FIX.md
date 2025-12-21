# Final Path Fix - Complete Solution

## Problem

After restructuring notebooks into subdirectories, two path issues emerged:

1. **Notebook paths**: References to `data/` and `templates/` from notebooks ✅ FIXED
2. **Library default paths**: The `overlay_ionosonde()` function in `hf_tdoa_lib.py` had a hardcoded relative path ✅ FIXED

## Solution

### 1. Notebook Paths (Previously Fixed)
Updated all notebooks to use `../data` and `../templates` instead of `data` and `templates`.

### 2. Library Default Path (Just Fixed)
Modified `hf_tdoa/hf_tdoa_lib.py` to automatically resolve the default ionosonde CSV path relative to the package installation location.

**Before:**
```python
def overlay_ionosonde(ax, csv_path='data/CSVs/2024-04-08_Austin_TX_Ionosonde_hmE_hmF2.csv', ...):
```

**After:**
```python
def overlay_ionosonde(ax, csv_path=None, ...):
    if csv_path is None:
        # Get the directory containing this package
        package_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        csv_path = os.path.join(package_dir, 'data', 'CSVs', '2024-04-08_Austin_TX_Ionosonde_hmE_hmF2.csv')
```

## How It Works

The `overlay_ionosonde()` function now:
1. Accepts `csv_path=None` as default
2. When `None`, dynamically calculates the path based on where the `hf_tdoa` package is installed
3. Uses `__file__` to find the package directory, then navigates to `../data/CSVs/...`

This works because:
- `__file__` = `/path/to/hf_tdoa/hf_tdoa_lib.py`
- `os.path.dirname(__file__)` = `/path/to/hf_tdoa/`
- `os.path.dirname(os.path.dirname(__file__))` = `/path/to/` (project root)
- `os.path.join(..., 'data', 'CSVs', 'file.csv')` = `/path/to/data/CSVs/file.csv`

## Benefits

1. **Works from any directory**: Notebooks can be in subdirectories, and the library still finds the data
2. **Portable**: Package can be installed anywhere on the system
3. **Backward compatible**: Notebooks can still override with explicit paths if needed
4. **No configuration needed**: Just works out of the box

## Testing

Run the test suite:
```bash
python test_figures.py
```

All three notebooks should now execute successfully and generate their figures!

## Changes Made

### Files Modified:
1. `hf_tdoa/hf_tdoa_lib.py` - Updated `overlay_ionosonde()` default path handling
2. `fig_11/fig_11.ipynb` - Updated paths to use `../data` and `../templates`
3. `fig_12/fig_12.ipynb` - Updated paths to use `../data` and `../templates`
4. `fig_13_14/fig_13_14.ipynb` - Updated paths to use `../data` and `../templates`

## Why This Approach?

Using `__file__` to dynamically resolve paths is a Python best practice for packages because:
- It works regardless of the current working directory
- It works whether the package is installed in development mode (`pip install -e .`) or normally
- It's more robust than relative paths
- It follows the "it just works" principle

The notebooks can now run from anywhere, and the library will always find the shared data directory!
