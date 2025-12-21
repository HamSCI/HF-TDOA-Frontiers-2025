# Notebook Path Updates

## Issue

After restructuring the project to place notebooks in subdirectories, the notebooks were failing because they used relative paths like `'data'` and `'templates'`, which no longer resolved correctly from the subdirectories.

## Solution

Updated all notebook cells to use parent directory references (`'../data'` and `'../templates'`) to correctly access shared resources.

## Changes Made

### fig_11/fig_11.ipynb
- Cell 3: Changed `base_dir = 'data'` to `base_dir = '../data'`
- Cell 3: Changed `template = os.path.join('templates', ...)` to `template = os.path.join('../templates', ...)`
- Cell 11: Changed `csv_path: 'data/CSVs/...'` to `csv_path: '../data/CSVs/...'`

### fig_12/fig_12.ipynb
- Cell 3: Changed `template = os.path.join('templates', ...)` to `template = os.path.join('../templates', ...)`
- Cell 4: Changed `base_dir = 'data'` to `base_dir = '../data'`
- Cell 14: Changed both CSV paths from `'data/CSVs/...'` to `'../data/CSVs/...'`

### fig_13_14/fig_13_14.ipynb
- Cell 3: Changed `base_dir = 'data'` to `base_dir = '../data'`
- Cell 3: Changed `template = os.path.join('templates', ...)` to `template = os.path.join('../templates', ...)`

## Directory Structure

```
fig_11-14/
├── data/           # Shared across all notebooks
├── templates/      # Shared across all notebooks
├── fig_11/
│   └── fig_11.ipynb      # Uses ../data and ../templates
├── fig_12/
│   └── fig_12.ipynb      # Uses ../data and ../templates
└── fig_13_14/
    └── fig_13_14.ipynb   # Uses ../data and ../templates
```

## Testing

Run the test suite to verify all notebooks execute correctly:
```bash
python test_figures.py
```

This will execute each notebook and verify that all figures are generated correctly.
