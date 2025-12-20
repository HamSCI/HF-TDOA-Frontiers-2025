# HF TDOA Analysis - Figures 11-14

This directory contains Jupyter notebooks and Python libraries for analyzing High Frequency (HF) Time Difference of Arrival (TDOA) measurements from chirp sounder data, as presented in Figures 11-14 of the Frontiers manuscript.

## Overview

The analysis uses cross-correlation techniques to detect chirp signals in WAV recordings and extract TDOA measurements between different ionospheric propagation modes. These TDOAs are then converted to ionospheric layer heights using a spherical Earth virtual height model and compared with ionosonde measurements.

## Installation

### Prerequisites

- [Conda](https://docs.conda.io/en/latest/) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html)

### Setup Environment

1. Create and activate the conda environment:

```bash
conda env create -f environment.yml
conda activate hf-tdoa
```

This will install all required dependencies including:
- Python 3.12
- NumPy, SciPy, Pandas
- Matplotlib, Cartopy
- Astropy, GeographicLib
- JupyterLab
- tqdm, pydub

2. Launch Jupyter Lab:

```bash
jupyter lab
```

## File Structure

```
.
├── README.md                    # This file
├── environment.yml              # Conda environment specification
├── hf_tdoa_lib.py              # Core analysis library
├── fig_11.ipynb                # Figure 11 analysis notebook
├── fig_12.ipynb                # Figure 12 analysis notebook
├── fig_13_14.ipynb             # Figures 13-14 analysis notebook
├── data/                       # WAV recordings and CSV data
│   ├── TX_WA5FRF_EL09nn-RX_N5DUP_EM02ch-40m/
│   ├── TX_WA5FRF_EL09nn-RX_AB5YO_EL09so-40m/
│   ├── TX_WA5FRF_EL09nn-RX_AB5YO_EL09so-60m/
│   ├── TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m/
│   └── CSVs/                   # Ionosonde and manual TDOA data
├── templates/                  # Reference chirp templates
└── eclipse_calculator/         # Solar position and eclipse calculations
```

## Core Library: hf_tdoa_lib.py

The `hf_tdoa_lib.py` module provides comprehensive functionality for HF TDOA analysis:

### Key Classes

#### `PathInfo`
Manages transmitter/receiver path information and calculations:
- Parses station callsigns and grid squares from filenames
- Calculates great circle distances and azimuths
- Computes path midpoints and layer heights
- Generates TDOA model coefficients (slope, intercept)

**Example:**
```python
import hf_tdoa_lib as tdoa

# Create PathInfo from a prefix string
path_info = tdoa.PathInfo('TX_WA5FRF_EL09nn-RX_N5DUP_EM02ch-40m')

# Get coordinates
tx_lat, tx_lon = path_info.get_tx_latlon()
rx_lat, rx_lon = path_info.get_rx_latlon()

# Calculate TDOA model for 2F2-1F2 mode
slope, intercept = path_info.calculate_TDOA_model('2F2-1F2')
```

### Core Functions

#### Chirp Detection
- `obtain_wav_list()` - Get sorted list of WAV files
- `load_wav()` - Load and normalize WAV files
- `find_chirps()` - Cross-correlate WAV files with template to find chirp locations

#### TDOA Analysis
- `filter()` - Bandpass filter signals for specific propagation modes
- `chirp_fft()` - Compute FFT of chirp signals
- `find_TDOAs()` - Extract TDOA measurements from beat frequencies
- `build_tdoa_config()` - Build configuration with model coefficients

#### Visualization
- `plot_TDOAs()` - Plot TDOA measurements over time
- `plot_hmf2()` - Plot layer heights with ionosonde comparison
- `plot_hmf2_subplot()` - Multi-panel layer height plots
- `overlay_ionosonde()` - Add ionosonde data to plots
- `overlay_tdoa_csv()` - Add manual TDOA measurements to plots

### Propagation Modes

Three propagation modes are pre-configured in `MODE_CONFIGS`:

| Mode | Description | Filter Limits (Hz) | Freq Range (Hz) | Line Style |
|------|-------------|-------------------|-----------------|------------|
| `2F2-1F2` | 2-hop F2 minus 1-hop F2 | 10-50 | 11-20 | `--` (dashed) |
| `1F2-1E` | 1-hop F2 minus 1-hop E | 2.5-30 | 5-12 | `-.` (dash-dot) |
| `2F2-1E` | 2-hop F2 minus 1-hop E | 20-30 | 22-30 | `:` (dotted) |

Each mode includes optimized filter limits, search windows, and plotting parameters.

## Jupyter Notebooks

### fig_11.ipynb
Analyzes the WA5FRF→N5DUP 40m path and generates Figure 11:
- Chirp detection via cross-correlation
- TDOA extraction for 2F2-1F2 mode
- Layer height comparison with Austin ionosonde
- Solar elevation and eclipse obscuration overlay
- Comparison with manual period analysis

**Key Steps:**
1. Load WAV files and find chirps using template matching
2. Process TDOAs for the 2F2-1F2 propagation mode
3. Calculate TDOA model coefficients automatically
4. Generate layer height plot with multiple overlays

### fig_12.ipynb
Generates Figure 12 with dual-band analysis (40m and 60m):
- Two-panel subplot comparing different frequency bands
- WA5FRF→AB5YO path on both 40m and 60m
- Demonstrates multi-panel plotting capabilities

### fig_13_14.ipynb
Generates Figures 13 and 14:
- Analysis of WA5FRF→AB5YO and WA5FRF→N6RFM paths
- Multi-mode TDOA processing (2F2-1F2, 1F2-1E, 2F2-1E)
- Comparison of different propagation modes

## Usage Example

### Basic Workflow

```python
import os
import hf_tdoa_lib as tdoa

# Setup plotting
tdoa.setup_plotting_style()

# Define data paths
base_dir = 'data'
data_set = 'TX_WA5FRF_EL09nn-RX_N5DUP_EM02ch-40m'
template = os.path.join('templates', 'N6RFM_10Hz_per_ms_template.wav')
data_dir = os.path.join(base_dir, data_set)

# Get WAV file list
wavlist = tdoa.obtain_wav_list(data_dir)

# Find chirps (10 Hz/ms sweep rate)
chirps = tdoa.find_chirps(wavlist, template, sweep_rate=10)

# Process TDOA for 2F2-1F2 mode
chirps = tdoa.find_TDOAs(chirps, mode_string='2F2-1F2')

# Build configuration with model coefficients
tdoa_dct = tdoa.build_tdoa_config(chirps, mode_strings=['2F2-1F2'])

# Plot layer heights
tdoa.plot_hmf2(chirps, tdoa_dct)
```

### Advanced: Solar and Eclipse Overlays

```python
# Get path midpoint for solar calculations
path_info = chirps.attrs['path_info']
solar_lat, solar_lon = path_info.get_midpoint()

# Plot with solar elevation and eclipse obscuration
model_coeffs = tdoa_dct['2F2-1F2']['model_coeffs']

tdoa.plot_hmf2(
    chirps, tdoa_dct,
    solar_lat=solar_lat,
    solar_lon=solar_lon,
    overlay_solar_elevation=True,
    overlay_eclipse=True,
    ionosonde_dct={'overlay_hmE': False},
    tdoa_csv_dct={
        'csv_path': 'data/CSVs/2024-04-08_TX_WA5FRF_EL09nn-RX_N5DUP_EM02ch-40m_TDOA.csv',
        'model_coeffs': model_coeffs
    }
)
```

### Multi-Mode Analysis

```python
# Process multiple propagation modes
for mode in ['2F2-1F2', '1F2-1E', '2F2-1E']:
    chirps = tdoa.find_TDOAs(chirps, mode_string=mode)

# Build config for all modes
tdoa_dct = tdoa.build_tdoa_config(chirps)

# Plot all modes
tdoa.plot_hmf2(chirps, tdoa_dct, ylim=(75, 450))
```

## Data Format

### WAV Files
WAV files contain chirp sounder recordings with the naming convention:
```
YYYYMMDD.HHMM-<PREFIX>-<STATION_INFO>.wav
```

Example: `20240408.1413-TX_WA5FRF_EL09nn-RX_N5DUP_EM02ch-40m.wav`

### CSV Files
TDOA CSV files contain manual measurements with columns:
- `utc`: Timestamp (ISO format)
- `manualBeatNote_TDOA_ms`: Manual period analysis TDOA (ms)
- `autoCorrelation_TDOA_ms`: Auto-correlation TDOA (ms)

Ionosonde CSV files contain:
- `UTC`: Timestamp
- `hmF2`: F2 layer height (km)
- `hmE`: E layer height (km)

## Physical Model

The analysis uses a spherical Earth virtual height model to relate TDOA measurements to ionospheric layer heights:

**TDOA Calculation:**
```
TDOA [ms] = (Beat Frequency [Hz]) / (Sweep Rate [Hz/ms])
```

**Layer Height Model:**
For each propagation mode, the library calculates a linear relationship:
```
h_F2 = slope × TDOA + intercept
```

The model coefficients are automatically computed from the path geometry using equations 4-7 in the manuscript.

## Troubleshooting

### Import Errors
If you encounter import errors with `eclipse_calculator`, ensure the conda environment is properly activated:
```bash
conda activate hf-tdoa
```

### Missing Data
If WAV files are missing, check that you're in the correct directory and the `data/` folder contains the required datasets.

### Plot Not Displaying
If plots don't display in Jupyter, ensure you have:
```python
%matplotlib inline
```
at the top of your notebook.

## References

For details on the HF TDOA method and physical model, see the associated Frontiers manuscript.

## Contact

For questions about the analysis code or data, please refer to the Frontiers manuscript authors.

## License

See the main repository for license information.
