[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.17326808.svg)](https://doi.org/10.5281/zenodo.17326808)

# HF-TDOA-Frontiers-2025

Data and software repository to support:

**Cerwin, Stephen A., Jesse T. McMahan, Alexandros S. Papadopoulos, Gerard N. Piccini, Nathaniel A. Frissell, Kristina V. Collins, Aidain Montare, Paul Bilberry, Samuel Blackshear, and David R. Themens (2025). "HamSCI HF Multipath Propagation Mode Analysis Using Amateur Radios and Audio Waveforms Sensitive to Time Difference of Arrival." submitted to _Frontiers in Astronomy and Space Sciences - Space Physics_.**

**Archive Information:**
- **Zenodo DOI:** [10.5281/zenodo.17326808](https://doi.org/10.5281/zenodo.17326808)
- **GitHub Repository:** [https://github.com/HamSCI/HF-TDOA-Frontiers-2025](https://github.com/HamSCI/HF-TDOA-Frontiers-2025)
- **Archive Date:** January 2026
- **Software Version:** 1.0.0

# HF TDOA Analysis - Complete Manuscript Figures

This directory contains Jupyter notebooks, Python libraries, and supporting materials for generating all figures in the Frontiers manuscript on HF Time Difference of Arrival (TDOA) measurements from chirp sounder data.

## Overview

The analysis uses cross-correlation techniques to detect chirp signals in WAV recordings and extract TDOA measurements between different ionospheric propagation modes. These TDOAs are then converted to ionospheric layer heights using a spherical Earth virtual height model and compared with ionosonde measurements.

This repository contains:

- **15 figures** (Figures 1-15) with supporting materials for manuscript publication
- **11 Jupyter notebooks** for computational figure generation and analysis
- **Complete HF TDOA analysis library** (`hf_tdoa` package) with core functionality
- **Automated figure reproducibility testing** to verify computational results
- **Validation notebooks** for ionosonde data quality assurance

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
- Testing tools (pillow, scikit-image, nbconvert)

2. Install the `hf_tdoa` package in development mode:

```bash
pip install -e .
```

This makes the `hf_tdoa` library accessible from anywhere on your system, allowing notebooks in subdirectories to import it easily.

3. Launch Jupyter Lab:

```bash
jupyter lab
```

## File Structure

```
.
├── README.md                    # This file
├── environment.yml              # Conda environment specification
├── setup.py                     # Package installation script
├── test_figures.py              # Figure reproducibility test script
├── TEST_README.md               # Testing documentation
├── fig_01/                      # Figure 1: System diagram
│   ├── fig_01.jpg
│   └── fig_01.pptx
├── fig_02/                      # Figure 2: Receiver setup
│   ├── fig_02.jpg
│   └── fig_02.pptx
├── fig_03/                      # Figure 3: Spectrogram analysis
│   ├── fig_03.ipynb
│   └── fig_03.jpg
├── fig_04/                      # Figure 4: Chirp detection
│   ├── fig_04.ipynb
│   └── fig_04.jpg
├── fig_05/                      # Figure 5: Multi-path detection
│   ├── fig_05.ipynb
│   └── fig_05.jpg
├── fig_06/                      # Figure 6: Virtual height geometry
│   ├── fig_06_VirtualHeightGeometry.jpg
│   └── fig_06_VirtualHeightGeometry.pptx
├── fig_07/                      # Figure 7: TDOA spreadsheet
│   ├── Expected Mode TDOA -Figure-B.xlsx
│   └── fig_07_spreadsheet.jpg
├── fig_08/                      # Figure 8: Same path, different frequencies
│   ├── fig_08.ipynb
│   └── fig_08.jpg
├── fig_09/                      # Figure 9: Same frequency, different paths
│   ├── fig_09.ipynb
│   └── fig_09.jpg
├── fig_10/                      # Figure 10: Statistical comparison
│   ├── fig_10.jpg
│   └── fig_10.pptx
├── fig_11/                      # Figure 11: TDOA analysis (WA5FRF→N5DUP)
│   ├── fig_11.ipynb
│   └── fig_11.jpg
├── fig_12/                      # Figure 12: Dual-band comparison
│   ├── fig_12.ipynb
│   └── fig_12.jpg
├── fig_13/                      # Figure 13: Manual vs Automated Scatter
│   ├── fig_13.ipynb
│   └── fig_13.jpg
├── fig_14/                      # Figure 14: Multi-mode TDOA
│   ├── fig_14.ipynb
│   └── fig_14.jpg
├── fig_15/                      # Figure 15: Etalon mode analysis
│   ├── fig_15.ipynb
│   ├── fig_15.jpg
│   ├── fig_15_EtalonModes.jpg
│   └── fig_15_EtalonModes.pptx
├── fig_supporting/              # Supporting validation figures
│   ├── Austin_Ionosonde_Data_Validation.ipynb
│   ├── validation_virtual_heights.jpg
│   └── validation_critical_frequencies.jpg
├── data/                        # WAV recordings and CSV data
│   ├── TX_WA5FRF_EL09nn-RX_N5DUP_EM02ch-40m/
│   ├── TX_WA5FRF_EL09nn-RX_AB5YO_EL09mn-40m/
│   ├── TX_WA5FRF_EL09nn-RX_AB5YO_EL09mn-60m/
│   ├── TX_WA5FRF_EL09nn-RX_N6RFM_EM12jw-40m/
│   ├── templates/               # Reference chirp templates
│   │   └── N6RFM_10Hz_per_ms_template.wav
│   └── CSVs/                    # Ionosonde and manual TDOA data
├── submitted_figures/           # Baseline figures for testing
│   ├── fig_03.jpg
│   ├── fig_05.jpg
│   ├── fig_08.jpg
│   ├── fig_09.jpg
│   ├── fig_11.jpg
│   ├── fig_12.jpg
│   ├── fig_13.jpg
│   ├── fig_14.jpg
│   └── fig_15.jpg
└── hf_tdoa/                     # HF TDOA analysis package
    ├── __init__.py              # Package initialization
    ├── hf_tdoa_lib.py           # Core TDOA analysis library
    ├── calcSun.py               # Solar position calculations
    ├── eclipse_calc.py          # Eclipse obscuration calculations
    ├── solarContext.py          # Solar context and overlays
    ├── locator.py               # Grid square conversions
    ├── geopack.py               # Geodetic calculations
    ├── maps.py                  # Mapping utilities
    ├── gen_lib.py               # General utilities
    └── rayTracePaths.py         # Ray tracing functions
```

## Core Package: hf_tdoa

The `hf_tdoa` package provides comprehensive functionality for HF TDOA analysis. All functions and classes are accessible via `import hf_tdoa` thanks to the pythonic package structure.

### Installation as a Package

The `hf_tdoa` library is installed as an editable package using `pip install -e .` during setup. This approach:

- Makes the library importable from any location on your system
- Allows notebooks in subdirectories to easily access the library
- Enables immediate reflection of code changes without reinstallation
- Follows Python best practices for local development

After installation, you can import the library from anywhere:

```python
import hf_tdoa as tdoa
```

The library remains editable, so any changes you make to files in the `hf_tdoa/` directory are immediately available without reinstalling.

### Key Classes

#### `PathInfo`
Manages transmitter/receiver path information and calculations:
- Parses station callsigns and grid squares from filenames
- Calculates great circle distances and azimuths
- Computes path midpoints and layer heights
- Generates TDOA model coefficients (slope, intercept)

**Example:**
```python
import hf_tdoa as tdoa

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

#### Scatter Plot Comparison

- `align_and_resample_data()` - Align TDOA measurements with ionosonde data
- `plot_scatter_comparison()` - Create scatter plot comparing TDOA vs ionosonde
- `create_scatter_plot_figure()` - Create individual and combined scatter plots
- `create_manual_vs_automated_scatter_figure()` - Create 2x2 comparison figure

### Propagation Modes

Four propagation modes are pre-configured in `MODE_CONFIGS`:

| Mode | Description | Filter Limits (Hz) | Freq Range (Hz) | Line Style |
|------|-------------|-------------------|-----------------|------------|
| `3F2-1F2` | 3-hop F2 minus 1-hop F2 | 10-50 | 25-50 | `-.` (dash-dot) |
| `2F2-1F2` | 2-hop F2 minus 1-hop F2 | 10-50 | 11-20 | `--` (dashed) |
| `1F2-1E` | 1-hop F2 minus 1-hop E | 2.5-30 | 5-12 | `-.` (dash-dot) |
| `2F2-1E` | 2-hop F2 minus 1-hop E | 20-30 | 22-30 | `:` (dotted) |

Each mode includes optimized filter limits, search windows, and plotting parameters.

## Figures Overview

The repository contains materials for generating all manuscript figures:

### Static Figures (Pre-generated)

- **Figures 1-2**: System diagrams and receiver setup (PowerPoint/JPG)
- **Figure 6**: Virtual height geometry diagram (PowerPoint/JPG)
- **Figures 7, 10**: Spreadsheet analyses and comparisons (Excel/PowerPoint/JPG)

### Jupyter Notebook-Generated Figures

Each computational figure is organized in its own subdirectory with a Jupyter notebook that generates the figure.

#### [fig_03/fig_03.ipynb](fig_03/fig_03.ipynb)

Spectrogram analysis showing chirp signals in the time-frequency domain.

#### [fig_04/fig_04.ipynb](fig_04/fig_04.ipynb)

Chirp detection using cross-correlation template matching.

#### [fig_05/fig_05.ipynb](fig_05/fig_05.ipynb)

Multi-path propagation mode detection and analysis.

#### [fig_08/fig_08.ipynb](fig_08/fig_08.ipynb)

**Same path, different frequencies:**

- Compares TDOA measurements at different frequencies on the same propagation path
- WA5FRF→N6RFM path analysis
- Demonstrates frequency-dependent ionospheric behavior

#### [fig_09/fig_09.ipynb](fig_09/fig_09.ipynb)

**Same frequency, different paths:**

- Compares TDOA measurements at the same frequency on different propagation paths
- Multiple path analysis showing spatial variations
- Demonstrates path-dependent ionospheric characteristics

#### [fig_11/fig_11.ipynb](fig_11/fig_11.ipynb)

**TDOA analysis for WA5FRF→N5DUP 40m path:**

- Chirp detection via cross-correlation
- TDOA extraction for 2F2-1F2 mode
- Layer height comparison with Austin ionosonde
- Solar elevation and eclipse obscuration overlay
- Comparison with manual period analysis

#### [fig_12/fig_12.ipynb](fig_12/fig_12.ipynb)

**Dual-band comparison (40m and 60m):**

- Two-panel subplot comparing different frequency bands
- WA5FRF→AB5YO path analysis
- Demonstrates frequency-dependent propagation

#### [fig_13/fig_13.ipynb](fig_13/fig_13.ipynb)

**Manual vs Automated TDOA Scatter Comparison:**

- Compares manual and automated TDOA height measurements with Austin ionosonde hmF2
- 2x2 layout: Manual analysis (top) vs Automated analysis (bottom)
- Each row includes scatter plot and statistics table with correlation metrics
- Demonstrates agreement between manual period/autocorrelation and automated methods
- Validates automated TDOA extraction across multiple receivers and propagation modes

#### [fig_14/fig_14.ipynb](fig_14/fig_14.ipynb)

**Multi-mode TDOA analysis:**

- WA5FRF→N6RFM path with multiple propagation modes
- Processes 2F2-1F2, 1F2-1E, and 2F2-1E modes
- Compares different ionospheric layer interactions

#### [fig_15/fig_15.ipynb](fig_15/fig_15.ipynb)

**Etalon mode analysis:**

- Detection and analysis of etalon propagation modes
- Multi-hop reflection pattern identification

#### [fig_supporting/Austin_Ionosonde_Data_Validation.ipynb](fig_supporting/Austin_Ionosonde_Data_Validation.ipynb)

**Ionosonde data validation:**

- Validates Austin ionosonde measurements
- Generates supporting validation plots for virtual heights and critical frequencies

## Usage Example

### Basic Workflow

```python
import os
import hf_tdoa as tdoa

# Setup plotting
tdoa.setup_plotting_style()

# Define data paths
base_dir = 'data'
data_set = 'TX_WA5FRF_EL09nn-RX_N5DUP_EM02ch-40m'
template = os.path.join(base_dir, 'templates', 'N6RFM_10Hz_per_ms_template.wav')
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

# Plot with solar elevation, eclipse obscuration, and manual analysis overlays
tdoa.plot_hmf2(
    chirps, tdoa_dct,
    solar_lat=solar_lat,
    solar_lon=solar_lon,
    overlay_solar_elevation=True,
    overlay_eclipse=True,
    ionosonde_dct={'overlay_hmE': False},
    tdoa_csv_dct={
        'csv_path_period': 'data/CSVs/2024-04-08_TX_WA5FRF_EL09nn-RX_N5DUP_EM02ch-40m_manual_period_analysis.csv',
        'csv_path_autocorr': 'data/CSVs/2024-04-08_TX_WA5FRF_EL09nn-RX_N5DUP_EM02ch-40m_manual_autocorrelation_analysis.csv'
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

### Template Files

Reference chirp templates are located in [data/templates/](data/templates/) and are used for cross-correlation matching:

- `N6RFM_10Hz_per_ms_template.wav` - Standard 10 Hz/ms chirp template

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

## Figure Reproducibility Testing

The repository includes an automated testing script to verify figure reproducibility. This script executes all Jupyter notebooks and compares the generated figures against baseline figures.

### Running Tests

```bash
python test_figures.py
```

The test script:

1. Executes each notebook in its subdirectory
2. Compares generated JPG figures to baseline figures in [submitted_figures/](submitted_figures/)
3. Reports differences using pixel-level and structural similarity metrics

See [TEST_README.md](TEST_README.md) for detailed testing documentation.

## Troubleshooting

### Import Errors

If you encounter import errors with `hf_tdoa`, ensure the conda environment is properly activated:

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

For questions about the analysis code or data, please contact Nathaniel Frissell at nathaniel.frissell@scranton.edu.

## License

Software in this repository is released under the GNU General Public License version 3.

## Funding

NAF acknowledges the support of NASA grants 80NSSC23K1322 and 80NSSC25K7026, NSF grants AGS-2230345 and AGS-2045755, and grants from Amateur Radio Digital Communications (ARDC). DRT contributions to this study are supported through US Office of Naval Research PRISM Grant N00014-23-S-B001. KC contributions are supported through the National Science Foundation NSF AGS-2432824.

## Acknowledgments

We are grateful to the amateur radio community who voluntarily produced the HF radio observations used in this paper. Beat note analyses were developed by the University of Scranton Fall 2024 Digital Signal Processing (DSP) class led by Alexander Papadopoulos, Gerard Piccini, and Nathaniel Frissell and by correlation techniques developed by Tom McMahan. Claude.ai and ChatGPT large language models (LLMs) were used to assist in the implementation, refactoring, and documentation of python analysis code. LLMs were not used in the development of the methodology or scientific conjectures in this project, or the writing of this manuscript. Data and software for this paper are available from <https://doi.org/10.5281/zenodo.17326808>. Austin, TX ionosonde data was accessed through the University of Massachusetts Lowell Global Ionosphere Radio Observatory (GIRO) at <https://giro.uml.edu/>.

## Citation

If you use this software or data in your research, please cite:

```bibtex
@software{cerwin2025hftdoa,
  author = {Cerwin, Stephen A. and McMahan, Jesse T. and Papadopoulos, Alexandros S. and Piccini, Gerard N. and Frissell, Nathaniel A. and Collins, Kristina V. and Montare, Aidain and Bilberry, Paul and Blackshear, Samuel and Themens, David R.},
  title = {HF-TDOA-Frontiers-2025: HF TDOA Analysis Software},
  year = {2025},
  publisher = {Zenodo},
  version = {1.0.0},
  doi = {10.5281/zenodo.17326808},
  url = {https://doi.org/10.5281/zenodo.17326808}
}
```

And the associated paper:

```bibtex
@article{cerwin2025hamsci,
  author = {Cerwin, Stephen A. and McMahan, Jesse T. and Papadopoulos, Alexandros S. and Piccini, Gerard N. and Frissell, Nathaniel A. and Collins, Kristina V. and Montare, Aidain and Bilberry, Paul and Blackshear, Samuel and Themens, David R.},
  title = {HamSCI HF Multipath Propagation Mode Analysis Using Amateur Radios and Audio Waveforms Sensitive to Time Difference of Arrival},
  journal = {Frontiers in Astronomy and Space Sciences - Space Physics},
  year = {2025},
  note = {submitted}
}
```
