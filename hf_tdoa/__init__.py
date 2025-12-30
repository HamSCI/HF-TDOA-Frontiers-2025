"""
HF TDOA Package

This package provides tools for analyzing High Frequency (HF) Time Difference of Arrival (TDOA)
measurements from chirp sounder data.
"""

# Eclipse calculator submodules
from .eclipse_calc import calculate_obscuration
from .calcSun import calculate_solarAzEl
from . import locator
from . import maps
from . import gen_lib
from . import geopack
from . import solarContext
from . import rayTracePaths

# Main HF TDOA analysis functionality
from .hf_tdoa_lib import (
    # Classes
    PathInfo,
    TDOAData,

    # Configuration
    MODE_CONFIGS,

    # Setup and utility functions
    setup_plotting_style,
    obtain_wav_list,
    load_wav,

    # Signal processing functions
    filter,
    chirp_fft,
    find_chirps,
    find_TDOAs,
    find_max,

    # Plotting functions
    plot_chirp_fft,
    plot_TDOAs,
    plot_hmf2,
    plot_hmf2_subplot,
    plot_tdoa_hmf2_subplot,
    title_from_pfx,

    # Configuration builder
    build_tdoa_config,

    # Data loading and overlay functions
    load_ionosonde_data,
    load_tdoa_csv,
    load_tdoa_csv_generic,
    save_tdoa_csv,
    overlay_ionosonde,
    overlay_austin_ionosonde,
    overlay_tdoa_csv,

    # Scatter plot comparison functions
    align_and_resample_data,
    plot_resampling_validation,
    plot_scatter_comparison,
    create_scatter_plot_figure,
)

__all__ = [
    # Eclipse calculator exports
    'calculate_obscuration',
    'calculate_solarAzEl',
    'locator',
    'maps',
    'gen_lib',
    'geopack',
    'solarContext',
    'rayTracePaths',

    # HF TDOA exports
    'PathInfo',
    'TDOAData',
    'MODE_CONFIGS',
    'setup_plotting_style',
    'obtain_wav_list',
    'load_wav',
    'filter',
    'chirp_fft',
    'find_chirps',
    'find_TDOAs',
    'find_max',
    'plot_chirp_fft',
    'plot_TDOAs',
    'plot_hmf2',
    'plot_hmf2_subplot',
    'plot_tdoa_hmf2_subplot',
    'title_from_pfx',
    'build_tdoa_config',
    'load_ionosonde_data',
    'load_tdoa_csv',
    'load_tdoa_csv_generic',
    'save_tdoa_csv',
    'overlay_ionosonde',
    'overlay_austin_ionosonde',
    'overlay_tdoa_csv',
    'align_and_resample_data',
    'plot_resampling_validation',
    'plot_scatter_comparison',
    'create_scatter_plot_figure',
]
