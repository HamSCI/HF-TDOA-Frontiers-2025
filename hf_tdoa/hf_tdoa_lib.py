"""
HF TDOA Analysis Library
=========================

This library provides tools for analyzing High-Frequency (HF) radio Time Difference
of Arrival (TDOA) measurements from WAV file recordings of chirp signals.

Physical Concept
----------------
When an HF radio signal propagates via the ionosphere, it can take multiple paths
by reflecting off different ionospheric layers (E and F2) with different numbers of
"hops" (reflections between ground and ionosphere). Each path arrives at a slightly
different time, creating interference patterns known as "beats."

By measuring the Time Difference of Arrival (TDOA) between these paths and knowing
the transmitter-receiver geometry, we can infer the virtual height of ionospheric
layers. This technique provides remote sensing of the ionosphere using standard
amateur radio equipment.

Key Workflow
------------
1. **Data Loading**: Load WAV file recordings containing chirp signals
2. **Chirp Detection**: Cross-correlate with a template to find chirps in the recordings
3. **Signal Processing**: Filter and extract the beat frequency envelope
4. **TDOA Extraction**: Measure beat frequencies and convert to TDOA values
5. **Layer Height Calculation**: Apply the spherical Earth model to convert TDOA to layer heights
6. **Visualization**: Compare with ionosonde measurements and solar conditions

Mode Nomenclature
-----------------
Modes are labeled as "{n_hops}{Layer}-{m_hops}{Layer}":
  - "2F2-1F2": 2-hop F2 layer minus 1-hop F2 layer
  - "1F2-1E":  1-hop F2 layer minus 1-hop E layer
  - "2F2-1E":  2-hop F2 layer minus 1-hop E layer
  - "3F2-1F2": 3-hop F2 layer minus 1-hop F2 layer

Main Classes
------------
PathInfo : Parse transmitter/receiver geometry and calculate propagation parameters

Main Functions
--------------
Data Loading:
  - obtain_wav_list() : Get sorted list of WAV files from directory
  - load_wav() : Load WAV file into pandas DataFrame
  - load_ionosonde_data() : Load ionosonde CSV for validation
  - load_tdoa_csv() : Load previously-computed TDOA measurements

Signal Processing:
  - filter() : Bandpass filter signal and extract envelope
  - chirp_fft() : Calculate FFT of chirp signal
  - find_max() : Find maximum peak in frequency spectrum

Analysis:
  - find_chirps() : Locate chirps via cross-correlation with template
  - find_TDOAs() : Extract TDOA measurements from chirp signals
  - build_tdoa_config() : Build configuration for plotting multiple modes

Visualization:
  - setup_plotting_style() : Configure matplotlib defaults
  - plot_TDOAs() : Plot raw TDOA measurements vs time
  - plot_hmf2() : Plot derived layer heights vs time
  - plot_hmf2_subplot() : Multi-panel plots for multiple datasets
  - plot_tdoa_hmf2_subplot() : Combined TDOA and layer height plots

For detailed function documentation, see individual function docstrings.

Physical Constants
------------------
- Earth radius: 6378 km
- E layer height: 105 km (fixed)
- F2 layer height: Variable, typically 225-375 km
- Speed of light: 300,000 km/s
"""

import os
import glob
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.dates as mdates
import datetime
from scipy import signal
from scipy.io import wavfile
import pandas as pd
from tqdm import tqdm


# ============================================================================
# Configuration Constants and Global Parameters
# ============================================================================

# Mode-specific configuration parameters
# Each mode includes filter limits, search limits, and plotting parameters
MODE_CONFIGS = {
    '3F2-1F2': {
        'filter_limts': [10, 50],               # Bandpass filter limits [Hz]
        'search_limits': [-0.1, 0.1, 25, 50],   # (start_offset, end_offset, min_freq, max_freq)
        'linestyle': '-.',
        'linewidth': 3,
        'color': 'mediumblue',
        'marker': 'v'
    },
    '2F2-1F2': {
        'filter_limts': [10, 50],               # Bandpass filter limits [Hz]
        'search_limits': [-0.1, 0.1, 11, 20],   # (start_offset, end_offset, min_freq, max_freq)
        'linestyle': '--',
        'linewidth': 3,
        'color': 'tab:green',
        'marker': 'o'
    },
    '1F2-1E': {
        'filter_limts': [2.5, 30],              # Bandpass filter limits [Hz]
        'search_limits': [-0.1, 0.1, 5, 12],    # (start_offset, end_offset, min_freq, max_freq)
        'linestyle': '-.',
        'linewidth': 1.5,
        'color': 'tab:blue',
        'marker': '*'
    },
    '2F2-1E': {
        'filter_limts': [20, 30],               # Bandpass filter limits [Hz]
        'search_limits': [-0.1, 0.1, 22, 30],   # (start_offset, end_offset, min_freq, max_freq)
        'linestyle': ':',
        'linewidth': 2.5,
        'color': 'tab:orange',
        'marker': '^'
    }
}


# ============================================================================
# Core Data Classes
# ============================================================================

class TDOAData:
    """
    Class for loading and analyzing TDOA data from CSV files.

    This class loads automated TDOA analysis results from CSV files and provides
    methods for plotting and analyzing the data.

    Attributes:
    -----------
    csv_path : str
        Path to the CSV file
    df : pd.DataFrame
        DataFrame containing the TDOA data with UTC index
    metadata : dict
        Dictionary containing metadata parsed from CSV header
    """

    def __init__(self, csv_path):
        """
        Initialize TDOAData by loading a CSV file.

        Parameters:
        -----------
        csv_path : str
            Path to the TDOA CSV file
        """
        self.csv_path = csv_path
        self.metadata = {}
        self.df = self._load_csv()

    def _load_csv(self):
        """
        Load CSV file and parse metadata from header.

        Returns:
        --------
        df : pd.DataFrame
            DataFrame with UTC index and TDOA data
        """
        # Parse metadata from header comments
        with open(self.csv_path, 'r') as f:
            for line in f:
                if not line.startswith('#'):
                    break
                # Parse key-value pairs from header
                if ':' in line:
                    line = line.lstrip('#').strip()
                    if line.startswith('TX:'):
                        self.metadata['tx_call'] = line.split(':', 1)[1].strip()
                    elif line.startswith('TX Grid Square:'):
                        self.metadata['tx_grid'] = line.split(':', 1)[1].strip()
                    elif line.startswith('RX:'):
                        self.metadata['rx_call'] = line.split(':', 1)[1].strip()
                    elif line.startswith('RX Grid Square:'):
                        self.metadata['rx_grid'] = line.split(':', 1)[1].strip()
                    elif line.startswith('Frequency:'):
                        self.metadata['frequency'] = line.split(':', 1)[1].strip()
                    elif line.startswith('Ground Range:'):
                        range_str = line.split(':', 1)[1].strip()
                        self.metadata['ground_range_km'] = float(range_str.split()[0])
                    elif line.startswith('Propagation Mode:'):
                        self.metadata['mode'] = line.split(':', 1)[1].strip()
                    elif line.startswith('Model:'):
                        self.metadata['model_equation'] = line.split(':', 1)[1].strip()
                    elif line.startswith('Generated:'):
                        self.metadata['generated'] = line.split(':', 1)[1].strip()

        # Load CSV data
        df = pd.read_csv(self.csv_path, comment='#', parse_dates=['utc'])
        df = df.set_index('utc')

        return df

    def plot(self, savefig=None, ylim_tdoa=None, ylim_height=None, figsize=(16, 12)):
        """
        Create a publication-quality plot of TDOA measurements and layer heights.

        Creates a two-panel figure with:
        - Top panel: TDOA measurements vs time
        - Bottom panel: Layer heights vs time

        Parameters:
        -----------
        savefig : str, optional
            If provided, saves the figure to this file path. High-resolution JPEG recommended.
            If None (default), figure is displayed but not saved to disk.
        ylim_tdoa : tuple, optional
            Y-axis limits for TDOA plot (min, max). If None, uses automatic scaling.
        ylim_height : tuple, optional
            Y-axis limits for height plot (min, max). If None, uses automatic scaling.
        figsize : tuple, optional
            Figure size (width, height). Default is (16, 12).

        Returns:
        --------
        fig : matplotlib.figure.Figure
            The created figure object
        """
        fig, (ax1, ax2) = plt.subplots(2, 1, figsize=figsize)

        # Plot TDOA measurements
        ax1.plot(self.df.index, self.df['mean_tdoa_ms'],
                marker='o', linestyle='-', linewidth=2, markersize=6,
                color='tab:blue', label=f"Mode: {self.metadata.get('mode', 'N/A')}")
        ax1.set_ylabel('TDOA [ms]')
        ax1.set_xlabel('Time UTC')
        ax1.legend(loc='best')
        ax1.grid(True, linestyle=':', alpha=0.7)
        if ylim_tdoa is not None:
            ax1.set_ylim(ylim_tdoa)

        # Format datetime axis
        _format_datetime_axis(ax1)

        # Add title with path information
        if 'tx_call' in self.metadata and 'rx_call' in self.metadata:
            title_left = f"TX: {self.metadata['tx_call']} ({self.metadata.get('tx_grid', 'N/A')})\n"
            title_left += f"RX: {self.metadata['rx_call']} ({self.metadata.get('rx_grid', 'N/A')})"
            ax1.set_title(title_left, loc='left', fontsize=14)

            if 'ground_range_km' in self.metadata and 'frequency' in self.metadata:
                title_right = f"Ground Range: {self.metadata['ground_range_km']:.1f} km\n"
                title_right += f"Band: {self.metadata['frequency']}"
                ax1.set_title(title_right, loc='right', fontsize=14)

        # Add center title with mode
        if 'mode' in self.metadata:
            # Get date from first timestamp
            date_str = self.df.index[0].strftime('%Y %b %d')
            ax1.set_title(f"TDOA\n{date_str}", loc='center', fontsize=18, fontweight='bold')

        # Plot layer heights
        ax2.plot(self.df.index, self.df['mean_tdoa_hgt_km'],
                marker='o', linestyle='-', linewidth=2, markersize=6,
                color='tab:green', label=f"Mode: {self.metadata.get('mode', 'N/A')}")
        ax2.set_ylabel('Layer Height [km]')
        ax2.set_xlabel('Time UTC')
        ax2.legend(loc='best')
        ax2.grid(True, linestyle=':', alpha=0.7)
        if ylim_height is not None:
            ax2.set_ylim(ylim_height)

        # Format datetime axis
        _format_datetime_axis(ax2)

        # Add title with path information
        if 'tx_call' in self.metadata and 'rx_call' in self.metadata:
            ax2.set_title(title_left, loc='left', fontsize=14)
            if 'ground_range_km' in self.metadata and 'frequency' in self.metadata:
                ax2.set_title(title_right, loc='right', fontsize=14)

        # Add center title
        if 'mode' in self.metadata:
            ax2.set_title(f"Ionospheric Layer Height\n{date_str}", loc='center',
                         fontsize=18, fontweight='bold')

        # Add subplot labels
        ax1.text(-0.08, 1.075, '(a)', transform=ax1.transAxes,
                fontsize=24, fontweight='bold', va='top', ha='left')
        ax2.text(-0.08, 1.075, '(b)', transform=ax2.transAxes,
                fontsize=24, fontweight='bold', va='top', ha='left')

        # Rotate x-axis labels
        for ax in [ax1, ax2]:
            for tick_label in ax.get_xticklabels():
                tick_label.set_rotation(45)
                tick_label.set_horizontalalignment('right')

        # Apply tight layout to prevent overlapping elements
        fig.tight_layout()


        # Save figure if filename provided
        if savefig is not None:
            fig.savefig(savefig, dpi=300, bbox_inches='tight', format='jpeg',
                       pil_kwargs={'quality': 95})
            print(f"Figure saved to: {savefig}")
            plt.close(fig)  # Close the figure to free memory
            return None

        return fig

    def __repr__(self):
        """String representation of TDOAData."""
        mode = self.metadata.get('mode', 'N/A')
        n_points = len(self.df)
        return (f"TDOAData(mode='{mode}', n_points={n_points}, "
                f"csv_path='{os.path.basename(self.csv_path)}')")


class PathInfo:
    """
    Parses and stores transmitter/receiver path information from a prefix string.

    The prefix string is expected to follow the format:
    <anything>-TX_CALL-TX_GRID-<anything>-RX_CALL-RX_GRID-<BAND>m

    Attributes:
    -----------
    pfx : str
        Original prefix string
    tx_call : str
        Transmitter callsign
    tx_grid : str
        Transmitter grid square
    rx_call : str
        Receiver callsign
    rx_grid : str
        Receiver grid square
    band : int
        Band in meters (e.g., 20, 40, 80)
    band_str : str
        Band frequency string (e.g., '14 MHz', '7 MHz')
    Re : float
        Earth radius in km (6378.0)
    E_layer_height : float
        Fixed E layer height in km (105.0)
    """

    def __init__(self, pfx):
        """
        Initialize PathInfo by parsing the prefix string.

        Parameters:
        -----------
        pfx : str
            Prefix string containing TX/RX station information
        """
        self.pfx = pfx
        self.Re = 6378.0  # Earth radius in km
        self.E_layer_height = 105.0  # Fixed E layer height in km
        self._parse_prefix()

    def _parse_prefix(self):
        """Parse the prefix string and extract TX/RX information."""
        pfx_parts = self.pfx.replace('_', '-').split('-')

        self.tx_call = pfx_parts[1]
        self.tx_grid = pfx_parts[2]
        self.rx_call = pfx_parts[4]
        self.rx_grid = pfx_parts[5]

        self.band = int(pfx_parts[6].replace('m', ''))

        self.band_str = self._get_band_str()

    def _get_band_str(self):
        """Convert band (in meters) to frequency string."""
        band_dct = {
            80: '3.5 MHz',
            60: '5 MHz',
            40: '7 MHz',
            20: '14 MHz',
            15: '21 MHz',
            10: '28 MHz'
        }
        return band_dct.get(self.band, f'{self.band} m')

    def _import_locator(self):
        """
        Import locator module directly to avoid astropy dependency issues.

        Locator imports geopack, so we need to load geopack first and inject it
        into the locator module's namespace to handle the relative import.

        Returns:
        --------
        module
            The locator module
        """
        import importlib.util
        import os
        import sys
        import types

        # Try to find locator.py in hf_tdoa package directory
        try:
            # Get the directory of this file (hf_tdoa package)
            this_dir = os.path.dirname(os.path.abspath(__file__))
            locator_path = os.path.join(this_dir, 'locator.py')
            geopack_path = os.path.join(this_dir, 'geopack.py')

            if not os.path.exists(locator_path):
                # Fallback: try relative import
                try:
                    from . import locator
                    return locator
                except ImportError:
                    from hf_tdoa import locator
                    return locator

            # First load geopack into a fake package structure
            geopack_spec = importlib.util.spec_from_file_location("hf_tdoa.geopack", geopack_path)
            geopack = importlib.util.module_from_spec(geopack_spec)

            # Create a fake hf_tdoa package if it doesn't exist
            if 'hf_tdoa' not in sys.modules:
                hf_tdoa_pkg = types.ModuleType('hf_tdoa')
                hf_tdoa_pkg.__path__ = [this_dir]
                hf_tdoa_pkg.__package__ = 'hf_tdoa'
                sys.modules['hf_tdoa'] = hf_tdoa_pkg

            # Add geopack to sys.modules with package name
            sys.modules['hf_tdoa.geopack'] = geopack
            geopack_spec.loader.exec_module(geopack)

            # Now load locator with package context - it can now do "from . import geopack"
            spec = importlib.util.spec_from_file_location("hf_tdoa.locator", locator_path)
            locator = importlib.util.module_from_spec(spec)
            sys.modules['hf_tdoa.locator'] = locator
            spec.loader.exec_module(locator)
            return locator
        except Exception:
            # Final fallback: try regular import
            try:
                from . import locator
                return locator
            except ImportError:
                from hf_tdoa import locator
                return locator

    def _import_geopack(self):
        """
        Import geopack module directly to avoid astropy dependency issues.

        Returns:
        --------
        module
            The geopack module
        """
        import importlib.util
        import os

        # Try to find geopack.py in hf_tdoa package
        try:
            # Get the directory of this file (hf_tdoa package)
            this_dir = os.path.dirname(os.path.abspath(__file__))
            geopack_path = os.path.join(this_dir, 'geopack.py')

            if not os.path.exists(geopack_path):
                # Fallback: try relative import
                try:
                    from . import geopack
                    return geopack
                except ImportError:
                    from hf_tdoa import geopack
                    return geopack

            # Load module directly to avoid __init__.py astropy import
            spec = importlib.util.spec_from_file_location("geopack", geopack_path)
            geopack = importlib.util.module_from_spec(spec)
            spec.loader.exec_module(geopack)
            return geopack
        except Exception:
            # Final fallback: try regular import
            try:
                from . import geopack
                return geopack
            except ImportError:
                from hf_tdoa import geopack
                return geopack

    def get_tx_latlon(self):
        """
        Get transmitter coordinates from grid square.

        Returns:
        --------
        tuple
            (latitude, longitude) in degrees
        """
        locator = self._import_locator()
        lat, lon = locator.gridsquare2latlon(self.tx_grid, position='center')
        return (_ensure_scalar(lat), _ensure_scalar(lon))

    def get_rx_latlon(self):
        """
        Get receiver coordinates from grid square.

        Returns:
        --------
        tuple
            (latitude, longitude) in degrees
        """
        locator = self._import_locator()
        lat, lon = locator.gridsquare2latlon(self.rx_grid, position='center')
        return (_ensure_scalar(lat), _ensure_scalar(lon))

    def get_midpoint(self):
        """
        Calculate the great circle midpoint between TX and RX stations.

        Returns:
        --------
        tuple
            (midpoint_lat, midpoint_lon) in degrees
        """
        locator = self._import_locator()
        return locator.gridsquare_midpoint(self.tx_grid, self.rx_grid)

    def get_path_azimuth(self):
        """
        Calculate the azimuth from TX to RX along the great circle path.

        Returns:
        --------
        float
            Azimuth in degrees (0-360, where 0 is North)
        """
        geopack = self._import_geopack()
        tx_lat, tx_lon = self.get_tx_latlon()
        rx_lat, rx_lon = self.get_rx_latlon()
        azm = geopack.greatCircleAzm(tx_lat, tx_lon, rx_lat, rx_lon)
        return azm

    def get_range_km(self):
        """
        Calculate the great circle distance between TX and RX stations.

        Returns:
        --------
        float
            Distance in kilometers
        """
        locator = self._import_locator()
        # Use the grid_range_km function from locator module
        dist_km = locator.grid_range_km(self.tx_grid, self.rx_grid, Re=6371.0)
        return dist_km

    def calculate_path_length(self, n_hops, layer_height):
        """
        Calculate the total path length for n-hop propagation at a given layer height.

        Uses the spherical Earth virtual height model from Section 2.4 of the manuscript.
        Equations 4-7:
        - P1 = 2√[2r(r + h)(1 - cos(D/2r)) + h²]
        - P2 = 4√[2r(r + h)(1 - cos(D/4r)) + h²]
        - P3 = 6√[2r(r + h)(1 - cos(D/6r)) + h²]
        - P4 = 8√[2r(r + h)(1 - cos(D/8r)) + h²]

        Parameters:
        -----------
        n_hops : int
            Number of hops (1, 2, 3, or 4)
        layer_height : float
            Virtual layer height in km

        Returns:
        --------
        float
            Total path length in km
        """
        D = self.get_range_km()  # Ground distance in km
        r = self.Re
        h = layer_height

        # Calculate path length using spherical Earth virtual height formula
        path_length = (2 * n_hops) * np.sqrt(
            2 * r * (r + h) * (1 - np.cos(D / (2 * n_hops * r))) + h**2
        )

        return path_length

    def calculate_TOF(self, n_hops, layer_height):
        """
        Calculate Time of Flight for n-hop propagation at a given layer height.

        Parameters:
        -----------
        n_hops : int
            Number of hops (1, 2, 3, or 4)
        layer_height : float
            Virtual layer height in km

        Returns:
        --------
        float
            Time of Flight in milliseconds
        """
        c = 300000.0  # Speed of light in km/s
        path_length = self.calculate_path_length(n_hops, layer_height)
        tof_ms = (path_length / c) * 1000  # Convert to milliseconds
        return tof_ms

    def calculate_TDOA(self, mode1, mode2, layer_height_F2=None, layer_height_E=None):
        """
        Calculate TDOA between two propagation modes.

        Mode nomenclature: "{n_hops}{Layer}-{m_hops}{Layer}"
        Examples:
        - "2F2-1F2": 2-hop F2 minus 1-hop F2
        - "1F2-1E": 1-hop F2 minus 1-hop E
        - "2F2-1E": 2-hop F2 minus 1-hop E

        Parameters:
        -----------
        mode1 : str
            First mode (e.g., "2F2", "1F2", "1E")
        mode2 : str
            Second mode (e.g., "1F2", "1E")
        layer_height_F2 : float, optional
            F2 layer height in km. If None, uses a default of 280 km.
        layer_height_E : float, optional
            E layer height in km. If None, uses self.E_layer_height (105 km).

        Returns:
        --------
        float
            TDOA in milliseconds (TOF_mode1 - TOF_mode2)
        """
        if layer_height_F2 is None:
            layer_height_F2 = 280.0  # Default F2 layer height
        if layer_height_E is None:
            layer_height_E = self.E_layer_height

        # Parse mode strings
        def parse_mode(mode):
            """Parse mode string like '2F2' into (n_hops=2, layer='F2')"""
            # Extract number of hops (first digit)
            n_hops = int(mode[0])
            # Extract layer (rest of string)
            layer = mode[1:]
            return n_hops, layer

        n_hops1, layer1 = parse_mode(mode1)
        n_hops2, layer2 = parse_mode(mode2)

        # Get appropriate layer heights
        h1 = layer_height_F2 if layer1 == 'F2' else layer_height_E
        h2 = layer_height_F2 if layer2 == 'F2' else layer_height_E

        # Calculate TOFs
        tof1 = self.calculate_TOF(n_hops1, h1)
        tof2 = self.calculate_TOF(n_hops2, h2)

        # Return TDOA
        return tof1 - tof2

    def calculate_TDOA_model(self, mode_string, h_min=225, h_max=375, n_points=100):
        """
        Calculate a linear TDOA model (slope, intercept) for a given propagation mode.

        The model fits a line: layer_height = slope * TDOA + intercept
        over the F2 layer height range.

        Parameters:
        -----------
        mode_string : str
            Mode string like "1F2-2F2", "1F2-1E", "2F2-1E"
        h_min : float, optional
            Minimum F2 layer height in km (default: 225)
        h_max : float, optional
            Maximum F2 layer height in km (default: 375)
        n_points : int, optional
            Number of points for linear fit (default: 100)

        Returns:
        --------
        tuple
            (slope, intercept) for the linear model
        """
        # Parse the mode string (e.g., "1F2-2F2" -> mode1="1F2", mode2="2F2")
        modes = mode_string.split('-')
        if len(modes) != 2:
            raise ValueError(f"Invalid mode string: {mode_string}. Expected format: 'mode1-mode2'")

        mode1, mode2 = modes

        # Generate array of F2 layer heights
        heights = np.linspace(h_min, h_max, n_points)

        # Calculate TDOAs for each height
        tdoas = np.array([self.calculate_TDOA(mode1, mode2, layer_height_F2=h) for h in heights])

        # Fit a linear model: height = slope * TDOA + intercept
        # Use polyfit with degree 1
        slope, intercept = np.polyfit(tdoas, heights, 1)

        return slope, intercept

    def __repr__(self):
        """String representation of PathInfo."""
        range_km = self.get_range_km()
        return (f"PathInfo(TX: {self.tx_call} ({self.tx_grid}), "
                f"RX: {self.rx_call} ({self.rx_grid}), "
                f"Range: {range_km:.1f} km, Band: {self.band_str})")


# ============================================================================
# Data Loading Functions
# ============================================================================

def obtain_wav_list(directory):
    """
    Obtains a sorted list of WAV files in a directory.

    Arguments:
    directory : str
        Path to directory containing WAV files.

    Returns:
    wav_list : list
        Sorted list of WAV file paths.
    """
    wav_list = glob.glob(os.path.join(directory,'*.wav'))
    wav_list.sort()
    return wav_list


def load_wav(fname, normalize=True):
    """
    Loads a WAV file and returns it as a pandas DataFrame.

    Arguments:
    fname : str
        Path to WAV file.
    normalize : bool, optional
        If True, normalizes the signal amplitude. Default is True.

    Returns:
    df : pd.DataFrame
        DataFrame with 'time' index and 'x' column containing signal data.
        Includes 'fs' (sampling frequency) in df.attrs.
    fs : float
        Sampling frequency in Hz.
    """
    # Load WAV file.
    fs, x0 = wavfile.read(fname)

    if normalize:
        x0 = x0/float(np.max(np.abs(x0)))

    # Compute time vector.
    N  = len(x0)             # Number of samples in signal
    k  = np.arange(len(x0))  # Integer time vector
    Ts = 1/fs                # Sampling Period
    t  = k*Ts

    df = pd.DataFrame({'time':t, 'x':x0})
    df = df.set_index('time')
    df.attrs['fs'] = fs

    return df, fs


def load_ionosonde_data(csv_path=None):
    """
    Loads ionosonde data from a CSV file.

    The CSV file should have columns: UTC, foF2, foF1, foE, hmF2, hmF1, hmE
    Header lines starting with # are automatically skipped.

    Arguments:
    csv_path : str, optional
        Path to the ionosonde CSV file.
        If None, uses default Austin TX ionosonde data.

    Returns:
    ionosonde_df : pd.DataFrame
        DataFrame containing ionosonde data with datetime index named 'UTC'.
        Columns include: foF2, foF1, foE, hmF2, hmF1, hmE (as available in the file).
    """
    # If no path provided, use default relative to package root
    if csv_path is None:
        # Get the directory containing this file (hf_tdoa package directory)
        package_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        csv_path = os.path.join(package_dir, 'data', 'CSVs', '2024-04-08_AU930_AustinTX_Ionosonde_ManualScaled.csv')

    # Load CSV file, skipping comment lines starting with #
    ionosonde_df = pd.read_csv(csv_path, comment='#', parse_dates=['UTC'])
    ionosonde_df = ionosonde_df.set_index('UTC')

    return ionosonde_df


def load_tdoa_csv(csv_path):
    """
    Loads TDOA data from a CSV file containing manual analysis results.

    Expected CSV format:
    - utc: timestamp in ISO 8601 format (e.g., '2024-04-08 14:13')
    - tdoa_ms: Manual analysis TDOA values in milliseconds (optional)
    - tdoa_hgt_km: Manual analysis layer heights already in kilometers (required)

    The CSV file should have pre-calculated layer heights in the tdoa_hgt_km column.
    This function is used for loading verified manual beatnote analysis and manual
    autocorrelation analysis CSV files.

    Arguments:
    csv_path : str
        Path to the CSV file containing TDOA data.

    Returns:
    tdoa_df : pd.DataFrame
        DataFrame containing:
        - UTC datetime index
        - manualBeatNote_TDOA_ms: Original TDOA values in ms (if available)
        - manualBeatNote_height_km: Layer heights from manual analysis
    """
    # Load CSV with ISO 8601 datetime parsing, skipping comment lines starting with #
    tdoa_df = pd.read_csv(
        csv_path,
        parse_dates=['utc'],
        date_format='%Y-%m-%d %H:%M',
        comment='#'
    )
    tdoa_df = tdoa_df.set_index('utc')

    # Check that we have the required column
    if 'tdoa_hgt_km' not in tdoa_df.columns:
        raise ValueError(f"CSV file {csv_path} must contain 'tdoa_hgt_km' column with pre-calculated layer heights")

    # Use the pre-calculated heights
    tdoa_df['manualBeatNote_height_km'] = tdoa_df['tdoa_hgt_km']

    # Keep the TDOA values if they exist
    if 'tdoa_ms' in tdoa_df.columns:
        tdoa_df['manualBeatNote_TDOA_ms'] = tdoa_df['tdoa_ms']

    return tdoa_df


def load_tdoa_csv_generic(csv_path):
    """
    Loads TDOA data from either manual or automated analysis CSV files.

    This is a generic loader that handles both:
    - Manual analysis CSVs with 'tdoa_hgt_km' column
    - Automated analysis CSVs with 'mean_tdoa_hgt_km' column

    Arguments:
    csv_path : str
        Path to the CSV file containing TDOA data.

    Returns:
    tdoa_df : pd.DataFrame
        DataFrame containing:
        - UTC datetime index
        - height_km: Layer heights (from either manual or automated analysis)
    """
    # Load CSV with ISO 8601 datetime parsing, skipping comment lines starting with #
    tdoa_df = pd.read_csv(
        csv_path,
        parse_dates=['utc'],
        date_format='%Y-%m-%d %H:%M',
        comment='#'
    )
    tdoa_df = tdoa_df.set_index('utc')

    # Determine which type of CSV this is and extract height column
    if 'tdoa_hgt_km' in tdoa_df.columns:
        # Manual analysis CSV
        tdoa_df['height_km'] = tdoa_df['tdoa_hgt_km']
    elif 'mean_tdoa_hgt_km' in tdoa_df.columns:
        # Automated analysis CSV
        tdoa_df['height_km'] = tdoa_df['mean_tdoa_hgt_km']
    else:
        raise ValueError(
            f"CSV file {csv_path} must contain either 'tdoa_hgt_km' (manual) "
            f"or 'mean_tdoa_hgt_km' (automated) column"
        )

    return tdoa_df[['height_km']]


def save_tdoa_csv(chirps, mode_string, tdoa_config, output_dir, date_str=None):
    """
    Saves TDOA measurements and derived layer heights to a CSV file.

    This function creates CSV files with automated TDOA analysis results that can be
    used to recreate figures. The CSV includes header comments with metadata and
    model coefficients.

    Arguments:
    chirps : pd.DataFrame
        DataFrame containing chirp data and TDOA measurements from find_TDOAs().
        Must have 'path_info' in attrs.
    mode_string : str
        Propagation mode identifier (e.g., '2F2-1F2', '1F2-1E', '2F2-1E').
    tdoa_config : dict
        Configuration dictionary for this mode from build_tdoa_config().
        Must contain 'model_coeffs' key with (slope, intercept) tuple.
    output_dir : str
        Directory path where CSV file will be saved. Created if it doesn't exist.
    date_str : str, optional
        Date string in YYYYmmdd format for filename. If None, extracted from first chirp timestamp.

    Returns:
    csv_path : str
        Path to the saved CSV file.

    Example:
    >>> chirps = tdoa.find_chirps(wavlist, template, sweep_rate=10)
    >>> chirps = tdoa.find_TDOAs(chirps, mode_string='2F2-1F2')
    >>> tdoa_dct = tdoa.build_tdoa_config(chirps)
    >>> tdoa.save_tdoa_csv(chirps, '2F2-1F2', tdoa_dct['2F2-1F2'], 'output/tdoa_calculations')
    """
    # Get path info from chirps
    path_info = chirps.attrs.get('path_info')
    if path_info is None:
        raise ValueError("chirps must have 'path_info' in attrs")

    # Extract date from first chirp if not provided
    if date_str is None:
        first_utc = chirps['utc'].iloc[0]
        date_str = first_utc.strftime('%Y%m%d')

    # Get model coefficients
    slope, intercept = tdoa_config['model_coeffs']

    # Get mean TDOA values
    mean_tdoas = chirps[f'{mode_string}_mean']

    # Calculate layer heights from TDOAs using the model
    mean_heights = (slope * mean_tdoas) + intercept

    # Create output directory if it doesn't exist
    os.makedirs(output_dir, exist_ok=True)

    # Build filename: {YYYYmmdd}_{mode_string}_{path_string}.csv
    # Path string includes TX-RX-band information
    path_string = f"TX_{path_info.tx_call}_{path_info.tx_grid}-RX_{path_info.rx_call}_{path_info.rx_grid}-{path_info.band}m"
    filename = f"{date_str}_{mode_string}_{path_string}.csv"
    csv_path = os.path.join(output_dir, filename)

    # Build header block with metadata
    header_lines = []
    header_lines.append("# Automated TDOA Analysis Results")
    header_lines.append("#")
    header_lines.append(f"# Generated: {datetime.datetime.now().strftime('%Y-%m-%d %H:%M:%S')}")
    header_lines.append("#")
    header_lines.append("# Propagation Path Information:")
    header_lines.append(f"# TX: {path_info.tx_call}")
    header_lines.append(f"# TX Grid Square: {path_info.tx_grid}")
    header_lines.append(f"# RX: {path_info.rx_call}")
    header_lines.append(f"# RX Grid Square: {path_info.rx_grid}")
    header_lines.append(f"# Frequency: {path_info.band_str}")
    header_lines.append(f"# Ground Range: {path_info.get_range_km():.1f} km")
    header_lines.append("#")
    header_lines.append(f"# Propagation Mode: {mode_string}")
    header_lines.append("#")
    header_lines.append("# TDOA Model Coefficients:")
    header_lines.append(f"# Model: mean_tdoa_hgt_km = {slope:.1f}*mean_tdoa_ms + {intercept:.1f}")
    header_lines.append("#")
    header_lines.append("# Columns:")
    header_lines.append("# utc - UTC timestamp (YYYY-mm-dd HH:MM)")
    header_lines.append("# mean_tdoa_ms - Mean TDOA in milliseconds")
    header_lines.append("# mean_tdoa_hgt_km - Mean layer height in kilometers (calculated from TDOA)")
    header_lines.append("#")

    # Create DataFrame with results
    output_df = pd.DataFrame({
        'utc': chirps['utc'],
        'mean_tdoa_ms': mean_tdoas,
        'mean_tdoa_hgt_km': mean_heights
    })

    # Format UTC column for output
    output_df['utc'] = output_df['utc'].dt.strftime('%Y-%m-%d %H:%M')

    # Write to CSV with header
    with open(csv_path, 'w') as f:
        # Write header lines
        for line in header_lines:
            f.write(line + '\n')

        # Write data (pandas will write column names)
        output_df.to_csv(f, index=False, float_format='%.2f')

    print(f"Saved TDOA CSV: {csv_path}")

    return csv_path


# ============================================================================
# Signal Processing Functions
# ============================================================================

def filter(wav_df, low_pass_freq=250, high_pass_freq=10):
    """
    Filters a WAV signal using envelope detection and bandpass filtering.

    Arguments:
    wav_df : pd.DataFrame
        DataFrame containing the signal with 'x' column.
    low_pass_freq : float, optional
        Low pass filter cutoff frequency in Hz. Default is 250.
    high_pass_freq : float, optional
        High pass filter cutoff frequency in Hz. Default is 10.

    Returns:
    env_2 : pd.DataFrame
        Filtered envelope signal.
    fs : float
        Sampling frequency in Hz.
    """
    x0  = wav_df
    fs  = wav_df.attrs['fs']

    # Square the signal.
    x1  = x0**2

    # Compute envelope using a rolling max function.
    Tc  = 1.2e-3 # Rolling Max Window [seconds]
    env = x1.rolling(int(Tc*fs), center=True).max()
    env = env.dropna() # Drop NaNs

    wp    = low_pass_freq
    ws    = 1.1*wp

    gpass =  3 # The maximum loss in the passband (dB).
    gstop = 40 # The minimum attenuation in the stopband (dB).

    N, Wn = signal.buttord(wp, ws, gpass, gstop, fs=fs)
    sos   = signal.butter(N, Wn, 'low', fs=fs, output='sos')

    env_1      = env.copy()
    env_1['x'] = signal.sosfiltfilt(sos, env['x'])

    if high_pass_freq != 0:
        wp    = high_pass_freq
        ws    = 0.8*wp

        gpass =  3 # The maximum loss in the passband (dB).
        gstop = 40 # The minimum attenuation in the stopband (dB).

        N, Wn = signal.buttord(wp, ws, gpass, gstop, fs=fs)
        sos   = signal.butter(N, Wn, 'high', fs=fs, output='sos')

        env_2      = env_1.copy()
        env_2['x'] = signal.sosfiltfilt(sos, env_1['x'])
    else:
        env_2 = env_1

    return env_2, fs


def chirp_fft(df, tlim=None):
    """
    Calculates the FFT of a chirp signal.

    Arguments:
    df : pd.DataFrame
        DataFrame containing the chirp signal with time index and 'x' column.
    tlim : tuple, optional
        Time limits (start, end) to select a portion of the signal for FFT calculation. If None, uses the entire signal.

    Returns:
    X_psd : np.ndarray
        Power Spectral Density of the chirp signal.
    f : np.ndarray
        Frequency vector corresponding to the PSD.
    """
    env    = df['x']
    tvec = df.index
    Ts   = tvec[1] - tvec[0]
    Fs   = 1/Ts

    if tlim is None:
        tlim = (0, np.max(tvec))

    tf  = np.logical_and(tvec >= tlim[0], tvec <  tlim[1])
    xt  = env[tf].copy()

    han_win = np.hanning(len(xt))
    x_han   = han_win*xt
    nfft    = len(x_han)
    if nfft < 2**16:
        nfft = 2**16

    X_psd   = np.abs(np.fft.fftshift(np.fft.fft(x_han, n=nfft)*Ts*2))**2
    f       = np.fft.fftshift(np.fft.fftfreq(nfft, Ts))

    return X_psd, f


def find_max(freq, X_psd, minfreq, maxfreq):
    """
    Finds the maximum peak in a power spectral density within a frequency range.

    Arguments:
    freq : np.ndarray
        Frequency vector.
    X_psd : np.ndarray
        Power spectral density values.
    minfreq : float
        Minimum frequency for peak search.
    maxfreq : float
        Maximum frequency for peak search.

    Returns:
    maximum_x : float
        Frequency of maximum peak (or NaN if none found).
    maximum_y : float
        Amplitude of maximum peak (or NaN if none found).
    local_peaks_x : np.ndarray
        Frequencies of all peaks in range.
    local_peaks_y : np.ndarray
        Amplitudes of all peaks in range.
    """
    peaks = signal.find_peaks(X_psd)[0]
    local_peaks = pd.DataFrame({'f':freq[peaks], 'X_psd':X_psd[peaks]})
    local_peaks = local_peaks.sort_values(by='X_psd', ascending=False)

    local_peaks_x = local_peaks.f[local_peaks.f>minfreq][local_peaks.f<maxfreq].values
    local_peaks_y = local_peaks.X_psd[local_peaks.f>minfreq][local_peaks.f<maxfreq].values

    if len(local_peaks_x) > 0:
        maximum_x = local_peaks_x[0]
        maximum_y = local_peaks_y[0]
    else:
        maximum_x = np.nan
        maximum_y = np.nan

    return maximum_x, maximum_y, local_peaks_x, local_peaks_y


# ============================================================================
# TDOA Analysis Functions
# ============================================================================

def find_chirps(wavlist, template, sweep_rate, pfx=None, plot_correlation=False):
    """
    Finds chirps in a list of WAV files using cross-correlation with a template.
    This functions returns the top 10 chirp locations and their correlation coefficients for each WAV file.

    The function displays detailed processing information including:
    - Path information (TX/RX stations, range, band)
    - Template file and sweep rate
    - Progress bar showing file-by-file processing
    - Number of chirps found per file

    Arguments:
    wavlist : list
        List of paths to WAV files to be analyzed.
    template : str
        Path to the WAV file containing the chirp template.
    sweep_rate : float
        Sweep rate in Hz/ms. Needed for TDOA calculations.
    pfx : str or PathInfo, optional
        String or PathInfo object that identifies the transmit and receive station.
        If None, the pfx will be derived from the first name in the wavlist.
        A PathInfo object is stored as chirps.attrs['path_info'].
    plot_correlation : bool, optional
        If True, plots the correlation results for each WAV file. Default is False.

    Returns:
    chirps : pd.DataFrame
        DataFrame containing:
            'utc' : UTC timestamps
            'data': WAV file data as a pd.DataFrame
            'x'   : indices of chirp starts
            'y'   : cross-correlation coefficients at chirp starts
            'file': path to the source WAV file
        The DataFrame also has attrs['path_info'] containing a PathInfo object.
    """
    if pfx is None:
        pfx = os.path.basename(wavlist[0])[17:-4]

    # Convert pfx to PathInfo object if it's a string
    if isinstance(pfx, str):
        path_info = PathInfo(pfx)
    else:
        path_info = pfx

    # Print processing information
    print(f"\n{'='*70}")
    print(f"Finding Chirps via Cross-Correlation")
    print(f"{'='*70}")
    print(f"  Path Info:      {path_info}")
    print(f"  Template:       {os.path.basename(template)}")
    print(f"  Sweep Rate:     {sweep_rate} Hz/ms")
    print(f"  Files to Process: {len(wavlist)}")
    print(f"  Chirps per File: Top 10")
    print(f"{'='*70}\n")

    sample_env, sample_fs = load_wav(template)
    template_data = sample_env['x']

    # Here we are effectively taking time components
    utcs     = []
    chirps_x = []
    chirps_y = []
    wav_data = []

    # Use tqdm progress bar for file processing
    pbar = tqdm(wavlist, desc="Finding chirps", unit="file", disable=plot_correlation, dynamic_ncols=True)

    for file in pbar:
        bname   = os.path.basename(file)
        utc_str = bname[:13]
        utc     = datetime.datetime.strptime(utc_str, '%Y%m%d.%H%M')

        wav_env, wav_fs = load_wav(file)

        y = signal.correlate(wav_env['x'], template_data, mode='same')
        x = np.arange(0, len(y), 1)

        peaks       = signal.find_peaks(y, distance=300)[0]
        local_peaks = pd.DataFrame({'x':x[peaks], 'y':y[peaks]})
        local_peaks = local_peaks.sort_values(by='y', ascending=False)

        cx = local_peaks['x'][0:10].sort_index().values # Index locations
        cy = local_peaks['y'][0:10].sort_index().values # Correlation Amplitudes

        utcs.append(utc)
        chirps_x.append(cx)
        chirps_y.append(cy)
        wav_data.append(wav_env)

        # Update progress bar with number of chirps found and max correlation
        if not plot_correlation:
            max_corr = np.max(cy) if len(cy) > 0 else 0
            pbar.set_postfix({'Chirps': len(cx), 'Max Corr': f'{max_corr:.2e}'})

        if plot_correlation:
            fig = plt.figure(figsize=(16,9))
            ax = fig.add_subplot(1,1,1)
            ax.plot(wav_env.index, wav_env['x'])
            for i in range(0, 10):
                x0 = cx[i]/wav_fs
                y0 = cy/float(np.max(np.abs(cy)))
                ax.plot([x0, x0], [-y0, y0], color='tab:orange')
            ax.set_xlim(cx[0]/sample_fs-2.5, cx[-1]/sample_fs+2.5)
            ax.set_xlabel('Time [s]')
            ax.set_ylabel('Corr. Coef.')
            ax.set_title(bname)
            plt.show()
            plt.close(fig)

    if not plot_correlation:
        pbar.close()

    chirps = pd.DataFrame({'utc':utcs, 'data':wav_data, 'x':chirps_x, 'y':chirps_y, 'file':wavlist})
    chirps.attrs['path_info'] = path_info
    chirps.attrs['pfx'] = path_info.pfx  # Keep for backward compatibility
    chirps.attrs['sweep_rate'] = sweep_rate

    print(f"✓ Completed chirp detection: {len(chirps)} files processed, {len(chirps)*10} total chirps found\n")

    return chirps


def find_TDOAs(wav_data, search_limits=None, filter_limts=None,
               mode_string='2F2-1F2', set_name=None, plot_fft=False, plot_only_one=True, save_fft_dir=None, **overrides):
    """
    Finds TDOAs in a set of WAV files given chirp start locations.

    TDOA [ms] = (Beat Frequency [Hz]) / (Sweep Rate [Hz/ms])

    The function displays detailed processing information including:
    - Mode being processed and all parameters
    - Progress bar showing file-by-file processing
    - Mean TDOA value for each file being processed

    Arguments:
    wav_data : pd.DataFrame
        DataFrame containing WAV file data and chirp start locations.
    search_limits : tuple, optional
        (start_offset, end_offset, min_freq, max_freq)
        This is used to range gate the chirp signal in time to look for an echo from a particular ionospheric region.
        If None, uses the default for the specified mode_string.
        start_offset : float
            Start time offset (in seconds) from chirp start to begin FFT analysis.
        end_offset : float
            End time offset (in seconds) from chirp start to end FFT analysis.
        min_freq : float
            Minimum frequency for FFT peak search.
        max_freq : float
            Maximum frequency for FFT peak search.
    filter_limts : tuple, optional
        Frequency limits for bandpass filtering (low_pass_freq, high_pass_freq).
        If None, uses the default for the specified mode_string.
    mode_string : str, optional
        Propagation mode identifier (e.g., '2F2-1F2', '1F2-1E', '2F2-1E').
        Default is '2F2-1F2'. Used to look up default parameters from MODE_CONFIGS.
    set_name : str, optional
        Name of the column to store TDOA results. If None, uses mode_string as set_name.
    plot_fft : bool, optional
        If True, plots the FFT for each chirp analyzed. Default is False.
    plot_only_one : bool, optional
        If True, processes only the first chirp in first WAV file. Used for debugging.
    save_fft_dir : str, optional
        Directory path to save FFT plots. If provided, all FFT plots are saved with unique
        filenames instead of being displayed. Creates directory if it doesn't exist.
        Filename format: {mode_string}_{file_index:03d}_{chirp_index:02d}.png
        If None, FFT plots are only shown if plot_fft=True.
    **overrides : dict, optional
        Additional keyword arguments to override default mode configuration parameters.

    Returns:
    wav_data : pd.DataFrame
        Updated DataFrame with TDOA results added.
    """
    # Get default configuration for the mode
    if mode_string not in MODE_CONFIGS:
        raise ValueError(f"Unknown mode_string: {mode_string}. Available modes: {list(MODE_CONFIGS.keys())}")

    config = MODE_CONFIGS[mode_string].copy()

    # Apply any overrides
    config.update(overrides)

    # Use provided parameters or fall back to config defaults
    if search_limits is None:
        search_limits = config['search_limits']
    if filter_limts is None:
        filter_limts = config['filter_limts']
    if set_name is None:
        set_name = mode_string
    sweep_rate = wav_data.attrs['sweep_rate'] # Sweep rate in Hz/ms

    # Create save directory if specified
    if save_fft_dir is not None:
        os.makedirs(save_fft_dir, exist_ok=True)

    # Print processing information
    print(f"\n{'='*70}")
    print(f"Processing Mode: {mode_string}")
    print(f"{'='*70}")
    print(f"  Filter Limits:  {filter_limts[0]:.1f} - {filter_limts[1]:.1f} Hz")
    print(f"  Search Window:  {search_limits[0]:.2f} to {search_limits[1]:.2f} s offset")
    print(f"  Freq Range:     {search_limits[2]:.1f} - {search_limits[3]:.1f} Hz")
    print(f"  Set Name:       {set_name}")
    print(f"  Sweep Rate:     {sweep_rate} Hz/ms")
    print(f"  Files to Process: {len(wav_data)}")
    print(f"{'='*70}\n")

    all_beats  = []
    mean_beats = []

    # Use tqdm progress bar for file processing
    pbar = tqdm(wav_data.iterrows(), total=len(wav_data),
                desc=f"Finding TDOAs ({mode_string})",
                unit="file", dynamic_ncols=True)

    for idx, row in pbar:
        maxes = []
        wav_df = row['data']
        fpath  = row['file']  # Now correctly has unique filename per row
        bname  = os.path.basename(fpath)
        wav_env, wav_fs = filter(wav_df, filter_limts[1], filter_limts[0])

        for chirp_idx, peak in enumerate(row['x']):
            peak = peak/wav_fs
            tlim = (peak+search_limits[0], peak+search_limits[1])
            minfreq = search_limits[2]
            maxfreq = search_limits[3]

            X_psd, freq = chirp_fft(wav_env, tlim=tlim)

            maximum_x, maximum_y, local_peaks_x, local_peaks_y = find_max(freq, X_psd, minfreq, maxfreq)

            # Determine if we should plot/save FFT
            # Only plot/save if save_fft_dir is set (plot_fft is legacy behavior)
            # When plot_only_one=True, only plot the first chirp (chirp_idx == 0)
            should_plot = save_fft_dir is not None and (not plot_only_one or chirp_idx == 0)
            if should_plot:
                # Generate unique filename identifier (without .png extension)
                # Format: {mode_string}_{wav_filename}_chirp{###}
                # Each chirp is uniquely identified by its WAV file and chirp index
                fname_noext = os.path.splitext(bname)[0]
                plot_filename = f"{mode_string}_{fname_noext}_chirp{chirp_idx:03d}"
                savefig_path = os.path.join(save_fft_dir, f"{plot_filename}.png")

                plot_chirp_fft(title=bname, tlim=tlim,
                           minfreq=minfreq, maxfreq=maxfreq, wav_df=wav_df, sweep_rate=sweep_rate,
                           env=wav_env['x'], tvec=wav_env.index, X_psd=X_psd, f=freq,
                           savefig=savefig_path, plot_filename=plot_filename,
                           mode_string=mode_string, chirp_number=chirp_idx)

            maxes.append(maximum_x)

        beat_arr = np.array(maxes).flatten()/sweep_rate  # TDOA in ms
        all_beats.append(beat_arr)
        # Handle case where beat_arr contains only NaN values
        if len(beat_arr) > 0 and not np.all(np.isnan(beat_arr)):
            mean_beats.append(np.nanmean(beat_arr))
        else:
            mean_beats.append(np.nan)

        # Update progress bar with current file's mean TDOA
        if len(beat_arr) > 0 and not np.all(np.isnan(beat_arr)):
            pbar.set_postfix({'Mean TDOA': f'{np.nanmean(beat_arr):.2f} ms'})

    pbar.close()
    wav_data[set_name]           = all_beats
    wav_data[f'{set_name}_mean'] = mean_beats
    print(f"✓ Completed processing {mode_string}: {len(wav_data)} files processed\n")

    return wav_data


def build_tdoa_config(chirps, mode_strings=None, **mode_overrides):
    """
    Build a TDOA configuration dictionary for use with plot_TDOAs and plot_hmf2.

    This function creates a dictionary with mode-specific parameters and TDOA model coefficients.
    It also prints detailed path information and calculated model coefficients.

    Arguments:
    chirps : pd.DataFrame
        DataFrame containing chirp data. Must have 'path_info' in attrs.
    mode_strings : list of str, optional
        List of mode strings to include in the configuration (e.g., ['2F2-1F2', '1F2-1E']).
        If None, automatically detects which modes have been processed by checking for
        columns with '_mean' suffix in the chirps DataFrame.
    **mode_overrides : dict, optional
        Override parameters for specific modes. Use mode_string as key with a dict of parameters.
        Example: build_tdoa_config(chirps, **{'2F2-1F2': {'color': 'red', 'linewidth': 5}})

    Returns:
    tdoa_dct : dict
        Dictionary with mode configurations including model_coeffs and plotting parameters.
    """
    path_info = chirps.attrs.get('path_info')
    if path_info is None:
        raise ValueError("chirps must have 'path_info' in attrs")

    if mode_strings is None:
        # Auto-detect which modes have been processed by checking for columns ending with '_mean'
        mode_strings = []
        for col in chirps.columns:
            if col.endswith('_mean'):
                mode_name = col[:-5]  # Remove '_mean' suffix
                if mode_name in MODE_CONFIGS:
                    mode_strings.append(mode_name)

    # Print detailed path information
    print(f"\n{'='*60}")
    print(f"Path Information")
    print(f"{'='*60}")
    print(f"  {path_info}")
    print()

    # Get TX and RX coordinates from gridsquares
    tx_lat, tx_lon = path_info.get_tx_latlon()
    rx_lat, rx_lon = path_info.get_rx_latlon()
    print(f"TX Location: {tx_lat:.3f}°N, {tx_lon:.3f}°E ({path_info.tx_grid})")
    print(f"RX Location: {rx_lat:.3f}°N, {rx_lon:.3f}°E ({path_info.rx_grid})")
    print()

    # Calculate the great circle midpoint between TX and RX
    mid_lat, mid_lon = path_info.get_midpoint()
    print(f"Path Midpoint: {mid_lat:.3f}°N, {mid_lon:.3f}°E")
    print()

    # Calculate the azimuth from TX to RX
    azimuth = path_info.get_path_azimuth()
    print(f"Path Azimuth (TX→RX): {azimuth:.1f}°")
    print(f"{'='*60}\n")

    tdoa_dct = {}
    for mode_string in mode_strings:
        if mode_string not in MODE_CONFIGS:
            raise ValueError(f"Unknown mode_string: {mode_string}. Available modes: {list(MODE_CONFIGS.keys())}")

        # Start with the base config
        config = MODE_CONFIGS[mode_string].copy()

        # Calculate model coefficients
        config['mode_string'] = mode_string
        config['model_coeffs'] = path_info.calculate_TDOA_model(mode_string)

        # Apply any mode-specific overrides
        if mode_string in mode_overrides:
            config.update(mode_overrides[mode_string])

        tdoa_dct[mode_string] = config

    # Print calculated TDOA model coefficients for verification
    print("Calculated TDOA Model Coefficients:")
    print("=" * 60)
    for set_name, params in tdoa_dct.items():
        slope, intercept = params['model_coeffs']
        mode_str = params['mode_string']
        print(f"{set_name:12s} ({mode_str:8s}): slope={slope:6.1f}, intercept={intercept:6.1f}")
    print("=" * 60)
    print()

    return tdoa_dct


# ============================================================================
# Plotting Utility Functions
# ============================================================================

def _ensure_scalar(value):
    """
    Convert numpy array to scalar if needed.

    Parameters:
    -----------
    value : float or np.ndarray
        Value that may be a scalar or numpy array

    Returns:
    --------
    float
        Scalar value
    """
    if not np.isscalar(value):
        return float(value.item())
    return float(value)


def _format_datetime_axis(ax, time_format='%H:%M', interval_minutes=15):
    """
    Format datetime x-axis with consistent styling.

    Parameters:
    -----------
    ax : matplotlib.axes.Axes
        The axes object to format
    time_format : str, optional
        strftime format string for time labels (default: '%H:%M')
    interval_minutes : int, optional
        Interval in minutes for major tick marks (default: 15)
    """
    myFmt = mdates.DateFormatter(time_format)
    ax.xaxis.set_major_formatter(myFmt)
    ax.xaxis.set_major_locator(mdates.MinuteLocator(interval=interval_minutes))


def setup_plotting_style():
    """Sets default style and font parameters for plots."""
    mpl.rcParams['font.size']          = 16
    mpl.rcParams['font.weight']        = 'bold'
    mpl.rcParams['axes.labelweight']   = 'bold'
    mpl.rcParams['axes.titleweight']   = 'bold'
    mpl.rcParams['figure.labelweight'] = 'bold'
    mpl.rcParams['figure.titleweight'] = 'bold'
    mpl.rcParams['axes.grid']          = True
    mpl.rcParams['grid.linestyle']     = ':'
    mpl.rcParams['figure.figsize']     = (12,9)
    mpl.rcParams['axes.xmargin']       = 0
    mpl.rcParams['axes.ymargin']       = 0.1


def title_from_pfx(ax, pfx, date=None, center_title=None):
    """
    Creates a formatted title for plots based on the prefix string.

    Arguments:
    ax : matplotlib.axes.Axes
        The axes object to add the title to.
    pfx : str or PathInfo
        Prefix string containing TX/RX station information, or a PathInfo object.
    date : datetime.datetime, optional
        Date to display in the center title.
    center_title : str, optional
        Custom text to display above the date in the center title (e.g., 'TDOA' or 'Ionospheric Layer Height').
        If provided, creates a two-line center title with this text on top and the date below.
    """
    # Accept either a string or a PathInfo object
    if isinstance(pfx, str):
        path_info = PathInfo(pfx)
    else:
        path_info = pfx

    title = f'TX: {path_info.tx_call} ({path_info.tx_grid})\nRX: {path_info.rx_call} ({path_info.rx_grid})'
    ax.set_title(title, loc='left')

    range_km = path_info.get_range_km()
    title = f'Ground Range: {range_km:.0f} km\nBand: {path_info.band_str}'
    ax.set_title(title, loc='right')

    if date is not None:
        date_str = date.strftime('%Y %b %d')
        if center_title is not None:
            # Two-line center title with custom text on top and date below
            center_text = f'{center_title}\n{date_str}'
        else:
            # Just the date
            center_text = date_str
        ax.set_title(center_text, loc='center', fontsize=26)


def plot_chirp_fft(title, tlim, minfreq, maxfreq, wav_df, sweep_rate, env, tvec, X_psd, f, savefig=None, plot_filename=None, mode_string=None, chirp_number=None):
    """
    Plots the FFT of a chirp signal along with its envelope and raw waveform.

    Arguments:
    title : str
        Title for the plot (WAV filename).
    tlim : tuple
        Time limits (start, end) for the chirp signal.
    minfreq : float
        Minimum frequency for FFT plot.
    maxfreq : float
        Maximum frequency for FFT plot.
    wav_df : pd.DataFrame or None
        DataFrame containing the raw WAV signal. If None, raw signal is not plotted.
    sweep_rate : float
        Sweep rate in Hz/ms for TDOA calculations.
    env : pd.Series
        Envelope of the chirp signal.
    tvec : np.ndarray
        Time vector corresponding to the envelope.
    X_psd : np.ndarray
        Power Spectral Density of the chirp signal.
    f : np.ndarray
        Frequency vector corresponding to the PSD.
    savefig : str, optional
        Path to save the figure. If None, figure is displayed but not saved.
    plot_filename : str, optional
        Filename identifier (without extension) to display in the title.
    mode_string : str, optional
        Mode string (e.g., '2F2-1F2') to display in the title.
    chirp_number : int, optional
        Chirp number to display in the title.
    """
    if maxfreq == 0:
        maxfreq = 10000
    flim = (0, maxfreq)

    nrows = 3
    nax   = 0
    fig = plt.figure(figsize=(15,15))

    nax += 1
    ax1  = fig.add_subplot(nrows, 1, nax)
    if wav_df is not None:
        xx = wav_df.index
        yy = wav_df['x']
        ax1.plot(xx, yy)
        ax1.set_xlabel('Time [s]')
        ax1.set_ylabel('Amplitude')
        ax1.axvspan(*tlim, color='red', label='Selected Chirp')
        ax1.legend(loc='upper right')

    nax += 1
    ax2  = fig.add_subplot(nrows, 1, nax)

    dc_offset = 0
    if wav_df is not None:
        tf  = np.logical_and(wav_df.index >= tlim[0], wav_df.index <  tlim[1])
        wav_xx = wav_df.index[tf]
        wav_yy = wav_df['x'][tf]
        ax2.plot(wav_xx, wav_yy, label='Raw WAV')
        dc_offset = np.mean(np.abs(wav_yy))

    env_xx = tvec
    env_yy = env + dc_offset
    ax2.plot(env_xx, env_yy, lw=3, color='red', label='Filtered Envelope')
    ax2.set_xlim(tlim)
    ax2.set_xlabel('Time [s]')
    ax2.set_ylabel('Amplitude')
    ax2.legend(loc='upper right')
    ax2.set_ylim(0, None)

    nax += 1
    ax3 = fig.add_subplot(nrows, 1, nax)
    ax3.plot(f, X_psd)

    maximum_x, maximum_y, local_peaks_x, local_peaks_y = find_max(f, X_psd, minfreq, maxfreq)

    ax3.plot([maximum_x, maximum_x], [0, maximum_y], linewidth=2)
    ax3.plot([0, maximum_x], [maximum_y, maximum_y], linewidth=2, color='red')
    ax3.plot([minfreq, minfreq], [0, maximum_y], linewidth=3, color='black')

    lbl = []
    lbl.append(f'$f_b = {maximum_x:0.2f}$ Hz')
    lbl.append(f'TDOA = {maximum_x/sweep_rate:0.2f} ms')
    ax3.scatter(maximum_x, maximum_y, linewidths=5, label='\n'.join(lbl), color='green')
    ax3.scatter(local_peaks_x, local_peaks_y, color='black')
    ax3.set_xlim(flim)
    ax3.set_ylim(0, None)
    ax3.set_xlabel('Frequency [Hz]')
    ax3.set_ylabel('$|X(f)|^2$')

    # Make x-labels with both beat frequency and TDOA values.
    xticks = ax3.get_xticks()
    ax3.set_xticks(xticks)
    xtls = []
    for xtk in xticks:
        tdoa = xtk/sweep_rate
        xtl = f'{xtk:0.1f}\n{tdoa:0.1f}'
        xtls.append(xtl)
    xtls[-1] = 'f [Hz]\nTDOA [ms]'
    ax3.set_xticklabels(xtls)
    ax3.set_xlabel('')

    ax3.legend()

    _title = []
    _title.append(title)  # WAV filename
    if mode_string is not None and chirp_number is not None:
        _title.append(f'Mode: {mode_string}  Chirp #{chirp_number:03d}')
    _title.append(f'Chirp Sweep Rate: {sweep_rate} Hz/ms')
    plt.suptitle('\n'.join(_title))

    plt.tight_layout()

    # Add subplot labels (a), (b), (c) after layout adjustment
    ypos = 1.12
    ax1.text(-0.08, ypos, '(a)', transform=ax1.transAxes,
            fontsize=30, fontweight='bold', va='top', ha='left')
    ax2.text(-0.08, ypos, '(b)', transform=ax2.transAxes,
            fontsize=30, fontweight='bold', va='top', ha='left')
    ax3.text(-0.08, ypos, '(c)', transform=ax3.transAxes,
            fontsize=30, fontweight='bold', va='top', ha='left')

    if savefig is not None:
        plt.savefig(savefig, dpi=150, bbox_inches='tight')
    else:
        plt.show()

    plt.close(fig)


# ============================================================================
# Main Plotting Functions - TDOA
# ============================================================================

def _plot_tdoa_on_axis(ax, chirps, tdoa_dct, ylim=(0, 3), show_date=True, legend_offset=-0.06):
    """
    Internal helper function to plot TDOA measurements on a given axis.

    This function is used by both plot_TDOAs() and plot_tdoa_hmf2_subplot() to avoid code duplication.

    Arguments:
    ax : matplotlib.axes.Axes
        The axes object to plot on.
    chirps : pd.DataFrame
        DataFrame containing chirp data and TDOA measurements.
    tdoa_dct : dict
        Dictionary containing TDOA set configurations with plotting parameters.
    ylim : tuple, optional
        Y-axis limits for the plot. Default is (0, 3).
    show_date : bool, optional
        If True, show date in the title. Default is True.
    legend_offset : float, optional
        Vertical offset for legend position as fraction of axis height. Default is -0.06.
        More negative values move the legend down.
    """
    times = chirps['utc']

    lgnds = []
    for set_name, params in tdoa_dct.items():
        TDOAs = chirps[f'{set_name}_mean']
        line, = ax.plot(times, TDOAs,
                label=f'Automated Analysis: {set_name}',
                marker=params.get('marker', 'o'),
                linestyle=params.get('linestyle'),
                linewidth=params.get('linewidth'),
                color=params.get('color'))

        lgnd = {}
        lgnd['line'] = line
        lgnd['xpos'] = 0.01
        lgnd['ypos'] = np.nanmin(TDOAs)
        lgnds.append(lgnd)

    ax.set_ylim(ylim)

    # Add a Legend Under Each Trace
    ylim_actual = ax.get_ylim()
    for lgnd in lgnds:
        line = lgnd['line']
        xpos = lgnd['xpos']
        ypos = lgnd['ypos']
        ypos = (ypos - ylim_actual[0]) / (ylim_actual[1] - ylim_actual[0])  # Convert to fraction of axis height
        ypos = ypos + legend_offset  # Nudge legend down
        legend = ax.legend(handles=[line], loc=(xpos, ypos), fontsize='small')
        ax.add_artist(legend)

    _format_datetime_axis(ax)

    ax.set_ylabel('TDOA [ms]')
    ax.set_xlabel('Time UTC')

    # Use path_info from attrs
    path_info = chirps.attrs['path_info']

    # Show date if requested
    date = times.iloc[0] if show_date else None
    title_from_pfx(ax, path_info, date, center_title='TDOA')


def plot_TDOAs(chirps, tdoa_dct, ylim=(0,3), savefig=None):
    """
    Plots TDOA measurements over time for multiple propagation modes.

    Arguments:
    chirps : pd.DataFrame
        DataFrame containing chirp data and TDOA measurements.
    tdoa_dct : dict
        Dictionary containing TDOA set configurations with plotting parameters.
    ylim : tuple, optional
        Y-axis limits for the plot. Default is (0, 3).
    savefig : str, optional
        If provided, saves the figure to this file path. High-resolution JPEG recommended.
        If None (default), figure is not saved to disk.
    """
    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(1, 1, 1)

    # Use the shared helper function to do all the plotting
    _plot_tdoa_on_axis(ax, chirps, tdoa_dct, ylim=ylim, show_date=True)

    fig.autofmt_xdate()
    plt.tight_layout()

    # Save figure if filename provided
    if savefig is not None:
        plt.savefig(savefig, dpi=300, bbox_inches='tight', format='jpeg', pil_kwargs={'quality': 95})
        print(f"Figure saved to: {savefig}")

    plt.show()
    plt.close(fig)


# ============================================================================
# Main Plotting Functions - Layer Heights (hmf2)
# ============================================================================

def _plot_hmf2_on_axis(ax, chirps, tdoa_dct, ylim=(75,450),
                       solar_lat=None, solar_lon=None, solar_start=None, solar_end=None,
                       overlay_solar_elevation=False, overlay_eclipse=False,
                       ionosonde_dct=None, tdoa_csv_dct=None, show_date=True,
                       legend_loc='best', legend_fontsize=None):
    """
    Internal helper function to plot HF TDOA layer heights on a given axis.

    This function is used by both plot_hmf2() and plot_hmf2_subplot() to avoid code duplication.

    Arguments:
    ax : matplotlib.axes.Axes
        The axes object to plot on.
    chirps : pd.DataFrame
        DataFrame containing chirp data and TDOA measurements.
    tdoa_dct : dict
        Dictionary containing TDOA set configurations with model coefficients and plotting parameters.
    ylim : tuple, optional
        Y-axis limits for the plot. Default is (75, 450).
    solar_lat : float, optional
        Latitude for solar calculations (degrees, +N/-S).
    solar_lon : float, optional
        Longitude for solar calculations (degrees, +E/-W).
    solar_start : datetime.datetime, optional
        Start time for solar calculations. If None, derived from plot x-axis limits.
    solar_end : datetime.datetime, optional
        End time for solar calculations. If None, derived from plot x-axis limits.
    overlay_solar_elevation : bool, optional
        If True, overlay solar elevation angle. Default is False.
    overlay_eclipse : bool, optional
        If True, overlay eclipse obscuration. Default is False.
    ionosonde_dct : dict, optional
        Dictionary containing parameters to pass to overlay_ionosonde().
    tdoa_csv_dct : dict, optional
        Dictionary containing parameters to pass to overlay_tdoa_csv().
        Supported keys: csv_path_period, csv_path_autocorr, label_period, label_autocorr,
        color_period, linestyle_period, color_autocorr, linestyle_autocorr.
    show_date : bool, optional
        If True, show date in the title. Default is True.
    legend_loc : str, optional
        Legend location. Default is 'best'. Can be 'upper left', 'upper right', 'lower left',
        'lower right', 'center', 'best', etc.
    legend_fontsize : int or str, optional
        Legend font size. If None, uses matplotlib default. Can be int or 'small', 'medium', 'large', etc.
    """
    times = chirps['utc']

    # Overlay TDOA CSV data first (Manual Period and Autocorrelation Analysis)
    if (tdoa_csv_dct is not None) and (tdoa_csv_dct is not False):
        if tdoa_csv_dct is True:
            tdoa_csv_dct = {}
        overlay_tdoa_csv(ax, **tdoa_csv_dct)

    # Plot TDOA-derived layer heights (Automated Analysis)
    # Sort modes by mean TDOA value (longest to shortest) for consistent legend ordering
    sorted_modes = sorted(tdoa_dct.items(),
                         key=lambda x: chirps[f'{x[0]}_mean'].mean(),
                         reverse=True)

    for set_name, params in sorted_modes:
        TDOAs = chirps[f'{set_name}_mean']
        coefs = params['model_coeffs']
        layer_heights = (coefs[0] * TDOAs) + coefs[1]

        ax.plot(times, layer_heights,
                label=f'Automated Analysis: {set_name}',
                marker=params.get('marker', 'o'),
                linestyle=params.get('linestyle'),
                linewidth=params.get('linewidth'),
                color=params.get('color'))

    # Overlay ionosonde data last (Austin hmF2)
    if (ionosonde_dct is not None) and (ionosonde_dct is not False):
        if ionosonde_dct is True:
            ionosonde_dct = {}
        overlay_ionosonde(ax, **ionosonde_dct)

    ax.set_ylim(ylim)

    _format_datetime_axis(ax)
    ax.set_ylabel('Layer Height [km]')
    ax.set_xlabel('Time UTC')

    # Add solar elevation and/or eclipse obscuration overlays if requested
    # Collect line objects to add to the legend
    solar_lines = []
    solar_labels = []

    if overlay_solar_elevation or overlay_eclipse:
        if solar_lat is None or solar_lon is None:
            print('WARNING: solar_lat and solar_lon must be provided for solar overlays.')
        else:
            try:
                from . import solarContext

                # Set default times based on x-axis limits if not provided
                if solar_start is None or solar_end is None:
                    xlim = ax.get_xlim()
                    if solar_start is None:
                        solar_start = mdates.num2date(xlim[0]).replace(tzinfo=None)
                    if solar_end is None:
                        solar_end = mdates.num2date(xlim[1]).replace(tzinfo=None)

                # Create solarTimeseries object
                solar_ts = solarContext.solarTimeseries(
                    sTime=solar_start,
                    eTime=solar_end,
                    lat=solar_lat,
                    lon=solar_lon
                )

                if overlay_solar_elevation:
                    line = solar_ts.overlaySolarElevation(ax)
                    if line is not None:
                        solar_lines.append(line)
                        solar_labels.append(line.get_label())

                if overlay_eclipse:
                    line = solar_ts.overlayEclipse(ax)
                    if line is not None:
                        solar_lines.append(line)
                        solar_labels.append(line.get_label())

            except (ImportError, AttributeError) as e:
                print(f'WARNING: Cannot overlay solar data. Error: {e}')

    # Create legend with all elements including solar overlays
    if solar_lines:
        # Get existing legend handles and labels from the main axis
        handles, labels = ax.get_legend_handles_labels()
        # Add solar overlay lines to the legend
        handles.extend(solar_lines)
        labels.extend(solar_labels)
        ax.legend(handles=handles, labels=labels, loc=legend_loc, fontsize=legend_fontsize)
    else:
        ax.legend(loc=legend_loc, fontsize=legend_fontsize)

    # Use path_info from attrs
    path_info = chirps.attrs['path_info']

    # Show date if requested
    date = times.iloc[0] if show_date else None
    title_from_pfx(ax, path_info, date, center_title='Ionospheric Layer Height')


def plot_hmf2(chirps, tdoa_dct, ylim=(75,450), xlim=None,
              solar_lat=None, solar_lon=None, solar_start=None, solar_end=None,
              overlay_solar_elevation=False, overlay_eclipse=False,
              ionosonde_dct=None, tdoa_csv_dct=None,
              legend_loc='best', legend_fontsize=None,
              savefig=None):
    """
    Plots layer heights derived from TDOAs and compares with ionosonde measurements.

    Arguments:
    chirps : pd.DataFrame
        DataFrame containing chirp data and TDOA measurements.
    tdoa_dct : dict
        Dictionary containing TDOA set configurations with model coefficients and plotting parameters.
    ylim : tuple, optional
        Y-axis limits for the plot. Default is (75, 450).
    xlim : tuple of datetime, optional
        X-axis limits for the plot as (start_datetime, end_datetime). Default is None (auto).
    solar_lat : float, optional
        Latitude for solar calculations (degrees, +N/-S). Required if overlay_solar_elevation or overlay_eclipse is True.
    solar_lon : float, optional
        Longitude for solar calculations (degrees, +E/-W). Required if overlay_solar_elevation or overlay_eclipse is True.
    solar_start : datetime.datetime, optional
        Start time for solar calculations. If None, derived from plot x-axis limits.
    solar_end : datetime.datetime, optional
        End time for solar calculations. If None, derived from plot x-axis limits.
    overlay_solar_elevation : bool, optional
        If True, overlay solar elevation angle. Default is False.
    overlay_eclipse : bool, optional
        If True, overlay eclipse obscuration. Default is False.
    ionosonde_dct : dict, optional
        Dictionary containing parameters to pass to overlay_ionosonde().
        If True, uses default parameters. If None, no ionosonde data is overlaid.
        Available parameters:
            - csv_path: Path to ionosonde CSV file
            - overlay_hmF2: If True, overlay hmF2 data (default: True)
            - overlay_hmE: If True, overlay hmE data (default: True)
            - label: Label prefix for ionosonde data (default: 'Austin Ionosonde')
            - hmF2_color: Color for hmF2 line (default: 'purple')
            - hmE_color: Color for hmE line (default: 'brown')
        Example: {'csv_path': 'path/to/file.csv', 'label': 'Boulder Ionosonde',
                  'overlay_hmF2': True, 'overlay_hmE': False, 'hmF2_color': 'red'}
        Default is None.
    tdoa_csv_dct : dict, optional
        Dictionary containing parameters to pass to overlay_tdoa_csv().
        If True, uses default parameters. If None, no TDOA CSV data is overlaid.
        Available parameters:
            - csv_path_period: Path to manual period analysis CSV file (optional)
            - csv_path_autocorr: Path to manual autocorrelation analysis CSV file (optional)
            - label_period: Label for period analysis (default: 'Manual Period Analysis')
            - label_autocorr: Label for autocorrelation (default: 'Manual Autocorrelation Analysis')
            - color_period: Color for period line (default: 'tab:blue')
            - linestyle_period: Linestyle for period (default: 'dotted')
            - color_autocorr: Color for autocorr line (default: 'tab:orange')
            - linestyle_autocorr: Linestyle for autocorr (default: 'dashdot')
        Example: {'csv_path_period': 'data/CSVs/period.csv',
                  'csv_path_autocorr': 'data/CSVs/autocorr.csv'}
        Default is None.
    legend_loc : str, optional
        Legend location. Default is 'best'. Can be 'upper left', 'upper right', 'lower left',
        'lower right', 'center', 'best', etc.
    legend_fontsize : int or str, optional
        Legend font size. If None, uses matplotlib default. Can be int or 'small', 'medium', 'large', etc.
    savefig : str, optional
        If provided, saves the figure to this file path. High-resolution JPEG recommended.
        If None (default), figure is not saved to disk.
    """
    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_subplot(1, 1, 1)

    # Use the shared helper function to do all the plotting
    _plot_hmf2_on_axis(
        ax, chirps, tdoa_dct, ylim=ylim,
        solar_lat=solar_lat, solar_lon=solar_lon,
        solar_start=solar_start, solar_end=solar_end,
        overlay_solar_elevation=overlay_solar_elevation,
        overlay_eclipse=overlay_eclipse,
        ionosonde_dct=ionosonde_dct,
        tdoa_csv_dct=tdoa_csv_dct,
        show_date=True,
        legend_loc=legend_loc,
        legend_fontsize=legend_fontsize
    )

    # Set x-axis limits if provided
    if xlim is not None:
        ax.set_xlim(xlim)

    fig.autofmt_xdate()
    plt.tight_layout()

    # Save figure if filename provided
    if savefig is not None:
        plt.savefig(savefig, dpi=300, bbox_inches='tight', format='jpeg', pil_kwargs={'quality': 95})
        print(f"Figure saved to: {savefig}")

    plt.show()
    plt.close(fig)


# ============================================================================
# Advanced Multi-Panel Plotting Functions
# ============================================================================

def plot_hmf2_subplot(chirps_list, tdoa_dct_list, subplot_labels=None, ylim=(200, 350), xlim=None,
                      solar_lat=None, solar_lon=None, solar_start=None, solar_end=None,
                      overlay_solar_elevation=False, overlay_eclipse=False,
                      ionosonde_dct=None, tdoa_csv_dct_list=None, figsize=(15, 16), savefig=None,
                      legend_loc='best', legend_fontsize=None):
    """
    Creates a multi-panel subplot figure with layer heights from multiple datasets.

    This is useful for creating figures like Figure 12 with AB5YO 40m and 60m.

    Arguments:
    chirps_list : list of pd.DataFrame
        List of chirps DataFrames, one for each subplot.
    tdoa_dct_list : list of dict
        List of TDOA configuration dictionaries, one for each subplot.
    subplot_labels : list of str, optional
        Labels for each subplot (e.g., ['(a)', '(b)']). If None, uses (a), (b), (c), etc.
    ylim : tuple, optional
        Y-axis limits for all plots. Default is (200, 350).
    xlim : tuple of datetime.datetime, optional
        X-axis limits (start_time, end_time) applied to all subplots.
        If None, uses the full time range from the data.
    solar_lat : float, optional
        Latitude for solar calculations (degrees, +N/-S).
    solar_lon : float, optional
        Longitude for solar calculations (degrees, +E/-W).
    solar_start : datetime.datetime, optional
        Start time for solar calculations.
    solar_end : datetime.datetime, optional
        End time for solar calculations.
    overlay_solar_elevation : bool, optional
        If True, overlay solar elevation angle. Default is False.
    overlay_eclipse : bool, optional
        If True, overlay eclipse obscuration. Default is False.
    ionosonde_dct : dict, optional
        Dictionary containing parameters to pass to overlay_ionosonde() for all subplots.
    tdoa_csv_dct_list : list of dict, optional
        List of TDOA CSV dictionaries, one for each subplot. If None, no CSV data is overlaid.
    figsize : tuple, optional
        Figure size (width, height). Default is (15, 16).
    savefig : str, optional
        If provided, saves the figure to this file path. High-resolution JPEG recommended.
        If None (default), figure is not saved to disk.
    legend_loc : str, optional
        Legend location for all subplots. Default is 'best'. Can be 'upper left', 'upper right',
        'lower left', 'lower right', 'center', 'best', etc.
    legend_fontsize : int or str, optional
        Legend font size for all subplots. If None, uses matplotlib default.
        Can be int or 'small', 'medium', 'large', etc.
    """
    import datetime

    n_plots = len(chirps_list)

    if subplot_labels is None:
        subplot_labels = [f'({chr(97 + i)})' for i in range(n_plots)]  # (a), (b), (c), ...

    if tdoa_csv_dct_list is None:
        tdoa_csv_dct_list = [None] * n_plots

    fig, axes = plt.subplots(n_plots, 1, figsize=figsize)

    # Ensure axes is iterable even if n_plots == 1
    if n_plots == 1:
        axes = [axes]

    for idx, (chirps, tdoa_dct, ax, label, tdoa_csv_dct) in enumerate(
            zip(chirps_list, tdoa_dct_list, axes, subplot_labels, tdoa_csv_dct_list)):

        # Use the shared helper function to do all the plotting
        _plot_hmf2_on_axis(
            ax, chirps, tdoa_dct, ylim=ylim,
            solar_lat=solar_lat, solar_lon=solar_lon,
            solar_start=solar_start, solar_end=solar_end,
            overlay_solar_elevation=overlay_solar_elevation,
            overlay_eclipse=overlay_eclipse,
            ionosonde_dct=ionosonde_dct,
            tdoa_csv_dct=tdoa_csv_dct,
            show_date=True,
            legend_loc=legend_loc,
            legend_fontsize=legend_fontsize
        )

        # Set x-axis limits if specified
        if xlim is not None:
            ax.set_xlim(xlim)

        # Add subplot label (a), (b), etc. using ax.text outside the axes
        # Position it in axes coordinates - outside upper left corner
        ax.text(-0.08, 1.075, label, transform=ax.transAxes,
                fontsize=30, fontweight='bold', va='top', ha='left')

        # Rotate x-axis labels
        for tick_label in ax.get_xticklabels():
            tick_label.set_rotation(45)
            tick_label.set_horizontalalignment('right')

    plt.tight_layout()
    plt.subplots_adjust(hspace=0.3)

    # Save figure if filename provided
    if savefig is not None:
        plt.savefig(savefig, dpi=300, bbox_inches='tight', format='jpeg', pil_kwargs={'quality': 95})
        print(f"Figure saved to: {savefig}")

    plt.show()
    plt.close(fig)


def plot_tdoa_hmf2_subplot(chirps, tdoa_dct, subplot_labels=None,
                           ylim_tdoa=(0, 5), ylim_hmf2=(200, 350), xlim=None,
                           legend_offset=-0.06,
                           solar_lat=None, solar_lon=None, solar_start=None, solar_end=None,
                           overlay_solar_elevation=False, overlay_eclipse=False,
                           ionosonde_dct=None, tdoa_csv_dct=None,
                           image_panel=None,
                           figsize=(15, 16), savefig=None):
    """
    Creates a two-panel (or three-panel) subplot figure with TDOA measurements (top) and layer heights (bottom).

    This is useful for creating figures that show both raw TDOA measurements and the
    derived layer heights in a single combined figure. Optionally includes an image panel.

    Arguments:
    chirps : pd.DataFrame
        DataFrame containing chirp data and TDOA measurements.
    tdoa_dct : dict
        Dictionary containing TDOA set configurations with model coefficients and plotting parameters.
    subplot_labels : list of str, optional
        Labels for each subplot (e.g., ['(a)', '(b)', '(c)']). If None, auto-generates based on number of panels.
    ylim_tdoa : tuple, optional
        Y-axis limits for TDOA plot. Default is (0, 5).
    ylim_hmf2 : tuple, optional
        Y-axis limits for layer height plot. Default is (200, 350).
    xlim : tuple of datetime.datetime, optional
        X-axis limits (start_time, end_time) applied to both TDOA and layer height plots.
        If None, uses the full time range from the data.
    legend_offset : float, optional
        Vertical offset for TDOA legend position as fraction of axis height. Default is -0.06.
        More negative values move the legend down. Use this to avoid overlap with data.
    solar_lat : float, optional
        Latitude for solar calculations (degrees, +N/-S).
    solar_lon : float, optional
        Longitude for solar calculations (degrees, +E/-W).
    solar_start : datetime.datetime, optional
        Start time for solar calculations.
    solar_end : datetime.datetime, optional
        End time for solar calculations.
    overlay_solar_elevation : bool, optional
        If True, overlay solar elevation angle on hmf2 plot. Default is False.
    overlay_eclipse : bool, optional
        If True, overlay eclipse obscuration on hmf2 plot. Default is False.
    ionosonde_dct : dict, optional
        Dictionary containing parameters to pass to overlay_ionosonde() for hmf2 plot.
    tdoa_csv_dct : dict, optional
        Dictionary containing parameters to pass to overlay_tdoa_csv() for hmf2 plot.
        See plot_hmf2() documentation for available parameters.
    image_panel : str, optional
        Path to an image file to include as an additional panel (e.g., 'etalon.png').
        If provided, creates a third panel below the layer heights plot.
    figsize : tuple, optional
        Figure size (width, height). Default is (15, 16).
    savefig : str, optional
        If provided, saves the figure to this file path. High-resolution JPEG recommended.
        If None (default), figure is not saved to disk.
    """
    import datetime
    from matplotlib import image as mpimg

    # Determine number of panels
    n_panels = 3 if image_panel is not None else 2

    if subplot_labels is None:
        subplot_labels = [f'({chr(97 + i)})' for i in range(n_panels)]  # (a), (b), (c), ...

    # Adjust figure size if we have 3 panels
    if n_panels == 3 and figsize == (15, 16):
        figsize = (15, 20)  # Make taller for 3 panels

    fig, axes = plt.subplots(n_panels, 1, figsize=figsize)

    # Ensure axes is always iterable
    if n_panels == 1:
        axes = [axes]

    ax_tdoa = axes[0]
    ax_hmf2 = axes[1]

    # Plot (a): TDOA measurements using the helper function
    _plot_tdoa_on_axis(ax_tdoa, chirps, tdoa_dct, ylim=ylim_tdoa, show_date=True, legend_offset=legend_offset)

    # Plot (b): Layer heights using the helper function
    _plot_hmf2_on_axis(
        ax_hmf2, chirps, tdoa_dct, ylim=ylim_hmf2,
        solar_lat=solar_lat, solar_lon=solar_lon,
        solar_start=solar_start, solar_end=solar_end,
        overlay_solar_elevation=overlay_solar_elevation,
        overlay_eclipse=overlay_eclipse,
        ionosonde_dct=ionosonde_dct,
        tdoa_csv_dct=tdoa_csv_dct,
        show_date=True
    )

    # Plot (c): Image panel if provided
    if image_panel is not None:
        ax_img = axes[2]
        img = mpimg.imread(image_panel)
        ax_img.imshow(img)
        ax_img.axis('off')  # Turn off axis for image panel

    # Set x-axis limits for both time series plots if specified
    if xlim is not None:
        ax_tdoa.set_xlim(xlim)
        ax_hmf2.set_xlim(xlim)

    # Rotate x-axis labels for time series subplots (not image panel)
    for ax in axes[:2]:  # Only first two panels have time axes
        for tick_label in ax.get_xticklabels():
            tick_label.set_rotation(45)
            tick_label.set_horizontalalignment('right')

    # Apply tight_layout first to calculate proper spacing
    plt.tight_layout()

    # Adjust spacing between subplots to prevent title overlap
    plt.subplots_adjust(hspace=0.4)

    # Add subplot labels AFTER layout adjustment
    ypos = 1.100
    ax_tdoa.text(-0.08, ypos, subplot_labels[0], transform=ax_tdoa.transAxes,
                fontsize=30, fontweight='bold', va='top', ha='left')

    ax_hmf2.text(-0.08, ypos, subplot_labels[1], transform=ax_hmf2.transAxes,
                fontsize=30, fontweight='bold', va='top', ha='left')

    if image_panel is not None:
        ax_img.text(-0.08, ypos, subplot_labels[2], transform=ax_img.transAxes,
                    fontsize=30, fontweight='bold', va='top', ha='left')

    # Save figure if filename provided
    if savefig is not None:
        plt.savefig(savefig, dpi=300, bbox_inches='tight', format='jpeg', pil_kwargs={'quality': 95})
        print(f"Figure saved to: {savefig}")

    plt.show()
    plt.close(fig)


# ============================================================================
# Data Overlay Functions for Validation
# ============================================================================

def overlay_ionosonde(ax, csv_path=None,
                      overlay_hmF2=True, overlay_hmE=True,
                      label='Austin Ionosonde',
                      hmF2_color='purple', hmE_color='brown'):
    """
    Overlays ionosonde data (hmF2 and hmE) on an existing axes.

    Arguments:
    ax : matplotlib.axes.Axes
        The axes object to plot on.
    csv_path : str, optional
        Path to the ionosonde CSV file. If None, uses default path relative to package.
    overlay_hmF2 : bool, optional
        If True, overlay hmF2 data. Default is True.
    overlay_hmE : bool, optional
        If True, overlay hmE data. Default is True.
    label : str, optional
        Label prefix for the ionosonde data. Default is 'Austin Ionosonde'.
    hmF2_color : str, optional
        Color for hmF2 line. Default is 'purple'.
    hmE_color : str, optional
        Color for hmE line. Default is 'brown'.
    """
    # If no path provided, use default relative to package root
    if csv_path is None:
        # Get the directory containing this file (hf_tdoa package directory)
        package_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        csv_path = os.path.join(package_dir, 'data', 'CSVs', '2024-04-08_AU930_AustinTX_Ionosonde_ManualScaled.csv')

    ionosonde_df = load_ionosonde_data(csv_path)

    if overlay_hmF2:
        ax.plot(ionosonde_df.index, ionosonde_df['hmF2'], color=hmF2_color,
                label=f"{label} hmF2", linewidth=2)

    if overlay_hmE:
        ax.plot(ionosonde_df.index, ionosonde_df['hmE'], color=hmE_color,
                label=f"{label} hmE", linewidth=2)


def overlay_austin_ionosonde(ax, csv_path='data/CSVs/2024-04-08_AU930_AustinTX_Ionosonde_ManualScaled.csv',
                             overlay_hmF2=True, overlay_hmE=True):
    """
    Overlays Austin, TX ionosonde data (hmF2 and hmE) on an existing axes.

    This is a convenience wrapper around overlay_ionosonde() with Austin-specific defaults.

    Arguments:
    ax : matplotlib.axes.Axes
        The axes object to plot on.
    csv_path : str, optional
        Path to the ionosonde CSV file. Default is 'data/CSVs/2024-04-08_AU930_AustinTX_Ionosonde_ManualScaled.csv'.
    overlay_hmF2 : bool, optional
        If True, overlay hmF2 data. Default is True.
    overlay_hmE : bool, optional
        If True, overlay hmE data. Default is True.
    """
    overlay_ionosonde(ax, csv_path=csv_path, overlay_hmF2=overlay_hmF2,
                     overlay_hmE=overlay_hmE, label='Austin Ionosonde',
                     hmF2_color='purple', hmE_color='brown')


def overlay_tdoa_csv(ax, csv_path_period=None, csv_path_autocorr=None,
                     label_period='Manual Period Analysis',
                     label_autocorr='Manual Autocorrelation Analysis',
                     color_period='tab:blue', linestyle_period='dotted',
                     color_autocorr='tab:orange', linestyle_autocorr='dashdot'):
    """
    Overlays manual TDOA analysis data from CSV files on an existing axes.

    This function loads verified manual TDOA measurements from CSV files
    (either manual period analysis or manual autocorrelation analysis)
    and plots the pre-calculated layer heights.

    Arguments:
    ax : matplotlib.axes.Axes
        The axes object to plot on.
    csv_path_period : str, optional
        Path to the manual period analysis CSV file. If None, period data is not plotted.
    csv_path_autocorr : str, optional
        Path to the manual autocorrelation analysis CSV file. If None, autocorr data is not plotted.
    label_period : str, optional
        Label for period analysis data. Default is 'Manual Period Analysis'.
    label_autocorr : str, optional
        Label for autocorrelation analysis data. Default is 'Manual Autocorrelation Analysis'.
    color_period : str, optional
        Color for period analysis line. Default is 'tab:blue'.
    linestyle_period : str, optional
        Linestyle for period analysis. Default is 'dotted'.
    color_autocorr : str, optional
        Color for autocorrelation analysis line. Default is 'tab:orange'.
    linestyle_autocorr : str, optional
        Linestyle for autocorrelation analysis. Default is 'dashdot'.
    """
    # Plot manual period analysis if CSV path provided
    if csv_path_period is not None:
        period_df = load_tdoa_csv(csv_path_period)
        if 'manualBeatNote_height_km' in period_df.columns:
            ax.plot(period_df.index, period_df['manualBeatNote_height_km'],
                    color=color_period, linestyle=linestyle_period,
                    label=label_period, linewidth=2)
            ax.scatter(period_df.index, period_df['manualBeatNote_height_km'],
                       color=color_period)

    # Plot manual autocorrelation analysis if CSV path provided
    if csv_path_autocorr is not None:
        autocorr_df = load_tdoa_csv(csv_path_autocorr)
        if 'manualBeatNote_height_km' in autocorr_df.columns:
            ax.plot(autocorr_df.index, autocorr_df['manualBeatNote_height_km'],
                    color=color_autocorr, linestyle=linestyle_autocorr,
                    label=label_autocorr, linewidth=2)
            ax.scatter(autocorr_df.index, autocorr_df['manualBeatNote_height_km'],
                       color=color_autocorr)


# ============================================================================
# Scatter Plot Comparison Functions
# ============================================================================

def align_and_resample_data(tdoa_df, ionosonde_df, resample_rule='1min', method='linear'):
    """
    Align TDOA data with ionosonde measurements using high-resolution ionosonde interpolation.

    This function resamples and interpolates the ionosonde data to a fine time resolution
    (1 minute by default), then matches each TDOA measurement to the nearest interpolated
    ionosonde value. This preserves the exact number of TDOA measurements while providing
    accurate ionosonde hmF2 values at the TDOA measurement times.

    Parameters:
    -----------
    tdoa_df : pd.DataFrame
        TDOA data with datetime index and height column
    ionosonde_df : pd.DataFrame
        Ionosonde data with datetime index and hmF2 column
    resample_rule : str, optional
        Resampling frequency for ionosonde interpolation (default: '1min')
    method : str, optional
        Interpolation method: 'linear', 'nearest', etc. (default: 'linear')

    Returns:
    --------
    aligned_df : pd.DataFrame
        DataFrame with columns: tdoa_height, hmF2, with aligned timestamps
        Contains exactly one row per original TDOA measurement
    tdoa_resampled : pd.DataFrame
        Copy of original TDOA data (for validation plotting)
    """
    # Resample ionosonde data to fine resolution (1 minute) and interpolate
    ionosonde_resampled = ionosonde_df[['hmF2']].resample(resample_rule).mean()
    if method == 'linear':
        ionosonde_resampled = ionosonde_resampled.interpolate(method='time')

    # Rename the TDOA height column to a standard name for merging
    tdoa_renamed = tdoa_df.copy()
    tdoa_renamed.columns = ['tdoa_height']

    # Drop any rows with NaN values in TDOA data (both NaN index and NaN values)
    # This is necessary because merge_asof requires non-null merge keys
    tdoa_renamed = tdoa_renamed.dropna()

    # Also drop rows where the index itself is NaT (Not a Time)
    tdoa_renamed = tdoa_renamed[tdoa_renamed.index.notna()]

    # Use merge_asof for time-tolerance matching
    # This matches each TDOA measurement to the nearest ionosonde value
    # Result will have exactly as many rows as TDOA measurements
    aligned_df = pd.merge_asof(
        tdoa_renamed.sort_index(),
        ionosonde_resampled.sort_index(),
        left_index=True,
        right_index=True,
        direction='nearest',
        tolerance=pd.Timedelta('2min')  # Allow up to 2 minutes difference
    )

    # Drop rows with NaN values (shouldn't be any with good ionosonde coverage)
    aligned_df = aligned_df.dropna()

    # Return original TDOA data for validation plotting
    return aligned_df, tdoa_df


def plot_resampling_validation(original_df, resampled_df, label, color,
                                 ax=None, savefig=None):
    """
    Plot original and resampled data to validate resampling fidelity.

    Parameters:
    -----------
    original_df : pd.DataFrame
        Original TDOA data with datetime index
    resampled_df : pd.DataFrame
        Resampled TDOA data
    label : str
        Label for the dataset
    color : str
        Color for plotting
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. If None, creates new figure.
    savefig : str, optional
        Path to save figure

    Returns:
    --------
    ax : matplotlib.axes.Axes
        The axes object
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(12, 6))

    # Plot original data
    ax.scatter(original_df.index, original_df.iloc[:, 0],
               color=color, alpha=0.5, s=50, label=f'{label} (Original)', zorder=2)

    # Plot resampled data
    ax.plot(resampled_df.index, resampled_df.iloc[:, 0],
            color=color, linewidth=2, linestyle='--',
            label=f'{label} (Resampled 5min)', zorder=3)
    ax.scatter(resampled_df.index, resampled_df.iloc[:, 0],
               color=color, s=100, marker='x', linewidth=2, zorder=4)

    ax.set_xlabel('UTC Time')
    ax.set_ylabel('Height (km)')
    ax.set_title(f'Resampling Validation: {label}')
    ax.legend(loc='best')
    ax.grid(True, alpha=0.3)

    # Format x-axis
    ax.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M'))
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=45, ha='right')

    if savefig:
        plt.savefig(savefig, dpi=300, bbox_inches='tight')
        print(f"Resampling validation plot saved to: {savefig}")

    return ax


def plot_scatter_comparison(tdoa_height, ionosonde_hmf2, label, color,
                             ax=None, marker='o', s=50, alpha=0.7,
                             show_1to1=True, show_stats=True, edgecolor='black', zorder=None,
                             linewidths=None):
    """
    Create scatter plot comparing TDOA heights with ionosonde hmF2.

    Parameters:
    -----------
    tdoa_height : pd.Series or np.array
        TDOA-derived heights
    ionosonde_hmf2 : pd.Series or np.array
        Ionosonde hmF2 values
    label : str
        Label for the dataset
    color : str
        Color for the scatter points
    ax : matplotlib.axes.Axes, optional
        Axes to plot on. If None, creates new figure.
    marker : str, optional
        Marker style (default: 'o')
    s : float, optional
        Marker size (default: 50)
    alpha : float, optional
        Transparency (default: 0.7)
    show_1to1 : bool, optional
        Show 1:1 reference line (default: True)
    show_stats : bool, optional
        Show correlation statistics (default: True)
    edgecolor : str, optional
        Edge color for the scatter points (default: 'black')
    zorder : int, optional
        Drawing order for the scatter points (default: None, which uses matplotlib default)
    linewidths : float, optional
        Line width for line-based markers like 'x', '+', etc. (default: None, uses 0.5 for edgecolors)

    Returns:
    --------
    ax : matplotlib.axes.Axes
        The axes object
    stats : dict
        Dictionary with correlation statistics
    """
    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 8))

    # Create scatter plot
    scatter_kwargs = {
        'color': color,
        'marker': marker,
        's': s,
        'alpha': alpha,
        'label': label,
        'edgecolors': edgecolor
    }

    # Use linewidths if specified, otherwise use linewidth
    # Note: linewidth and linewidths are aliases in matplotlib, can't use both
    if linewidths is not None:
        scatter_kwargs['linewidths'] = linewidths
    else:
        scatter_kwargs['linewidth'] = 0.5

    if zorder is not None:
        scatter_kwargs['zorder'] = zorder

    ax.scatter(ionosonde_hmf2, tdoa_height, **scatter_kwargs)

    # Calculate statistics
    valid_mask = ~np.isnan(tdoa_height) & ~np.isnan(ionosonde_hmf2)
    tdoa_valid = np.array(tdoa_height)[valid_mask]
    iono_valid = np.array(ionosonde_hmf2)[valid_mask]

    stats = {}
    if len(tdoa_valid) > 0:
        correlation = np.corrcoef(iono_valid, tdoa_valid)[0, 1]
        rmse = np.sqrt(np.mean((tdoa_valid - iono_valid)**2))
        bias = np.mean(tdoa_valid - iono_valid)

        # Calculate percent errors using ionosonde as reference (accepted values)
        rmse_percent = (rmse / np.mean(iono_valid)) * 100
        bias_percent = (bias / np.mean(iono_valid)) * 100

        stats = {
            'correlation': correlation,
            'rmse': rmse,
            'bias': bias,
            'rmse_percent': rmse_percent,
            'bias_percent': bias_percent,
            'n_points': len(tdoa_valid)
        }

    # Add 1:1 reference line
    if show_1to1:
        lims = [
            np.min([ax.get_xlim()[0], ax.get_ylim()[0]]),
            np.max([ax.get_xlim()[1], ax.get_ylim()[1]]),
        ]
        ax.plot(lims, lims, 'k--', alpha=0.5, linewidth=2, label='1:1 Line', zorder=1)

    return ax, stats


def create_scatter_plot_figure(datasets, ionosonde_df,
                                 individual_plots=True, combined_plot=True,
                                 output_dir=None):
    """
    Create scatter plots comparing multiple TDOA datasets with ionosonde hmF2.

    Parameters:
    -----------
    datasets : list of dict
        List of dataset dictionaries with keys:
        - 'df': DataFrame with TDOA data
        - 'label': Label for the dataset
        - 'color': Color for plotting
        - 'marker': Marker style
    ionosonde_df : pd.DataFrame
        Ionosonde data with hmF2 column
    individual_plots : bool, optional
        Create individual scatter plots for each dataset (default: True)
    combined_plot : bool, optional
        Create combined scatter plot with all datasets (default: True)
    output_dir : str, optional
        Directory to save plots. If None, displays only.

    Returns:
    --------
    results : dict
        Dictionary with statistics for each dataset
    """
    results = {}

    # Create individual plots
    if individual_plots:
        for dataset in datasets:
            fig, ax = plt.subplots(figsize=(8, 8))

            # Align and resample data
            aligned_df, resampled_df = align_and_resample_data(
                dataset['df'], ionosonde_df
            )

            # Create scatter plot
            ax, stats = plot_scatter_comparison(
                aligned_df['tdoa_height'],
                aligned_df['hmF2'],
                label=dataset['label'],
                color=dataset['color'],
                marker=dataset.get('marker', 'o'),
                ax=ax
            )

            # Add labels and title
            ax.set_xlabel('Austin Ionosonde hmF2 (km)', fontweight='bold')
            ax.set_ylabel('HF TDOA Height (km)', fontweight='bold')
            ax.set_title(f"TDOA vs Ionosonde hmF2\n{dataset['label']}",
                        fontweight='bold')
            ax.legend(loc='upper right', fontsize='small')
            ax.grid(True, alpha=0.3)

            # Set axis limits and aspect ratio
            ax.set_xlim(200, 400)
            ax.set_ylim(200, 400)
            ax.set_aspect('equal')

            # Add statistics text
            if stats:
                stats_text = (f"n = {stats['n_points']}\n"
                             f"r = {stats['correlation']:.3f}\n"
                             f"RMSE = {stats['rmse']:.1f} km\n"
                             f"Bias = {stats['bias']:.1f} km")
                ax.text(0.05, 0.95, stats_text, transform=ax.transAxes,
                       verticalalignment='top', fontsize=12,
                       bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))

            # Save or show
            if output_dir:
                os.makedirs(output_dir, exist_ok=True)
                filename = f"scatter_{dataset['label'].replace(' ', '_').replace('/', '-')}.jpg"
                filepath = os.path.join(output_dir, filename)
                plt.savefig(filepath, dpi=300, bbox_inches='tight')
                print(f"Saved: {filepath}")
            else:
                plt.show()

            results[dataset['label']] = stats

    # Create combined plot with statistics table
    if combined_plot:
        # Create figure with two panels: (a) scatter plot, (b) statistics table
        fig = plt.figure(figsize=(20, 8))
        gs = fig.add_gridspec(1, 2, width_ratios=[1, 1.2], wspace=0.3)

        # Panel (a): Scatter plot
        ax_scatter = fig.add_subplot(gs[0, 0])

        for idx, dataset in enumerate(datasets):
            # Align and resample data
            aligned_df, resampled_df = align_and_resample_data(
                dataset['df'], ionosonde_df
            )

            # Create scatter plot
            ax_scatter, stats = plot_scatter_comparison(
                aligned_df['tdoa_height'],
                aligned_df['hmF2'],
                label=dataset['label'],
                color=dataset['color'],
                marker=dataset.get('marker', 'o'),
                ax=ax_scatter,
                show_1to1=False  # Don't show 1:1 line
            )

            results[dataset['label']] = stats

        # Add labels and title for scatter plot
        ax_scatter.set_xlabel('Austin Ionosonde hmF2 (km)', fontweight='bold')
        ax_scatter.set_ylabel('HF TDOA Height (km)', fontweight='bold')
        ax_scatter.set_title('(a) HF TDOA vs Ionosonde hmF2',
                            fontweight='bold', fontsize=16, loc='left')
        ax_scatter.legend(loc='upper right', fontsize='small')
        ax_scatter.grid(True, alpha=0.3)

        # Set axis limits and aspect ratio
        ax_scatter.set_xlim(200, 400)
        ax_scatter.set_ylim(200, 400)
        ax_scatter.set_aspect('equal')

        # Panel (b): Statistics table with symbols
        ax_table = fig.add_subplot(gs[0, 1])
        ax_table.axis('off')
        ax_table.set_title('(b) Correlation Statistics',
                          fontweight='bold', fontsize=16, loc='left', pad=20)

        # Prepare table data (we'll add symbols manually)
        table_data = []
        for dataset in datasets:
            label = dataset['label']
            if label in results and results[label]:
                stats = results[label]
                table_data.append([
                    '',  # Placeholder for symbol
                    label,
                    f"{stats['n_points']}",
                    f"{stats['correlation']:.3f}",
                    f"{stats['rmse_percent']:.1f}%",
                    f"{stats['bias_percent']:.1f}%"
                ])

        # Create table with extra column for symbols
        table = ax_table.table(
            cellText=table_data,
            colLabels=['', 'Dataset', 'n', 'r', 'RMSE (%)', 'Bias (%)'],
            cellLoc='left',
            loc='upper center',
            bbox=[0, 0, 1, 0.9]
        )

        # Style the table
        table.auto_set_font_size(False)
        table.set_fontsize(11)
        table.scale(1, 2.5)

        # Set custom column widths: [symbol, dataset, n, r, RMSE, Bias]
        col_widths = [0.08, 0.35, 0.10, 0.12, 0.15, 0.15]
        for row in range(len(datasets) + 1):  # +1 for header
            for col in range(6):
                cell = table[(row, col)]
                cell.set_width(col_widths[col])

        # Color header row
        for i in range(6):
            cell = table[(0, i)]
            cell.set_facecolor('#4472C4')
            cell.set_text_props(weight='bold', color='white')

        # Alternate row colors
        for idx, dataset in enumerate(datasets):
            row_idx = idx + 1
            for j in range(6):
                cell = table[(row_idx, j)]
                if row_idx % 2 == 0:
                    cell.set_facecolor('#E7E6E6')
                else:
                    cell.set_facecolor('white')

        # Add symbols using scatter plot in table coordinates
        # The table is positioned at bbox=[0, 0, 1, 0.9], so we calculate relative positions
        num_rows = len(datasets) + 1  # +1 for header
        row_height = 0.9 / num_rows

        for idx, dataset in enumerate(datasets):
            row_idx = idx + 1
            # Calculate y position (from top, accounting for header)
            y_pos = 0.9 - (row_idx + 0.5) * row_height
            # X position in first column (centered)
            x_pos = 0.04  # Small offset for first column center

            # Add marker
            ax_table.plot(
                x_pos, y_pos,
                marker=dataset.get('marker', 'o'),
                color=dataset['color'],
                markersize=12,
                markeredgecolor='black',
                markeredgewidth=0.5,
                linestyle='',
                transform=ax_table.transAxes,
                clip_on=False
            )

        # Save or show
        if output_dir:
            os.makedirs(output_dir, exist_ok=True)
            filepath = os.path.join(output_dir, 'scatter_all_datasets.jpg')
            plt.savefig(filepath, dpi=300, bbox_inches='tight')
            print(f"Saved: {filepath}")
        else:
            plt.show()

    return results


def create_manual_vs_automated_scatter_figure(manual_datasets, automated_datasets, ionosonde_df, filepath=None):
    """
    Create a 2-row figure comparing manual and automated TDOA analysis methods.

    Layout:
    - Row 1: (a) Manual analysis scatter plot with integrated statistics table
    - Row 2: (b) Automated analysis scatter plot with integrated statistics table

    Each row contains:
    - Left: Scatter plot of TDOA heights vs ionosonde hmF2
    - Right: Statistics table with correlation metrics (n, r, RMSE %, Bias %)

    Parameters
    ----------
    manual_datasets : list of dict
        List of manual TDOA datasets, each with 'df', 'label', 'color', and 'marker'
    automated_datasets : list of dict
        List of automated TDOA datasets, each with 'df', 'label', 'color', and 'marker'
    ionosonde_df : pandas.DataFrame
        Ionosonde data with hmF2 column and datetime index
    filepath : str, optional
        Full path to save the output figure (including filename)

    Returns
    -------
    dict
        Dictionary with 'manual' and 'automated' keys, each containing results for datasets
    """
    import matplotlib.pyplot as plt
    from matplotlib.gridspec import GridSpec

    # Table font size configuration
    TABLE_FONTSIZE = 14

    results = {'manual': {}, 'automated': {}}

    # Create figure with 2x2 grid
    fig = plt.figure(figsize=(24, 16))
    gs = GridSpec(2, 2, figure=fig, width_ratios=[1, 1], height_ratios=[1, 1],
                  hspace=0.3, wspace=0.05)

    # ========================================================================
    # Top Row: Manual Analysis
    # ========================================================================

    # Panel (a): Manual scatter plot
    ax_manual_scatter = fig.add_subplot(gs[0, 0])

    for idx, dataset in enumerate(manual_datasets):
        # Align and resample data
        aligned_df, resampled_df = align_and_resample_data(
            dataset['df'], ionosonde_df
        )

        # Create scatter plot
        ax_manual_scatter, stats = plot_scatter_comparison(
            aligned_df['tdoa_height'],
            aligned_df['hmF2'],
            label=dataset['label'],
            color=dataset['color'],
            marker=dataset.get('marker', 'o'),
            edgecolor=dataset.get('edgecolor', 'black'),
            zorder=dataset.get('zorder', None),
            linewidths=dataset.get('linewidths', None),
            ax=ax_manual_scatter,
            show_1to1=False
        )

        results['manual'][dataset['label']] = stats

    # Style manual scatter plot
    ax_manual_scatter.set_xlabel('Austin Ionosonde hmF2 (km)', fontweight='bold', fontsize=14)
    ax_manual_scatter.set_ylabel('HF TDOA Height (km)', fontweight='bold', fontsize=14)
    # Use text instead of set_title for better positioning control
    ax_manual_scatter.text(-0.15, 1.05, '(a) Manual TDOA Heights vs Austin Ionosonde hmF2',
                          transform=ax_manual_scatter.transAxes,
                          fontweight='bold', fontsize=24, va='bottom', ha='left')
    # ax_manual_scatter.legend(loc='upper right', fontsize='small')  # Legend removed
    ax_manual_scatter.grid(True, alpha=0.3)
    ax_manual_scatter.set_xlim(200, 350)
    ax_manual_scatter.set_ylim(200, 350)
    ax_manual_scatter.set_aspect('equal')

    # Panel (b): Manual statistics table
    ax_manual_table = fig.add_subplot(gs[0, 1])
    ax_manual_table.axis('off')
    # ax_manual_table.set_title('(b) Manual Analysis: Correlation Statistics',
    #                           fontweight='bold', fontsize=16, loc='left', pad=20)

    # Prepare manual table data
    manual_table_data = []
    for dataset in manual_datasets:
        label = dataset['label']
        if label in results['manual'] and results['manual'][label]:
            stats = results['manual'][label]
            manual_table_data.append([
                '',  # Placeholder for symbol
                label,
                f"{stats['n_points']}",
                f"{stats['correlation']:.3f}",
                f"{stats['rmse_percent']:.1f}%",
                f"{stats['bias_percent']:.1f}%"
            ])

    # Create manual table
    manual_table = ax_manual_table.table(
        cellText=manual_table_data,
        colLabels=['', 'Dataset', 'n', 'r', 'RMSE (%)', 'Bias (%)'],
        cellLoc='left',
        loc='upper left',
        bbox=[-0.2, 0.05, 1.15, 0.9]
    )

    # Style manual table
    manual_table.auto_set_font_size(False)
    manual_table.set_fontsize(TABLE_FONTSIZE)
    manual_table.scale(1, 2.5)

    # Set column widths
    col_widths = [0.08, 0.45, 0.10, 0.12, 0.15, 0.15]
    for row in range(len(manual_datasets) + 1):
        for col in range(6):
            cell = manual_table[(row, col)]
            cell.set_width(col_widths[col])

    # Color header
    for i in range(6):
        cell = manual_table[(0, i)]
        cell.set_facecolor('#4472C4')
        cell.set_text_props(weight='bold', color='white', fontsize=TABLE_FONTSIZE + 2)

    # Alternate row colors
    for idx, dataset in enumerate(manual_datasets):
        row_idx = idx + 1
        for j in range(6):
            cell = manual_table[(row_idx, j)]
            if row_idx % 2 == 0:
                cell.set_facecolor('#E7E6E6')
            else:
                cell.set_facecolor('white')

    # Add symbols to manual table
    num_rows_manual = len(manual_datasets) + 1
    row_height_manual = 0.9 / num_rows_manual
    table_bottom = 0.05  # Match the bbox bottom position
    for idx, dataset in enumerate(manual_datasets):
        row_idx = idx + 1
        y_pos = (table_bottom + 0.9) - (row_idx + 0.5) * row_height_manual
        x_pos = -0.16  # Adjusted for shifted table (was 0.04)

        # Use linewidths for line-based markers if specified
        markeredgewidth = dataset.get('linewidths', 0.5) if dataset.get('marker') in ['x', '+', '*', '1', '2', '3', '4'] else 0.5

        ax_manual_table.plot(
            x_pos, y_pos,
            marker=dataset.get('marker', 'o'),
            color=dataset['color'],
            markersize=12,
            markeredgecolor=dataset.get('edgecolor', 'black'),
            markeredgewidth=markeredgewidth,
            linestyle='',
            transform=ax_manual_table.transAxes,
            clip_on=False
        )

    # ========================================================================
    # Bottom Row: Automated Analysis
    # ========================================================================

    # Panel (c): Automated scatter plot
    ax_auto_scatter = fig.add_subplot(gs[1, 0])

    for idx, dataset in enumerate(automated_datasets):
        # Align and resample data
        aligned_df, resampled_df = align_and_resample_data(
            dataset['df'], ionosonde_df
        )

        # Create scatter plot
        ax_auto_scatter, stats = plot_scatter_comparison(
            aligned_df['tdoa_height'],
            aligned_df['hmF2'],
            label=dataset['label'],
            color=dataset['color'],
            marker=dataset.get('marker', 'o'),
            edgecolor=dataset.get('edgecolor', 'black'),
            zorder=dataset.get('zorder', None),
            linewidths=dataset.get('linewidths', None),
            ax=ax_auto_scatter,
            show_1to1=False
        )

        results['automated'][dataset['label']] = stats

    # Style automated scatter plot
    ax_auto_scatter.set_xlabel('Austin Ionosonde hmF2 (km)', fontweight='bold', fontsize=14)
    ax_auto_scatter.set_ylabel('HF TDOA Height (km)', fontweight='bold', fontsize=14)
    # Use text instead of set_title for better positioning control
    ax_auto_scatter.text(-0.15, 1.05, '(b) Automated TDOA Heights vs Austin Ionosonde hmF2',
                        transform=ax_auto_scatter.transAxes,
                        fontweight='bold', fontsize=24, va='bottom', ha='left')
    # ax_auto_scatter.legend(loc='upper right', fontsize='small')  # Legend removed
    ax_auto_scatter.grid(True, alpha=0.3)
    ax_auto_scatter.set_xlim(200, 350)
    ax_auto_scatter.set_ylim(200, 350)
    ax_auto_scatter.set_aspect('equal')

    # Panel (d): Automated statistics table
    ax_auto_table = fig.add_subplot(gs[1, 1])
    ax_auto_table.axis('off')
    # ax_auto_table.set_title('(d) Automated Analysis: Correlation Statistics',
    #                         fontweight='bold', fontsize=16, loc='left', pad=20)

    # Prepare automated table data
    auto_table_data = []
    for dataset in automated_datasets:
        label = dataset['label']
        if label in results['automated'] and results['automated'][label]:
            stats = results['automated'][label]
            auto_table_data.append([
                '',  # Placeholder for symbol
                label,
                f"{stats['n_points']}",
                f"{stats['correlation']:.3f}",
                f"{stats['rmse_percent']:.1f}%",
                f"{stats['bias_percent']:.1f}%"
            ])

    # Create automated table
    auto_table = ax_auto_table.table(
        cellText=auto_table_data,
        colLabels=['', 'Dataset', 'n', 'r', 'RMSE (%)', 'Bias (%)'],
        cellLoc='left',
        loc='upper left',
        bbox=[-0.2, 0.05, 1.15, 0.9]
    )

    # Style automated table
    auto_table.auto_set_font_size(False)
    auto_table.set_fontsize(TABLE_FONTSIZE)
    auto_table.scale(1, 2.5)

    # Set column widths
    for row in range(len(automated_datasets) + 1):
        for col in range(6):
            cell = auto_table[(row, col)]
            cell.set_width(col_widths[col])

    # Color header
    for i in range(6):
        cell = auto_table[(0, i)]
        cell.set_facecolor('#4472C4')
        cell.set_text_props(weight='bold', color='white', fontsize=TABLE_FONTSIZE + 2)

    # Alternate row colors
    for idx, dataset in enumerate(automated_datasets):
        row_idx = idx + 1
        for j in range(6):
            cell = auto_table[(row_idx, j)]
            if row_idx % 2 == 0:
                cell.set_facecolor('#E7E6E6')
            else:
                cell.set_facecolor('white')

    # Add symbols to automated table
    num_rows_auto = len(automated_datasets) + 1
    row_height_auto = 0.9 / num_rows_auto
    table_bottom = 0.05  # Match the bbox bottom position
    for idx, dataset in enumerate(automated_datasets):
        row_idx = idx + 1
        y_pos = (table_bottom + 0.9) - (row_idx + 0.5) * row_height_auto
        x_pos = -0.16  # Adjusted for shifted table (was 0.04)

        # Use linewidths for line-based markers if specified
        markeredgewidth = dataset.get('linewidths', 0.5) if dataset.get('marker') in ['x', '+', '*', '1', '2', '3', '4'] else 0.5

        ax_auto_table.plot(
            x_pos, y_pos,
            marker=dataset.get('marker', 'o'),
            color=dataset['color'],
            markersize=12,
            markeredgecolor=dataset.get('edgecolor', 'black'),
            markeredgewidth=markeredgewidth,
            linestyle='',
            transform=ax_auto_table.transAxes,
            clip_on=False
        )

    # Save or show
    if filepath:
        # Create parent directory if it doesn't exist
        parent_dir = os.path.dirname(filepath)
        if parent_dir:
            os.makedirs(parent_dir, exist_ok=True)
        plt.savefig(filepath, dpi=300, bbox_inches='tight')
        print(f"Saved: {filepath}")
    else:
        plt.show()

    return results
