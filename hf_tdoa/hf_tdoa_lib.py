"""
HF TDOA Analysis Library

This library provides functions for analyzing HF TDOA (Time Difference of Arrival) data
from WAV file recordings of chirp signals.
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
# Utility Functions for Reducing Code Duplication
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
    pbar = tqdm(wavlist, desc="Finding chirps", unit="file", disable=plot_correlation)

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

    chirps = pd.DataFrame({'utc':utcs, 'data':wav_data, 'x':chirps_x, 'y':chirps_y, 'file':file})
    chirps.attrs['path_info'] = path_info
    chirps.attrs['pfx'] = path_info.pfx  # Keep for backward compatibility
    chirps.attrs['sweep_rate'] = sweep_rate

    print(f"✓ Completed chirp detection: {len(chirps)} files processed, {len(chirps)*10} total chirps found\n")

    return chirps


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


def plot_chirp_fft(title, tlim, minfreq, maxfreq, wav_df, sweep_rate, env, tvec, X_psd, f):
    """
    Plots the FFT of a chirp signal along with its envelope and raw waveform.

    Arguments:
    title : str
        Title for the plot.
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
    """
    if maxfreq == 0:
        maxfreq = 10000
    flim = (0, maxfreq)

    nrows = 3
    nax   = 0
    fig = plt.figure(figsize=(15,15))

    nax += 1
    ax  = fig.add_subplot(nrows, 1, nax)
    if wav_df is not None:
        xx = wav_df.index
        yy = wav_df['x']
        ax.plot(xx, yy)
        ax.set_xlabel('Time [s]')
        ax.set_ylabel('Amplitude')
        ax.axvspan(*tlim, color='red', label='Selected Chirp')
        ax.legend(loc='upper right')

    nax += 1
    ax  = fig.add_subplot(nrows, 1, nax)

    dc_offset = 0
    if wav_df is not None:
        tf  = np.logical_and(wav_df.index >= tlim[0], wav_df.index <  tlim[1])
        wav_xx = wav_df.index[tf]
        wav_yy = wav_df['x'][tf]
        ax.plot(wav_xx, wav_yy, label='Raw WAV')
        dc_offset = np.mean(np.abs(wav_yy))

    env_xx = tvec
    env_yy = env + dc_offset
    ax.plot(env_xx, env_yy, lw=3, color='red', label='Filtered Envelope')
    ax.set_xlim(tlim)
    ax.set_xlabel('Time [s]')
    ax.set_ylabel('Amplitude')
    ax.legend(loc='upper right')
    ax.set_ylim(0, None)

    nax += 1
    ax = fig.add_subplot(nrows, 1, nax)
    ax.plot(f, X_psd)

    maximum_x, maximum_y, local_peaks_x, local_peaks_y = find_max(f, X_psd, minfreq, maxfreq)

    ax.plot([maximum_x, maximum_x], [0, maximum_y], linewidth=2)
    ax.plot([0, maximum_x], [maximum_y, maximum_y], linewidth=2, color='red')
    ax.plot([minfreq, minfreq], [0, maximum_y], linewidth=3, color='black')

    lbl = []
    lbl.append(f'$f_b = {maximum_x:0.2f}$ Hz')
    lbl.append(f'TDOA = {maximum_x/sweep_rate:0.2f} ms')
    ax.scatter(maximum_x, maximum_y, linewidths=5, label='\n'.join(lbl), color='green')
    ax.scatter(local_peaks_x, local_peaks_y, color='black')
    ax.set_xlim(flim)
    ax.set_ylim(0, 10*10**-5)
    ax.set_xlabel('Frequency [Hz]')
    ax.set_ylabel('$|X(f)|^2$')

    # Make x-labels with both beat frequency and TDOA values.
    xticks = ax.get_xticks()
    ax.set_xticks(xticks)
    xtls = []
    for xtk in xticks:
        tdoa = xtk/sweep_rate
        xtl = f'{xtk:0.1f}\n{tdoa:0.1f}'
        xtls.append(xtl)
    xtls[-1] = 'f [Hz]\nTDOA [ms]'
    ax.set_xticklabels(xtls)
    ax.set_xlabel('')

    ax.legend()

    _title = []
    _title.append(title)
    _title.append(f'Chirp Sweep Rate: {sweep_rate} Hz/ms')
    plt.suptitle('\n'.join(_title))

    plt.tight_layout()
    plt.show()
    plt.close(fig)


def find_TDOAs(wav_data, search_limits=None, filter_limts=None,
               mode_string='2F2-1F2', set_name=None, plot_fft=False, only_one=True, **overrides):
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
    only_one : bool, optional
        If True, processes only the first chirp in first WAV file. Used for debugging.
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
                unit="file", disable=only_one)

    for file_num, row in pbar:
        maxes = []
        wav_df = row['data']
        fpath  = row['file']
        bname  = os.path.basename(fpath)
        wav_env, wav_fs = filter(wav_df, filter_limts[1], filter_limts[0])

        for peak in wav_data['x'][file_num]:
            peak = peak/wav_fs
            tlim = (peak+search_limits[0], peak+search_limits[1])
            minfreq = search_limits[2]
            maxfreq = search_limits[3]

            X_psd, freq = chirp_fft(wav_env, tlim=tlim)

            maximum_x, maximum_y, local_peaks_x, local_peaks_y = find_max(freq, X_psd, minfreq, maxfreq)

            if plot_fft:
                plot_chirp_fft(title=bname, tlim=tlim,
                           minfreq=minfreq, maxfreq=maxfreq, wav_df=wav_df, sweep_rate=sweep_rate,
                           env=wav_env['x'], tvec=wav_env.index, X_psd=X_psd, f=freq)

            maxes.append(maximum_x)
            if only_one:
                break

        beat_arr = np.array(maxes).flatten()/sweep_rate  # TDOA in ms
        if only_one:
            break
        else:
            all_beats.append(beat_arr)
            # Handle case where beat_arr contains only NaN values
            if len(beat_arr) > 0 and not np.all(np.isnan(beat_arr)):
                mean_beats.append(np.nanmean(beat_arr))
            else:
                mean_beats.append(np.nan)

            # Update progress bar with current file's mean TDOA
            if len(beat_arr) > 0 and not np.all(np.isnan(beat_arr)):
                pbar.set_postfix({'Mean TDOA': f'{np.nanmean(beat_arr):.2f} ms'})

    if not only_one:
        pbar.close()
        wav_data[set_name]           = all_beats
        wav_data[f'{set_name}_mean'] = mean_beats
        print(f"✓ Completed processing {mode_string}: {len(wav_data)} files processed\n")

    return wav_data


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


def build_tdoa_config(chirps, mode_strings=None, **mode_overrides):
    """
    Build a TDOA configuration dictionary for use with plot_TDOAs and plot_hmf2.

    This function creates a dictionary with mode-specific parameters and TDOA model coefficients.

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

    return tdoa_dct


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
                label=set_name,
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
        legend = ax.legend(handles=[line], loc=(xpos, ypos))
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


def load_ionosonde_data(csv_path):
    """
    Loads ionosonde data from a CSV file.

    Arguments:
    csv_path : str
        Path to the ionosonde CSV file. Expected columns: UTC, hmE, hmF2.

    Returns:
    ionosonde_df : pd.DataFrame
        DataFrame containing ionosonde data with datetime index.
    """
    ionosonde_df = pd.read_csv(csv_path, parse_dates=['UTC'])
    return ionosonde_df


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
        csv_path = os.path.join(package_dir, 'data', 'CSVs', '2024-04-08_Austin_TX_Ionosonde_hmE_hmF2.csv')

    ionosonde_df = load_ionosonde_data(csv_path)

    if overlay_hmF2:
        ax.plot(ionosonde_df['UTC'], ionosonde_df['hmF2'], color=hmF2_color,
                label=f"{label} hmF2", linewidth=2)

    if overlay_hmE:
        ax.plot(ionosonde_df['UTC'], ionosonde_df['hmE'], color=hmE_color,
                label=f"{label} hmE", linewidth=2)


def overlay_austin_ionosonde(ax, csv_path='data/CSVs/2024-04-08_Austin_TX_Ionosonde_hmE_hmF2.csv',
                             overlay_hmF2=True, overlay_hmE=True):
    """
    Overlays Austin, TX ionosonde data (hmF2 and hmE) on an existing axes.

    This is a convenience wrapper around overlay_ionosonde() with Austin-specific defaults.

    Arguments:
    ax : matplotlib.axes.Axes
        The axes object to plot on.
    csv_path : str, optional
        Path to the ionosonde CSV file. Default is 'data/CSVs/2024-04-08_Austin_TX_Ionosonde_hmE_hmF2.csv'.
    overlay_hmF2 : bool, optional
        If True, overlay hmF2 data. Default is True.
    overlay_hmE : bool, optional
        If True, overlay hmE data. Default is True.
    """
    overlay_ionosonde(ax, csv_path=csv_path, overlay_hmF2=overlay_hmF2,
                     overlay_hmE=overlay_hmE, label='Austin Ionosonde',
                     hmF2_color='purple', hmE_color='brown')


def load_tdoa_csv(csv_path, model_coeffs):
    """
    Loads TDOA data from a CSV file and applies the HF TDOA model to convert to layer heights.

    The CSV file should have columns:
    - utc: timestamp in ISO format (e.g., '2024-04-08 14:13:00')
    - manualBeatNote_TDOA_ms: Manual Period Analysis TDOA values in milliseconds
    - autoCorrelation_TDOA_ms: Auto-Correlated Analysis TDOA values in milliseconds

    Arguments:
    csv_path : str
        Path to the CSV file containing TDOA data.
    model_coeffs : tuple
        Tuple of (slope, intercept) for the TDOA model.
        Layer height = slope * TDOA_ms + intercept

    Returns:
    tdoa_df : pd.DataFrame
        DataFrame containing:
        - UTC datetime index
        - manualBeatNote_TDOA_ms: Original TDOA values in ms
        - autoCorrelation_TDOA_ms: Original TDOA values in ms
        - manualBeatNote_height_km: Layer heights from manual analysis
        - autoCorrelation_height_km: Layer heights from auto-correlation
    """
    # Load CSV with datetime parsing
    tdoa_df = pd.read_csv(csv_path, parse_dates=['utc'])
    tdoa_df = tdoa_df.set_index('utc')

    # Apply TDOA model to convert TDOA to layer heights
    slope, intercept = model_coeffs
    tdoa_df['manualBeatNote_height_km'] = slope * tdoa_df['manualBeatNote_TDOA_ms'] + intercept
    tdoa_df['autoCorrelation_height_km'] = slope * tdoa_df['autoCorrelation_TDOA_ms'] + intercept

    return tdoa_df


def overlay_tdoa_csv(ax, csv_path, model_coeffs,
                     overlay_manual=True, overlay_autocorr=True,
                     manual_label='Manual Period Analysis',
                     autocorr_label='Auto-Correlated Analysis',
                     manual_color='tab:blue', manual_linestyle='dotted',
                     autocorr_color='tab:orange', autocorr_linestyle='dashdot'):
    """
    Overlays TDOA data from CSV files on an existing axes.

    This function loads TDOA measurements from CSV files and applies the HF TDOA model
    to convert them to layer heights before plotting.

    Arguments:
    ax : matplotlib.axes.Axes
        The axes object to plot on.
    csv_path : str
        Path to the CSV file containing TDOA data.
    model_coeffs : tuple
        Tuple of (slope, intercept) for the TDOA model.
    overlay_manual : bool, optional
        If True, overlay manual period analysis data. Default is True.
    overlay_autocorr : bool, optional
        If True, overlay auto-correlated analysis data. Default is True.
    manual_label : str, optional
        Label for manual analysis data. Default is 'Manual Period Analysis'.
    autocorr_label : str, optional
        Label for auto-correlation data. Default is 'Auto-Correlated Analysis'.
    manual_color : str, optional
        Color for manual analysis line. Default is 'tab:blue'.
    manual_linestyle : str, optional
        Linestyle for manual analysis. Default is 'dotted'.
    autocorr_color : str, optional
        Color for auto-correlation line. Default is 'tab:orange'.
    autocorr_linestyle : str, optional
        Linestyle for auto-correlation. Default is 'dashdot'.
    """
    tdoa_df = load_tdoa_csv(csv_path, model_coeffs)

    if overlay_manual:
        ax.plot(tdoa_df.index, tdoa_df['manualBeatNote_height_km'],
                color=manual_color, linestyle=manual_linestyle,
                label=manual_label, linewidth=2)
        ax.scatter(tdoa_df.index, tdoa_df['manualBeatNote_height_km'],
                   color=manual_color)

    if overlay_autocorr:
        ax.plot(tdoa_df.index, tdoa_df['autoCorrelation_height_km'],
                color=autocorr_color, linestyle=autocorr_linestyle,
                label=autocorr_label, linewidth=2)
        ax.scatter(tdoa_df.index, tdoa_df['autoCorrelation_height_km'],
                   color=autocorr_color)


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
    show_date : bool, optional
        If True, show date in the title. Default is True.
    legend_loc : str, optional
        Legend location. Default is 'best'. Can be 'upper left', 'upper right', 'lower left',
        'lower right', 'center', 'best', etc.
    legend_fontsize : int or str, optional
        Legend font size. If None, uses matplotlib default. Can be int or 'small', 'medium', 'large', etc.
    """
    times = chirps['utc']

    # Plot TDOA-derived layer heights
    for set_name, params in tdoa_dct.items():
        TDOAs = chirps[f'{set_name}_mean']
        coefs = params['model_coeffs']
        layer_heights = (coefs[0] * TDOAs) + coefs[1]

        ax.plot(times, layer_heights,
                label=set_name,
                marker=params.get('marker', 'o'),
                linestyle=params.get('linestyle'),
                linewidth=params.get('linewidth'),
                color=params.get('color'))

    ax.set_ylim(ylim)

    _format_datetime_axis(ax)
    ax.set_ylabel('Layer Height [km]')
    ax.set_xlabel('Time UTC')

    # Overlay ionosonde data if requested
    if (ionosonde_dct is not None) and (ionosonde_dct is not False):
        if ionosonde_dct is True:
            ionosonde_dct = {}
        overlay_ionosonde(ax, **ionosonde_dct)

    # Overlay TDOA CSV data if requested
    if (tdoa_csv_dct is not None) and (tdoa_csv_dct is not False):
        if tdoa_csv_dct is True:
            tdoa_csv_dct = {}
        overlay_tdoa_csv(ax, **tdoa_csv_dct)

    # Add solar elevation and/or eclipse obscuration overlays if requested
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
                    solar_ts.overlaySolarElevation(ax)

                if overlay_eclipse:
                    solar_ts.overlayEclipse(ax)

            except (ImportError, AttributeError) as e:
                print(f'WARNING: Cannot overlay solar data. Error: {e}')

    ax.legend(loc=legend_loc, fontsize=legend_fontsize)

    # Use path_info from attrs
    path_info = chirps.attrs['path_info']

    # Show date if requested
    date = times.iloc[0] if show_date else None
    title_from_pfx(ax, path_info, date, center_title='Ionospheric Layer Height')


def plot_hmf2(chirps, tdoa_dct, ylim=(75,450),
              solar_lat=None, solar_lon=None, solar_start=None, solar_end=None,
              overlay_solar_elevation=False, overlay_eclipse=False,
              ionosonde_dct=None, tdoa_csv_dct=None, savefig=None):
    """
    Plots layer heights derived from TDOAs and compares with ionosonde measurements.

    Arguments:
    chirps : pd.DataFrame
        DataFrame containing chirp data and TDOA measurements.
    tdoa_dct : dict
        Dictionary containing TDOA set configurations with model coefficients and plotting parameters.
    ylim : tuple, optional
        Y-axis limits for the plot. Default is (75, 450).
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
            - csv_path: Path to TDOA CSV file (REQUIRED)
            - model_coeffs: Tuple of (slope, intercept) for TDOA model (REQUIRED)
            - overlay_manual: If True, overlay manual period analysis (default: True)
            - overlay_autocorr: If True, overlay auto-correlated analysis (default: True)
            - manual_label: Label for manual analysis (default: 'Manual Period Analysis')
            - autocorr_label: Label for auto-correlation (default: 'Auto-Correlated Analysis')
            - manual_color: Color for manual line (default: 'tab:blue')
            - manual_linestyle: Linestyle for manual (default: 'dotted')
            - autocorr_color: Color for autocorr line (default: 'tab:orange')
            - autocorr_linestyle: Linestyle for autocorr (default: 'dashdot')
        Example: {'csv_path': 'data/CSVs/file.csv', 'model_coeffs': (142, 36.1)}
        Default is None.
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
        show_date=True
    )

    fig.autofmt_xdate()
    plt.tight_layout()

    # Save figure if filename provided
    if savefig is not None:
        plt.savefig(savefig, dpi=300, bbox_inches='tight', format='jpeg', pil_kwargs={'quality': 95})
        print(f"Figure saved to: {savefig}")

    plt.show()
    plt.close(fig)


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
            show_date=True
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
    ax_tdoa.text(-0.08, 1.075, subplot_labels[0], transform=ax_tdoa.transAxes,
                fontsize=30, fontweight='bold', va='top', ha='left')

    ax_hmf2.text(-0.08, 1.075, subplot_labels[1], transform=ax_hmf2.transAxes,
                fontsize=30, fontweight='bold', va='top', ha='left')

    if image_panel is not None:
        ax_img.text(-0.08, 1.075, subplot_labels[2], transform=ax_img.transAxes,
                    fontsize=30, fontweight='bold', va='top', ha='left')

    # Save figure if filename provided
    if savefig is not None:
        plt.savefig(savefig, dpi=300, bbox_inches='tight', format='jpeg', pil_kwargs={'quality': 95})
        print(f"Figure saved to: {savefig}")

    plt.show()
    plt.close(fig)
