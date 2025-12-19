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

    Arguments:
    wavlist : list
        List of paths to WAV files to be analyzed.
    template : str
        Path to the WAV file containing the chirp template.
    sweep_rate : float
        Sweep rate in Hz/ms. Needed for TDOA calculations.
    pfx : str, optional
        String that identifies the transmit and receive station. If None, the pfx will
        be derived from the first name in the wavlist. This is stored as chirps.attrs['pfx'].
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
    """
    if pfx is None:
        pfx = os.path.basename(wavlist[0])[17:-4]

    sample_env, sample_fs = load_wav(template)
    template = sample_env['x']

    # Here we are effectively taking time components
    utcs     = []
    chirps_x = []
    chirps_y = []
    wav_data = []

    for file in wavlist:
        bname   = os.path.basename(file)
        utc_str = bname[:13]
        utc     = datetime.datetime.strptime(utc_str, '%Y%m%d.%H%M')

        wav_env, wav_fs = load_wav(file)

        y = signal.correlate(wav_env['x'], template, mode='same')
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
        else:
            print(file, end = "\r")

    chirps = pd.DataFrame({'utc':utcs, 'data':wav_data, 'x':chirps_x, 'y':chirps_y, 'file':file})
    chirps.attrs['pfx'] = pfx
    chirps.attrs['sweep_rate'] = sweep_rate
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


def find_TDOAs(wav_data, search_limits, filter_limts,
               set_name='TDOAs', plot_fft=False, only_one=True):
    """
    Finds TDOAs in a set of WAV files given chirp start locations.

    TDOA [ms] = (Beat Frequency [Hz]) / (Sweep Rate [Hz/ms])

    Arguments:
    wav_data : pd.DataFrame
        DataFrame containing WAV file data and chirp start locations.
    search_limits : tuple, (start_offset, end_offset, min_freq, max_freq)
        This is used to range gate the chirp signal in time to look for an echo from a particular ionospheric region.
        start_offset : float
            Start time offset (in seconds) from chirp start to begin FFT analysis.
        end_offset : float
            End time offset (in seconds) from chirp start to end FFT analysis.
        min_freq : float
            Minimum frequency for FFT peak search.
        max_freq : float
            Maximum frequency for FFT peak search.
    filter_limts : tuple
        Frequency limits for bandpass filtering (low_pass_freq, high_pass_freq).
    set_name : str, optional
        Name of the column to store TDOA results. Default is 'TDOAs'.
    plot_fft : bool, optional
        If True, plots the FFT for each chirp analyzed. Default is False.
    only_one : bool, optional
        If True, processes only the first chirp in first WAV file. Used for debugging.

    Returns:
    wav_data : pd.DataFrame
        Updated DataFrame with TDOA results added.
    """
    sweep_rate = wav_data.attrs['sweep_rate'] # Sweep rate in Hz/ms

    all_beats  = []
    mean_beats = []
    for file_num, row in wav_data.iterrows():
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
            print(fpath, end = "\r")
            if only_one:
                break

        beat_arr = np.array(maxes).flatten()/sweep_rate  # TDOA in ms
        if only_one:
            break
        else:
            all_beats.append(beat_arr)
            mean_beats.append(np.nanmean(beat_arr))

    if not only_one:
        wav_data[set_name]           = all_beats
        wav_data[f'{set_name}_mean'] = mean_beats
    return wav_data


def title_from_pfx(ax, pfx, date=None):
    """
    Creates a formatted title for plots based on the prefix string.

    Arguments:
    ax : matplotlib.axes.Axes
        The axes object to add the title to.
    pfx : str
        Prefix string containing TX/RX station information.
    date : datetime.datetime, optional
        Date to display in the center title.
    """
    pfx = pfx.replace('_', '-').split('-')
    tx_call = pfx[1]
    tx_grid = pfx[2]

    rx_call = pfx[4]
    rx_grid = pfx[5]

    range_km = float(pfx[6].lower().replace('km', ''))
    band     = int(pfx[7].replace('m', ''))

    band_dct = {}
    band_dct[80]  = '3.5 MHz'
    band_dct[60]  = '5 MHz'
    band_dct[40]  = '7 MHz'
    band_dct[20]  = '14 MHz'
    band_dct[15]  = '21 MHz'
    band_dct[10]  = '28 MHz'
    band_str = band_dct.get(band, f'{band} m')

    title = f'TX: {tx_call} ({tx_grid})\nRX: {rx_call} ({rx_grid})'
    ax.set_title(title, loc='left')

    title = f'Ground Range: {range_km} km\nBand: {band_str}'
    ax.set_title(title, loc='right')

    if date is not None:
        date_str = date.strftime('%Y %b %d')
        ax.set_title(date_str, loc='center', fontsize=34)


def plot_TDOAs(chirps, tdoa_dct, ylim=(0,3)):
    """
    Plots TDOA measurements over time for multiple propagation modes.

    Arguments:
    chirps : pd.DataFrame
        DataFrame containing chirp data and TDOA measurements.
    tdoa_dct : dict
        Dictionary containing TDOA set configurations with plotting parameters.
    ylim : tuple, optional
        Y-axis limits for the plot. Default is (0, 3).
    """
    times = chirps['utc']
    fig   = plt.figure(figsize=(16, 9))
    ax    = fig.add_subplot(1, 1, 1)

    lgnds = []
    for set_name, params in tdoa_dct.items():
        TDOAs = chirps[f'{set_name}_mean']
        line, = ax.plot(times, TDOAs,
                label=set_name,
                marker = params.get('marker', 'o'),
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
    ylim = ax.get_ylim()
    for lgnd in lgnds:
        line = lgnd['line']
        xpos = lgnd['xpos']
        ypos = lgnd['ypos']
        ypos = (ypos-ylim[0]) / (ylim[1]-ylim[0]) # Convert to fraction of axis height
        ypos = ypos - 0.06 # Nudge legend down a bit
        lgnd = ax.legend(handles=[line], loc=(xpos, ypos))
        ax.add_artist(lgnd)

    myFmt = mdates.DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(myFmt)
    ax.xaxis.set_major_locator(mdates.MinuteLocator(interval=15))

    ax.set_ylabel('TDOA [ms]')
    ax.set_xlabel('Time UTC')
    fig.autofmt_xdate()

    pfx = chirps.attrs.get('pfx', '')
    title_from_pfx(ax, pfx, times[0])

    plt.tight_layout()
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


def overlay_austin_ionosonde(ax, csv_path='data/CSVs/2024-04-08_Austin_TX_Ionosonde_hmE_hmF2.csv'):
    """
    Overlays Austin, TX ionosonde data (hmF2 and hmE) on an existing axes.

    Arguments:
    ax : matplotlib.axes.Axes
        The axes object to plot on.
    csv_path : str, optional
        Path to the ionosonde CSV file. Default is 'data/CSVs/2024-04-08_Austin_TX_Ionosonde_hmE_hmF2.csv'.
    """
    ionosonde_df = load_ionosonde_data(csv_path)

    ax.plot(ionosonde_df['UTC'], ionosonde_df['hmF2'], color='purple',
            label="Austin Ionosonde hmF2", linewidth=2)

    ax.plot(ionosonde_df['UTC'], ionosonde_df['hmE'], color='brown',
            label="Austin Ionosonde hmE", linewidth=2)


def plot_hmf2(chirps, tdoa_dct, csv_path='data/CSVs/2024-04-08_Austin_TX_Ionosonde_hmE_hmF2.csv', ylim=(75,450)):
    """
    Plots layer heights derived from TDOAs and compares with ionosonde measurements.

    Arguments:
    chirps : pd.DataFrame
        DataFrame containing chirp data and TDOA measurements.
    tdoa_dct : dict
        Dictionary containing TDOA set configurations with model coefficients and plotting parameters.
    csv_path : str, optional
        Path to the ionosonde CSV file. Default is 'data/CSVs/2024-04-08_Austin_TX_Ionosonde_hmE_hmF2.csv'.
    ylim : tuple, optional
        Y-axis limits for the plot. Default is (75, 450).
    """
    times = chirps['utc']

    fig   = plt.figure(figsize=(16, 9))
    ax    = fig.add_subplot(1, 1, 1)

    for set_name, params in tdoa_dct.items():
        TDOAs = chirps[f'{set_name}_mean']
        coefs = params['model_coeffs']
        layer_heights = (coefs[0] * TDOAs) + coefs[1]

        line, = ax.plot(times, layer_heights,
                label=set_name,
                marker = params.get('marker', 'o'),
                linestyle=params.get('linestyle'),
                linewidth=params.get('linewidth'),
                color=params.get('color'))

    ax.set_ylim(ylim)

    myFmt = mdates.DateFormatter('%H:%M')
    ax.xaxis.set_major_formatter(myFmt)
    ax.xaxis.set_major_locator(mdates.MinuteLocator(interval=15))
    ax.set_ylabel('Layer Height [km]')
    ax.set_xlabel('Time UTC')

    overlay_austin_ionosonde(ax, csv_path)

    ax.legend()
    fig.autofmt_xdate()

    pfx = chirps.attrs.get('pfx', '')
    title_from_pfx(ax, pfx, times[0])

    plt.tight_layout()
    plt.show()
