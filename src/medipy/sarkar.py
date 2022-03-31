'''
sarkar.py
QUelle
'''

# R-R interval in ms
rr = (rr / sf) * 1000
rri = np.diff(rr)

# Interpolate and compute HR


def interp_cubic_spline(rri, sf_up=4):
    """
    Interpolate R-R intervals using cubic spline.
    Taken from the `hrv` python package by Rhenan Bartels.

    Parameters
    ----------
    rri : np.array
        R-R peak interval (in ms)
    sf_up : float
        Upsampling frequency.

    Returns
    -------
    rri_interp : np.array
        Upsampled/interpolated R-R peak interval array
    """
    rri_time = np.cumsum(rri) / 1000.0
    time_rri = rri_time - rri_time[0]
    time_rri_interp = np.arange(0, time_rri[-1], 1 / float(sf_up))
    tck = splrep(time_rri, rri, s=0)
    rri_interp = splev(time_rri_interp, tck, der=0)
    return rri_interp


sf_up = 4
rri_interp = interp_cubic_spline(rri, sf_up)
hr = 1000 * (60 / rri_interp)
print('Mean HR: %.2f bpm' % np.mean(hr))

# Detrend and normalize
edr = detrend(hr)
edr = (edr - edr.mean()) / edr.std()

# Find respiratory peaks
resp_peaks, _ = find_peaks(edr, height=0, distance=sf_up)

# Convert to seconds
resp_peaks = resp_peaks
resp_peaks_diff = np.diff(resp_peaks) / sf_up

# Plot the EDR waveform
plt.plot(edr, '-')
plt.plot(resp_peaks, edr[resp_peaks], 'o')
_ = plt.title('ECG derived respiration')
