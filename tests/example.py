import math
from medipy import hamilton
from medipy import lipponen
from matplotlib import pyplot as plt
from scipy.misc import electrocardiogram

# Load ECG Data using scipy
ecg = electrocardiogram()
sample_rate = 360

# Detect R-Peaks
ecg_preprocessed = hamilton.preprocessing(ecg, sample_rate)
r_peaks = hamilton.detect_r_peaks(ecg, ecg_preprocessed, sample_rate, least_distance=0.2, th_coefficient=0.189, th_search_back=0.3)
r_peaks_refined = hamilton.refine_r_peaks(ecg, r_peaks, refine_samples=30)

# Calculate RR-Intervals
rr_intervals = [math.ceil(r_peaks[k] - r_peaks[k - 1]) for k in range(1, len(r_peaks))]

# Detect RR-Errors
rr_errors, ectopic_intervals, missed_intervals, extra_intervals, long_short_intervals = lipponen.rr_interval_error_detection(rr_intervals)

# Print Results
print("R-Peaks: ", len(r_peaks))
print("RR-Intervals: ", len(rr_intervals))
print("RR-Errors: ", sum(rr_errors))
print("Ectopic Intervals: ", sum(ectopic_intervals))
print("Missed Intervals: ", sum(missed_intervals))
print("Extra Intervals: ", sum(extra_intervals))
print("Long Short Intervals: ", sum(long_short_intervals))

# Plot
fig = plt.figure()
plt.plot(ecg)
plt.plot(r_peaks, ecg[r_peaks], 'rx')
plt.plot(r_peaks_refined, ecg[r_peaks_refined], 'go')
plt.show()
