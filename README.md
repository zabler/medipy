# medipy
Code Collection for Medical Data Processing and Analysis in Python

<!--![Version](https://img.shields.io/pypi/v/your-package-name.svg) -->
[![License](https://img.shields.io/badge/License-MIT-blue.svg)](LICENSE)

This Python package provides a broad collection of code for processing and analyzing medical data (like ECG, EEG, ACC, EMG, GYRO, ..). Every modul reflects a model, algorithm or concept published in scientific journal (or similiar) and is named after its main author. The publications are cited in the docstring of each modul.

## Implemented Code
These are the current modules which are implemented in this package:
- R-Peak-Detection for 1D ECG Signals by [Hamilton](src/medipy/hamilton.py)
- RR-Interval-Error-Detection for 1D ECG Signals by [Lipponen](src/medipy/lipponen.py)

## Installation
You can install the package using pip:

```bash
pip install medipy
```

## Usage
You can use any function of this package by either importing the complete package
```python
import medipy as mp
r_peaks= mp.hamilton.detect_r_peaks(ecg, sampling_rate)
```
or single modules

```python
from medipy import hamilton
r_peaks= hamilton.detect_r_peaks(ecg, sampling_rate)
```

## Example
```python
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
```

## Contribution 
If you want to support the project or contribute to the code, feel free to contact.

## License
This project is licensed under the MIT License - see the [License](LICENSE) file for details.



