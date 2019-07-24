import scipy.signal as signal
from scipy.signal import butter, lfilter
import numpy as np
import math
import peakutils as pu
from time import gmtime, strftime



class rdectoren():

    def __init__(self, ecgsignal, samplingfrequency):
        '''Init der BaseClass Objektvariablen'''
        self.signal = ecgsignal
        self.fs = samplingfrequency
        self.dt = (1 / self.fs) * 1000
        self.rpeaks = []
        self.rrintervall = []


class hamilton(rdectoren):

    def __init__(self, ecgsignal, samplingfrequency):
        super().__init__(ecgsignal, samplingfrequency)

    def hamiltonalgorithmus(self):

        f1 = 8 / self.fs
        f2 = 16 / self.fs

        b, a = signal.butter(1, [f1 * 2, f2 * 2], btype='bandpass')

        filtered_ecg = signal.lfilter(b, a, self.signal)

        diff = abs(np.diff(filtered_ecg))

        b = np.ones(int(0.08 * self.fs))
        b = b / int(0.08 * self.fs)
        a = [1]

        ma = signal.lfilter(b, a, diff)

        ma[0:len(b) * 2] = 0

        peaks, _ = signal.find_peaks(ma, distance=(0.25 * self.fs))

        n_pks = []
        n_pks_ave = 0.0
        s_pks = []
        s_pks_ave = 0.0
        QRS = []
        RR = []
        RR_ave = 0.0

        th = 0.0

        i = 0
        idx = []
        for peak in peaks:

            if ma[peak] > th:
                QRS.append(peak)
                idx.append(i)
                s_pks.append(ma[peak])
                if len(n_pks) > 8:
                    s_pks.pop(0)
                s_pks_ave = np.mean(s_pks)

                if RR_ave != 0.0:
                    if QRS[-1] - QRS[-2] > 1.5 * RR_ave:
                        missed_peaks = peaks[idx[-2] + 1:idx[-1]]
                        for missed_peak in missed_peaks:
                            if missed_peak - peaks[idx[-2]] > int(0.360 * self.fs) and ma[missed_peak] > 0.5 * th:
                                QRS.append(missed_peak)
                                QRS.sort()
                                break

                if len(QRS) > 2:
                    RR.append(QRS[-1] - QRS[-2])
                    if len(RR) > 8:
                        RR.pop(0)
                    RR_ave = int(np.mean(RR))

            else:
                n_pks.append(ma[peak])
                if len(n_pks) > 8:
                    n_pks.pop(0)
                n_pks_ave = np.mean(n_pks)

            th = n_pks_ave + 0.45 * (s_pks_ave - n_pks_ave)

            i += 1

        QRS.pop(0)

        window = int(0.1 * self.fs)

        rpeaks = hamilton.searchBack(self, QRS, window)

        return rpeaks

    def searchBack(self, detected_peaks, search_samples):

        r_peaks = []
        window = search_samples

        for i in detected_peaks:
            if i < window:
                section = self.signal[:i]
                r_peaks.append(np.argmax(section))
            else:
                section = self.signal[i - window:i]
                r_peaks.append(np.argmax(section) + i - window)
        return np.array(r_peaks)

    def detect(self):
        self.rpeaks = hamilton.hamiltonalgorithmus(self)
        for k in range(1, len(self.rpeaks)):
            self.rrintervall.append(
                math.ceil((self.rpeaks[k] - self.rpeaks[k - 1]) * self.dt))
        return self.rpeaks, self.rrintervall


class peakutils(rdectoren):
    def __init__(self, ecgsignal, samplingfrequency):
        super().__init__(ecgsignal, samplingfrequency)

    def peakutilsalgorithmus(self):
        indexes = pu.indexes(self.signal, thres=0.6, min_dist=350)
        return indexes

    def detect(self):
        self.rpeaks = peakutils.peakutilsalgorithmus(self)
        for k in range(1, len(self.rpeaks)):
            self.rrintervall.append(
                math.ceil((self.rpeaks[k] - self.rpeaks[k - 1]) * self.dt))
        return self.rpeaks, self.rrintervall


class skipi(rdectoren):
    def __init__(self, ecgsignal, samplingfrequency):
        super().__init__(ecgsignal, samplingfrequency)

    def skipialgorithmus(self):
        return signal.find_peaks(self.signal, height=0.7, distance=100)[0]

    def detect(self):
        self.rpeaks = skipi.skipialgorithmus(self)
        for k in range(1, len(self.rpeaks)):
            self.rrintervall.append(
                math.ceil((self.rpeaks[k] - self.rpeaks[k - 1]) * self.dt))
        return self.rpeaks, self.rrintervall


class ownr(rdectoren):
    def __init__(self, ecgsignal, samplingfrequency):
        super().__init__(ecgsignal, samplingfrequency)

    def ownralgorithmus(self):
        '''
        Selbst Entwickelter Algorithmus zur Detektion von R-Zacken
        '''
        # WEITERE PROBLEME LAST UND FIRST PEAK VON BENACHBARTEN REGISTERN MÜSSEN GLOBAL ÜBEREINSTIMMEN
        # Zeitliche Schwellwerte
        # Preporecessing HP TP..
        # Settingsparameter einrichten

        # Basierend auf aktuellen Werten
        used = 0
        stop = 0
        indexes = []

        while not stop:
            register, used, stop, currentindex = ownr.getnewvalues(self, used)
            threshold_up, threshold_down = ownr.calculatethresholds(
                self, register)
            indexes.extend(ownr.detectrzacken(self, register, threshold_up,
                                              threshold_down, currentindex))
            ownr.checker(self, indexes)

        return np.asarray(indexes)

    def getnewvalues(self, used):
        numberofvals = 30 * self.fs
        k = used * numberofvals  # FEHLER Updaste of used, makes k big
        if (k + numberofvals) >= len(self.signal):  # FEHHLER
            register = self.signal[k:-1]
            stop = 1
        else:
            register = self.signal[k:k + numberofvals]
            stop = 0
        used += 1
        return register, used, stop, k

    def calculatethresholds(self, register):
        # Werte aus akutellem Register berechnen
        meanreg = np.asarray(register).mean()
        maxreg = np.asarray(register).max()
        threshold_up = 0.6 * (maxreg - meanreg)
        threshold_down = 0.3 * (maxreg - meanreg)
        return threshold_up, threshold_down
        # Timethresholds

    def detectrzacken(self, register, threshold_up, threshold_down, currentindex):
        treshval = []
        maxofaround = []
        counter = currentindex
        peaks = []

        # For Schleife in Register
        for val in register:
            if val > threshold_up:
                treshval.append([counter, 1, val])
            else:
                treshval.append([counter, 0, val])
            counter += 1

        # For Schleife mit markierten Werte die höher als TreshUp liegen
        for i in treshval:
            if i[1] == 1:
                maxofaround.append(i)

                # Falls Ende erreicht
                if i[0] == treshval[-1][0]:
                    if maxofaround:  # Und wenn überhaupt Werte größere Schwellwert vorhanden sind
                        currentmax = maxofaround[0]
                        for t in maxofaround:
                            if t[2] > currentmax[2]:
                                currentmax = t
                        peaks.append(currentmax[0])
                        # allpeaks.append(currentmax[0])
                        # Liste Leeren
                        maxofaround = []

            else:  # Wenn keine Werte mehr über Schwellwert
                if maxofaround:  # Und wenn überhaupt Werte größer Schwellwert vorhanden sind
                    currentmax = maxofaround[0]
                    for t in maxofaround:
                        if t[2] > currentmax[2]:
                            currentmax = t
                    peaks.append(currentmax[0])
                    # allpeaks.append(currentmax[0])
                    # Liste Leeren
                    maxofaround = []
        return peaks  # globaleSamplewerte

    def checker(self, indexliste):
        '''
        Checkt ob es bei einem Registerübergang zu doppleter Detektion gekommen ist
        Falls JA: Aussortieren des Detektionswertes der einen geringeren Signalwert hat
        '''

        regtransit = math.ceil(self.fs * 0.4)  # 100ms Übergangsbereich
        if len(indexliste) < 2:
            return
        if indexliste[-1] - indexliste[-2] <= regtransit:
            if self.signal[indexliste[-1]] > self.signal[indexliste[-2]]:
                del indexliste[-2]
            else:
                del indexliste[-1]

    def detect(self):
        self.rpeaks = ownr.ownralgorithmus(self)
        for k in range(1, len(self.rpeaks)):
            self.rrintervall.append(
                math.ceil((self.rpeaks[k] - self.rpeaks[k - 1]) * self.dt))
        return self.rpeaks, self.rrintervall


class pantompkins(rdectoren):

    def __init__(self, ecgsignal, samplingfrequency):
        super().__init__(ecgsignal, samplingfrequency)
        # Configuration parameters.
        self.ecg_data_raw = self.signal
        self.signal_frequency = self.fs

        self.filter_lowcut = 0.0001
        self.filter_highcut = 15.0
        self.filter_order = 1

        # Change proportionally when adjusting frequency (in samples). !!!
        self.integration_window = 15

        self.findpeaks_limit = 0.35 #0.35
        # Change proportionally when adjusting frequency (in samples).!!
        self.findpeaks_spacing = 50 #50

        # Change proportionally when adjusting frequency (in samples).!!
        self.refractory_period = 120
        self.qrs_peak_filtering_factor = 0.125
        self.noise_peak_filtering_factor = 0.125
        self.qrs_noise_diff_weight = 0.25

        # Measured and calculated values.
        self.filtered_ecg_measurements = None
        self.differentiated_ecg_measurements = None
        self.squared_ecg_measurements = None
        self.integrated_ecg_measurements = None
        self.detected_peaks_indices = None
        self.detected_peaks_values = None

        self.qrs_peak_value = 0.0
        self.noise_peak_value = 0.0
        self.threshold_value = 0.0

        # Detection results.
        self.qrs_peaks_indices = np.array([], dtype=int)
        self.noise_peaks_indices = np.array([], dtype=int)

        # Final ECG data and QRS detection results array - samples with detected QRS are marked with 1 value.
        self.ecg_data_detected = None

    def detect_peaks(self):
        """
        Method responsible for extracting peaks from loaded ECG measurements data through measurements processing.
        """
        # Extract measurements from loaded ECG data.
        #ecg_measurements = self.ecg_data_raw[:, 1]

        # Measurements filtering - 0-15 Hz band pass filter.
        self.filtered_ecg_measurements = self.bandpass_filter(self.signal, lowcut=self.filter_lowcut,
                                                            highcut=self.filter_highcut, signal_freq=self.signal_frequency,
                                                            filter_order=self.filter_order)
        self.filtered_ecg_measurements[:5] = self.filtered_ecg_measurements[5]

        # Derivative - provides QRS slope information.
        self.differentiated_ecg_measurements = np.ediff1d(self.filtered_ecg_measurements)

        # Squaring - intensifies values received in derivative.
        self.squared_ecg_measurements = self.differentiated_ecg_measurements ** 2

        # Moving-window integration.
        self.integrated_ecg_measurements = np.convolve(self.squared_ecg_measurements, np.ones(self.integration_window))

        # Fiducial mark - peak detection on integrated measurements.
        self.detected_peaks_indices = self.findpeaks(data=self.integrated_ecg_measurements,
                                                    limit=self.findpeaks_limit,
                                                    spacing=self.findpeaks_spacing)

        self.detected_peaks_values = self.integrated_ecg_measurements[self.detected_peaks_indices]

    def detect_qrs(self):
        """
        Method responsible for classifying detected ECG measurements peaks either as noise or as QRS complex (heart beat).
        """
        for detected_peak_index, detected_peaks_value in zip(self.detected_peaks_indices, self.detected_peaks_values):

            try:
                last_qrs_index = self.qrs_peaks_indices[-1]
            except IndexError:
                last_qrs_index = 0

            # After a valid QRS complex detection, there is a 200 ms refractory period before next one can be detected.
            if detected_peak_index - last_qrs_index > self.refractory_period or not self.qrs_peaks_indices.size:
                # Peak must be classified either as a noise peak or a QRS peak.
                # To be classified as a QRS peak it must exceed dynamically set threshold value.
                if detected_peaks_value > self.threshold_value:
                    self.qrs_peaks_indices = np.append(self.qrs_peaks_indices, detected_peak_index)

                    # Adjust QRS peak value used later for setting QRS-noise threshold.
                    self.qrs_peak_value = self.qrs_peak_filtering_factor * detected_peaks_value + \
                                        (1 - self.qrs_peak_filtering_factor) * self.qrs_peak_value
                else:
                    self.noise_peaks_indices = np.append(self.noise_peaks_indices, detected_peak_index)

                    # Adjust noise peak value used later for setting QRS-noise threshold.
                    self.noise_peak_value = self.noise_peak_filtering_factor * detected_peaks_value + \
                                            (1 - self.noise_peak_filtering_factor) * self.noise_peak_value

                # Adjust QRS-noise threshold value based on previously detected QRS or noise peaks value.
                self.threshold_value = self.noise_peak_value + \
                                    self.qrs_noise_diff_weight * (self.qrs_peak_value - self.noise_peak_value)

        # Create array containing both input ECG measurements data and QRS detection indication column.
        # We mark QRS detection with '1' flag in 'qrs_detected' log column ('0' otherwise).
        measurement_qrs_detection_flag = np.zeros([len(self.ecg_data_raw[:, 1]), 1])
        measurement_qrs_detection_flag[self.qrs_peaks_indices] = 1
        self.ecg_data_detected = np.append(self.ecg_data_raw, measurement_qrs_detection_flag, 1)

    def bandpass_filter(self, data, lowcut, highcut, signal_freq, filter_order):
        """
        Method responsible for creating and applying Butterworth filter.
        :param deque data: raw data
        :param float lowcut: filter lowcut frequency value
        :param float highcut: filter highcut frequency value
        :param int signal_freq: signal frequency in samples per second (Hz)
        :param int filter_order: filter order
        :return array: filtered data
        """
        nyquist_freq = 0.5 * signal_freq
        low = lowcut / nyquist_freq
        high = highcut / nyquist_freq
        b, a = butter(filter_order, [low, high], btype="band")
        y = lfilter(b, a, data)
        return y

    def findpeaks(self, data, spacing=1, limit=None):
        """
        Janko Slavic peak detection algorithm and implementation.
        https://github.com/jankoslavic/py-tools/tree/master/findpeaks
        Finds peaks in `data` which are of `spacing` width and >=`limit`.
        :param ndarray data: data
        :param float spacing: minimum spacing to the next peak (should be 1 or more)
        :param float limit: peaks should have value greater or equal
        :return array: detected peaks indexes array
        """
        lene = data.size
        x = np.zeros(lene + 2 * spacing)
        x[:spacing] = data[0] - 1.e-6
        x[-spacing:] = data[-1] - 1.e-6
        x[spacing:spacing + lene] = data
        peak_candidate = np.zeros(lene)
        peak_candidate[:] = True
        for s in range(spacing):
            start = spacing - s - 1
            h_b = x[start: start + lene]  # before
            start = spacing
            h_c = x[start: start + lene]  # central
            start = spacing + s + 1
            h_a = x[start: start + lene]  # after
            peak_candidate = np.logical_and(peak_candidate, np.logical_and(h_c > h_b, h_c > h_a))

        ind = np.argwhere(peak_candidate)
        ind = ind.reshape(ind.size)
        if limit is not None:
            ind = ind[data[ind] > limit]
        return ind

    def pantompkinsalgorithmus(self):
        pantompkins.detect_peaks(self)
        #pantompkins.detect_qrs()
        # aus der geflaggten information eine rpeaks liste machen 
        return self.detected_peaks_indices

    def detect(self):
        self.rpeaks = pantompkins.pantompkinsalgorithmus(self)
        for k in range(1, len(self.rpeaks)):
            self.rrintervall.append(math.ceil((self.rpeaks[k] - self.rpeaks[k - 1]) * self.dt))
        return self.rpeaks, self.rrintervall
        


    
