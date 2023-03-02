'''
lipponen.py
'''

import math
import numpy as np


def rr_interval_error_detection(rr_intervals_list):
    '''
    This method detects errors of a given rr interval list with an modified version of Lipponen and Tarvainen  algorithm
    '''
    # Preparation

    rr_intervals = np.array(rr_intervals_list)
    rr_errors = np.zeros(len(rr_intervals))
    alpha = 5.2
    drr_intervals = np.diff(rr_intervals)
    drr_intervals = np.insert(drr_intervals, 0, 0)

    medrr_intervals = []
    mrr_intervals = []
    for index, rr_interval in enumerate(rr_intervals):
        if index < 10 or index > (len(rr_intervals) - 10):
            local_medrr = np.median(rr_intervals[0:20])
            medrr_intervals.append(local_medrr)
            mrr_intervals.append(rr_interval - local_medrr)
        else:
            local_medrr = np.median(rr_intervals[index - 10:index + 11])
            medrr_intervals.append(local_medrr)
            mrr_intervals.append(rr_interval - local_medrr)

    medrr_intervals = np.array(medrr_intervals)
    mrr_intervals = np.array(mrr_intervals)

    threshold_1 = []
    threshold_2 = []
    for index, rr_interval in enumerate(rr_intervals):
        if index < 45 or index > (len(rr_intervals) - 45):
            local_drr_intervals = drr_intervals[0:91]
            threshold_1.append(((np.quantile(abs(local_drr_intervals), 0.75) - np.quantile(abs(local_drr_intervals), 0.25)) / 2) * alpha)
            local_mrr_intervals = np.array(mrr_intervals[0:91])
            threshold_2.append(((np.quantile(abs(local_mrr_intervals), 0.75) - np.quantile(abs(local_mrr_intervals), 0.25)) / 2) * alpha)
        else:
            local_drr_intervals = drr_intervals[index - 45:index + 46]
            threshold_1.append(((np.quantile(abs(local_drr_intervals), 0.75) - np.quantile(abs(local_drr_intervals), 0.25)) / 2) * alpha)
            local_mrr_intervals = np.array(mrr_intervals[index - 45:index + 46])
            threshold_2.append(((np.quantile(abs(local_mrr_intervals), 0.75) - np.quantile(abs(local_mrr_intervals), 0.25)) / 2) * alpha)

    threshold_1 = np.array(threshold_1)
    threshold_2 = np.array(threshold_2)

    # Normalization
    drr_intervals_n = np.divide(drr_intervals, threshold_1)
    mrr_intervals_n = np.divide(mrr_intervals, threshold_2)

    # Constants
    const_1 = 0.13
    const_2 = 0.17
    index = 1
    ectopic_intervals = np.zeros(len(rr_intervals))
    missed_intervals = np.zeros(len(rr_intervals))
    extra_intervals = np.zeros(len(rr_intervals))
    long_short_intervals = np.zeros(len(rr_intervals))

    # Detektionsdurchlauf
    while index < len(drr_intervals_n) - 2:
        if abs(drr_intervals_n[index]) > 1:
            s11 = drr_intervals_n[index]
            if s11 > 0:
                s12 = max(drr_intervals_n[index - 1], drr_intervals_n[index + 1])
            else:
                s12 = min(drr_intervals_n[index - 1], drr_intervals_n[index + 1])
            eq1 = ((s11 > 1) and (s12 < -const_1 * s11 - const_2))
            eq2 = ((s11 < -1) and (s12 > -const_1 * s11 + const_2))
            if eq1:  # Ectopic NPN
                ectopic_intervals[index] = 1
                rr_errors[index] = 1
                if drr_intervals_n[index - 1] > drr_intervals_n[index + 1]:
                    ectopic_intervals[index - 1] = 1
                    rr_errors[index - 1] = 1
                    index += 1
                else:
                    ectopic_intervals[index + 1] = 1
                    rr_errors[index + 1] = 1
                    index += 2
                continue
            elif eq2:  # Ectopic PNP
                ectopic_intervals[index] = 1
                rr_errors[index] = 1
                if drr_intervals_n[index - 1] < drr_intervals_n[index + 1]:
                    ectopic_intervals[index - 1] = 1
                    rr_errors[index - 1] = 1
                    index += 1
                else:
                    ectopic_intervals[index + 1] = 1
                    rr_errors[index + 1] = 1
                    index += 2
                continue
            elif abs(mrr_intervals_n[index]) > 1:
                s21 = drr_intervals_n[index]
                if s21 >= 0:
                    s22 = min(drr_intervals_n[index + 1], drr_intervals_n[index + 2])
                else:
                    s22 = max(drr_intervals_n[index + 1], drr_intervals_n[index + 2])
                eq3 = (s21 < -1 and s22 > 1)
                eq4a = ((mrr_intervals_n[index]) < -1)
                eq4b = ((mrr_intervals_n[index]) > 1)
                eq5 = (s21 > 1 and s22 < -1)
                if eq3 and eq4a:
                    eq6 = abs(rr_intervals[index] + rr_intervals[index + 1] - medrr_intervals[index]) * 0.8 < threshold_2[index]
                    if eq6:  # Extra
                        extra_intervals[index:index + 2] = 1
                        rr_errors[index:index + 2] = 1
                        index += 2
                        continue
                    else:  # Short
                        long_short_intervals[index] = 1
                        rr_errors[index] = 1
                        index += 1
                        continue
                elif eq5 and eq4b:
                    eq7 = rr_intervals[index] / medrr_intervals[index] > 2
                    if eq7:  # Missed
                        weight = math.floor(np.divide(rr_intervals[index], medrr_intervals[index]))
                        missed_intervals[index] = int(weight)  # missed intervals have bigger impact on the error
                        rr_errors[index] = int(weight)
                        index += 1
                        continue
                    else:  # Long
                        long_short_intervals[index] = 1
                        rr_errors[index] = 1
                        index += 1
                        continue
                else:
                    index += 1
            else:
                index += 1
        else:
            index += 1

    return rr_errors, ectopic_intervals, missed_intervals, extra_intervals, long_short_intervals
