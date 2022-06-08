'''
varon.py
'''
import numpy as np


def hear_rate_decomposition(hrv, edr, frequency):
    '''
    Using ARMAX/OSP
    m first try const for complete data set, then vars for different windows
    frequency is int
    # Assume 5Hz
    # Because frequency of the input and output is the same
    # For every new edr and hrv val, there is a new model to create
    # Average Model order is 1 to 3s

    # Create empty vector same lentgh
    # hrv_x
    # hrv_y
    '''
    # Check length of both signals
    if len(hrv) != len(edr):
        print('Signals have not the same length')

    # Limits
    maximum_order = 10 * frequency  # Aussume min Resp. Freq. 0.1 Hz - > 10s; Widjaja et al. Average Model 1s-3s

    # Calculate the model order m, as the minimum value between MDL and AIC.
    model_order = maximum_order

    # Create the matrix v_spann = [X0 , X1 , . . . , XM]
    v_spann = []
    for delays in range(model_order):
        x_vector = edr[delays:len(edr) - model_order + delays]
        v_spann.append(x_vector)
    v_spann = np.array(v_spann)

    # Compute the projection matrix P using P = V (V.T V)−1 V.T
    p0_matrix = np.linalg.inv(np.matmul(v_spann, v_spann.T))
    p1_matrix = np.matmul(v_spann, p0_matrix)
    p2_matrix = np.matmul(p1_matrix, v_spann.T)

    # Calculate the respiratory component Yx = PY
    hrv_x = np.matmul(p2_matrix, hrv)

    # Calculate the residual component Yy = Y − Yx
    hrv_y = hrv - hrv_x

    return hrv_y, hrv_y
