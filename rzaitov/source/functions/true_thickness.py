import numpy as np

def true_thickness(T, top_enu, base_enu):
    '''
    T – coordinate transformation from ENU to SDP

    top and base are 1 x 3 arrays defining the location
    of top and base points in an ENU coordinate system.
    For each one of these arrays, the first, second
    and third entries are the E, N and U coordinates
    '''
    top_sdp = np.dot(T, top_enu)
    base_sdp = np.dot(T, base_enu)

    # compute the thickness of the unit. Eq. 5.12
    return np.abs(base_sdp[2] - top_sdp[2])