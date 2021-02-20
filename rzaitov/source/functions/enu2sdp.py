import numpy as np

def enu2sdp(strike, dip):
    '''
    NOTE: strike and dip should be input in radians
    '''
    # make the transformation matrix from ENU coordinates
    # to SDP coordinates. Eq. 5.10
    sin_str = np.sin(strike)
    cos_str = np.cos(strike)
    sin_dip = np.sin(dip)
    cos_dip = np.cos(dip)
    return np.array([
        [sin_str, cos_str, 0],
        [cos_str*cos_dip, -sin_str*cos_dip, -sin_dip],
        [-cos_str*sin_dip, sin_str*sin_dip, -cos_dip]
    ])