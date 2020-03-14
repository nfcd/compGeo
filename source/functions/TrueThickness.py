import numpy as np

def TrueThickness(strike,dip,top,base):
    '''
    TrueThickness calculates the thickness (t) of a unit 
    given the strike (strike) and dip (dip) of the unit, 
    and points at its top (top) and base (base) 

    top and base are 1 x 3 arrays defining the location 
    of top and base points in an ENU coordinate system. 
    For each of these arrays, the first entry is the E 
    coordinate, the second entry the N coordinate, 
    and the third entry the U coordinate

    NOTE: strike and dip should be input in radians
    '''
    # make the transformation matrix from ENU coordinates 
    # to SDP coordinates. Eq. 5.10
    sinStr = np.sin(strike)
    cosStr = np.cos(strike)
    sinDip = np.sin(dip)
    cosDip = np.cos(dip)
    a = np.array([[sinStr, cosStr, 0],
    [cosStr*cosDip, -sinStr*cosDip, -sinDip],
    [-cosStr*sinDip, sinStr*sinDip, -cosDip]])
    
    # transform the top and base points 
    # from ENU to SDP coordinates. Eq. 5.4 
    topn = np.zeros(3)
    basen = np.zeros(3)
    for i in range(0,3,1):
        for j in range(0,3,1):
            topn[i] = a[i,j]*top[j] + topn[i]
            basen[i] = a[i,j]*base[j] + basen[i]
    
    # compute the thickness of the unit. Eq. 5.12
    t = np.abs(basen[2] - topn[2])
    
    return t