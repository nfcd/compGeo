from uncertainties import umath as umath
import numpy as np
import uncertainties as unc

def TrueThicknessU(strike,dip,top,base):
    '''
    TrueThicknessU calculates the thickness (t) of a unit 
    given the strike (strike) and dip (dip) of the unit, 
    and points at its top (top) and base (base).
    strike and dip, as well as the points have
    uncertainties.

    top and base are 1 x 3 arrays defining the location 
    of top and base points in an ENU coordinate system. 
    For each one of these arrays, the first, second
    and third entries are the E, N and U coordinates. 
    These coordinates have uncertainties

    NOTE: strike and dip should be input in radians and
            they have uncertainties in radians. The
            returned thickness have also uncertainties
    '''
    # make the transformation matrix from ENU coordinates 
    # to SDP coordinates. Eq. 5.10
    sinStr = umath.sin(strike)
    cosStr = umath.cos(strike)
    sinDip = umath.sin(dip)
    cosDip = umath.cos(dip)
    a = np.array([[sinStr, cosStr, unc.ufloat(0,0)],
    [cosStr*cosDip, -sinStr*cosDip, -sinDip],
    [-cosStr*sinDip, sinStr*sinDip, -cosDip]])
    
    # transform the top and base points 
    # from ENU to SDP coordinates. Eq. 5.4 
    topn = np.array([unc.ufloat(0,0), unc.ufloat(0,0), unc.ufloat(0,0)])
    basen = np.array([unc.ufloat(0,0), unc.ufloat(0,0), unc.ufloat(0,0)])
    for i in range(0,3,1):
        for j in range(0,3,1):
            topn[i] = a[i,j]*top[j] + topn[i]
            basen[i] = a[i,j]*base[j] + basen[i]
    
    # compute the thickness of the unit. Eq. 5.12
    t = np.abs(basen[2] - topn[2])
    
    return t