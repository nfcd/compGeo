import math
from uncertainties import umath as umath
from ZeroTwoPi import ZeroTwoPi as ZeroTwoPi

def CartToSphU(cn,ce,cd):
    '''
    CartToSphU converts from Cartesian to spherical coordinates 

    CartToSphU(cn,ce,cd) returns the trend (trd)
    and plunge (plg) of a line for input north (cn), 
    east (ce), and down (cd) direction cosines.
    Notice that the input direction cosines have uncertainties

    NOTE: Trend and plunge are returned in radians and they
    		have uncertainties in radians

    CartToSphU uses function ZeroTwoPi
    It also uses the uncertainties package from
    Eric O. Lebigot
    
    Based on Python function CartToSph
    '''
    pi = math.pi
    # Plunge 
    plg = umath.asin(cd) # Eq. 4.13a
    
    #Trend
    #If north direction cosine is zero, trend is east or west
    #Choose which one by the sign of the east direction cosine
    if cn == 0.0:
        if ce < 0.0:
            trd = 3.0/2.0*pi # Eq. 4.14d, trend is west
        else:
            trd = pi/2.0 # Eq. 4.14c, trend is east
    # Else
    else:
        trd = umath.atan(ce/cn) # Eq. 4.14a
        if cn < 0.0:
            #Add pi 
            trd = trd+pi # Eq. 4.14b
        # Make sure trend is between 0 and 2*pi
        trd = ZeroTwoPi(trd)
    
    return trd, plg