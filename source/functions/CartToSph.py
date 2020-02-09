import math
from ZeroTwoPi import ZeroTwoPi as ZeroTwoPi

def CartToSph(cn,ce,cd):
    '''
    CartToSph converts from Cartesian to spherical coordinates 

    CartToSph(cn,ce,cd) returns the trend (trd)
    and plunge (plg) of a line for input north (cn), 
    east (ce), and down (cd) direction cosines

    NOTE: Trend and plunge are returned in radians

    CartToSph uses function ZeroTwoPi
    
    Python function translated from the Matlab function 
    CartToSph in Allmendinger et al. (2012)
    '''
    pi = math.pi
    # Plunge 
    plg = math.asin(cd) # Eq. 4.13a
    
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
        trd = math.atan(ce/cn) # Eq. 4.14a
        if cn < 0.0:
            #Add pi 
            trd = trd+pi # Eq. 4.14b
        # Make sure trend is between 0 and 2*pi
        trd = ZeroTwoPi(trd)
    
    return trd, plg