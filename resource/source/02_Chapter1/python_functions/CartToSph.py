import math
from ZeroTwoPi import ZeroTwoPi as ZeroTwoPi

def CartToSph(cn,ce,cd):
    '''
    CartToSph converts from cartesian to spherical coordinates 

    CartToSph(cn,ce,cd) returns the trend (trd)
    and plunge (plg) of a line for input north (cn), east (ce), 
    and down (cd) direction cosines

    NOTE: Trend and plunge are returned in radians

    CartToSph uses function ZeroTwoPi
    
    Python function translated from the Matlab function CartToSph in the book:
    Structural Geology Algorithms by Allmendinger, Cardozo, & Fisher, 2012
    '''
    # Plunge 
    plg = math.asin(cd)
    
    #Trend
    #If north direction cosine is zero, trend is east or west
    #Choose which one by the sign of the east direction cosine
    if cn == 0.0:
        if ce < 0.0:
            trd = 3.0/2.0*math.pi # trend is west
        else:
            trd = math.pi/2.0 # trend is east
    # Else
    else:
        trd = math.atan(ce/cn)
        if cn < 0.0:
            #Add pi 
            trd = trd+math.pi
        # Make sure trd is between 0 and 2*pi
        trd = ZeroTwoPi(trd)
    return trd, plg
