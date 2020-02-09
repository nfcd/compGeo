import math

def ZeroTwoPi(a):
    '''
    This function makes sure input azimuth (a)
    is within 0  to 2*pi (b)
    
    NOTE: Azimuths a and b are input/output in radians
    
    Python function translated from the Matlab function 
    ZeroTwoPi in Allmendinger et al. (2012)
    '''
    b=a
    twopi = 2*math.pi
    if b < 0:
        b += twopi
    elif b >= twopi:
        b -= twopi
    
    return b