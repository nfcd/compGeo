import math

def ZeroTwoPi(a):
    '''
    a = input azimuth
    b = output azimuth from 0 to 2*pi
    
    NOTE: Azimuths a and b are input/output in radians
    
    Python function translated from the Matlab function ZeroTwoPi in the book:
    Structural Geology Algorithms by Allmendinger, Cardozo, & Fisher, 2012
    '''
    b=a
    twopi = 2*math.pi
    if b < 0:
        b += twopi
    elif b >= twopi:
        b -= twopi
    return b
