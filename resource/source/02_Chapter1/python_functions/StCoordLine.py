# -*- coding: utf-8 -*-
import math

def StCoordLine(trd, plg, sttype):
    '''
    StCoordLine computes the coordinates of a line 
    in an equal angle or equal area stereonet of unit radius
    
    trd         = trend of line
    plg         = plunge of line
    sttype      = An integer indicating the type of stereonet. 0 for equal angle
                  and 1 for equal area
    xp and yp   = Coordinates of the line in the stereonet plot

    NOTE: trend and plunge should be entered in radians

    StCoordLine uses function ZeroTwoPi
    Python function translated from the Matlab function StCoordLine
    of Structural Geology Algorithms by Allmendinger, Cardozo, & Fisher, 2011
    '''
    
    # Take care of negative plunges
    if plg < 0:
        trd = ZeroTwoPi(trd+math.pi)
        plg = -plg
        
    # Some constants
    piS4 = math.pi/4
    s2 = math.sqrt(2)
    plgS2 = plg/2
        
    # Equal angle stereonet: From Equation 1.5 above
    # Also see Pollard and Fletcher (2005), eq.2.72
    if sttype == 0:
        xp = math.tan(piS4 - plgS2)*math.sin(trd)
        yp = math.tan(piS4 - plgS2)*math.cos(trd)
    # Equal area stereonet: From Equation 1.6 above
    # Also see Pollard and Fletcher (2005), eq.2.90
    elif sttype == 1:
        xp = s2*math.sin(piS4 - plgS2)*math.sin(trd)
        yp = s2*math.sin(piS4 - plgS2)*math.cos(trd)
    
    return xp, yp