# -*- coding: utf-8 -*-
import math

def SphToCart(trd,plg,k):
    '''
    SphToCart converts from spherical to cartesian coordinates 

    SphToCart(trd,plg,k) returns the north (cn), 
    east (ce), and down (cd) direction cosines of a line.

    k is an integer to tell whether the trend and plunge of a line 
    (k = 0) or strike and dip of a plane in right hand rule 
    (k = 1) are being sent in the trd and plg slots. In this 
    last case, the direction cosines of the pole to the plane 
    are returned
 
    NOTE: Angles should be entered in radians
    
    Python function translated from the Matlab function SphTpCart
    of Structural Geology Algorithms by Allmendinger, Cardozo, & Fisher, 2011
    '''
    # If line (see Table 2.1)
    if k == 0:
        cd = math.sin(plg)
        ce = math.cos(plg) * math.sin(trd)
        cn = math.cos(plg) * math.cos(trd)
    # Else pole to plane (see Table 2.1)
    elif k == 1:
        cd = math.cos(plg)
        ce = -math.sin(plg) * math.cos(trd)
        cn = math.sin(plg) * math.sin(trd)
    return cn, ce, cd