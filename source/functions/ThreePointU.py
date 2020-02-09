import numpy as np
from uncertainties import umath as umath
from CartToSphU import CartToSphU as CartToSphU
from Pole import Pole as Pole

def ThreePointU(p1,p2,p3):
    '''
    ThreePointU calculates the strike (strike) and dip (dip) 
    of a plane given the east (E), north (N), and up (U) 
    coordinates of three non-collinear points on the plane
    
    p1, p2 and p3 are 1 x 3 vectors defining the location 
    of the points in an ENU coordinate system. For each one 
    of these vectors the first entry is the E coordinate, 
    the second entry the N coordinate, and the third entry 
    the U coordinate. These coordinates have uncertainties
    
    NOTE: strike and dip are returned in radians and they 
    follow the right-hand rule format. The returned 
    uncertainties are also in radians

    ThreePointU uses functions CartToSphU and Pole
    It also uses the uncertainties package from
    Eric O. Lebigot
    '''
    # if points are given as lists, 
    # they must be converted to np.arrays
    p1 = np.array(p1)
    p2 = np.array(p2)
    p3 = np.array(p3)
    
    # make vectors v (p1 - p3) and u (p2 - p3)
    v = p1 - p2
    u = p2 - p3
    
    # take the cross product of v and u
    vcu0 = v[1]*u[2] - v[2]*u[1] 
    vcu1 = v[2]*u[0] - v[0]*u[2]
    vcu2 = v[0]*u[1] - v[1]*u[0]
    
    # magnitude of the vector
    r = umath.sqrt(vcu0*vcu0+vcu1*vcu1+vcu2*vcu2)
    
    # unit vector
    uvcu0 = vcu0/r
    uvcu1 = vcu1/r
    uvcu2 = vcu2/r
    
    # make the pole vector in NED coordinates
    pole0 = uvcu1
    pole1 = uvcu0
    pole2 = -uvcu2
    
    # Make sure pole point downwards
    if pole2 < 0:
        pole0 *= -1
        pole1 *= -1
        pole2 *= -1
        
    # find the trend and plunge of the pole
    trd, plg = CartToSphU(pole0,pole1,pole2)
    
    # find strike and dip of plane
    strike, dip = Pole(trd, plg, 0)
    
    return strike , dip