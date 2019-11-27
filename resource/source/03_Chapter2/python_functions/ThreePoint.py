import numpy as np
from numpy import linalg as la

# Note: str is changed to strike, because str has a meaning in python and can
# therefore not be used as variable name

def ThreePoint(p1,p2,p3):
    '''
    Three points calculate the strike (strike) and dip (dip) of a plane given the
    east (E), north (N), and up (U) coordinates of three non-collinear points
    on the plane
    
    p1, p2 and p3 are 1 x 3 vectors defining the location of the points in an
    ENU coordinate system. For each one of these vectors the first entry is
    the E coordinate, the second entry the N coordinate, and the third entry
    the U coordinate

    NOTE: strike and dip are returned in radians and following the right-hand
    rule convention

    ThreePoint uses functions CartToSph and ZeroTwoPi
    
    Python function translated from the Matlab function ThreePoints in the book:
    Structural Geology Algorithms by Allmendinger, Cardozo, & Fisher, 2012
    '''
    # if points are given as lists, they must be converted to np.arrays
    p1 = np.array(p1)
    p2 = np.array(p2)
    p3 = np.array(p3)
    # make vectors v (p1 - p3) and u (p2 - p3)
    v = p1 - p2
    u = p2 - p3
    
    # take the cross product of v and u
    vcu = np.cross(v,u)
    
    # make this vector a unit vector by dividing it by its magnitude
    mvcu = la.norm(vcu) # magnitude of the vector
    uvcu = vcu/mvcu # unit vector
    
    # make the pole vector in NED coordinates
    p = [uvcu[2], uvcu[1], -uvcu[3]]
    
    # if pole is pointing upwards make it point downwards
    if p[3] < 0:
        p = p * (-1)
        
    # find the trend and plunge of the pole
    # use function CartToSph
    trd, plg = CartToSph(p[1],p[2],p[3])
    
    # find strike and dip of plane
    # make sure strike is between 0 and 360
    strike = ZeroTwoPi(trd + np.pi/2.0); 
    dip = (np.pi/2.0 - plg); 
    
    return strike , dip