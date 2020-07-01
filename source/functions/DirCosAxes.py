import numpy as np

from SphToCart import SphToCart as SphToCart

def DirCosAxes(tX1,pX1,tX3):
    '''
    DirCosAxes calculates the direction cosines of a right 
    handed, orthogonal X1X2X3 cartesian coordinate system 
    of any orientation with respect to NED

    USE: dC = DirCosAxes(tX1,pX1,tX3)

    tX1 = trend of X1
    pX1 = plunge of X1
    tX3 = trend of X3
    dC = 3 x 3 matrix containing the direction cosines 
         of X1 (row 1), X2 (row 2), and X3 (row 3)

    Note: Input angles should be in radians

    DirCosAxes uses function SphToCart

    Python function translated from the Matlab function
    DirCosAxes in Allmendinger et al. (2012)
    '''
    
    # Some constants
    east = np.pi/2.0
    west = 3.0*east
    # tolerance for near zero values
    tol = 1e-6 
    
    # Initialize matrix of direction cosines
    dC = np.zeros((3,3))
    
    # Direction cosines of X1
    dC[0,0],dC[0,1],dC[0,2] = SphToCart(tX1,pX1,0)
    
    # Calculate plunge of axis 3
    # If axis 1 is horizontal
    if abs(pX1) < tol:
        dt = abs(tX1-tX3)
        if abs(dt - east) < tol and abs(dt - west) < tol:  # does and really work in the intended way?
            pX3 = 0.0
        else:
            pX3 = east
    # Else
    else:
        # From dot product, theta = 90, cos(theta) = 0
        pX3 = np.arctan(-(dC[0,0]*np.cos(tX3)+dC[0,1]*np.sin(tX3))/dC[0,2])

    # Direction cosines of X3
    dC[2,0],dC[2,1],dC[2,2] = SphToCart(tX3,pX3,0)
    
    # X2 is the cross product of X3 and X1
    dC[1,:] = np.cross(dC[2,:],dC[0,:])
    # Make it a unit vector
    dC[1,:] = dC[1,:]/np.linalg.norm(dC[1,:])
    
    return dC