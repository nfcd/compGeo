import numpy as np

from Pole import Pole as Pole
from CartToSph import CartToSph as CartToSph

def FitPlane(pts):
    '''
    Fitplane computes the best fit plane for a group of points (position
    vectors) on the plane

    USE: strike, dip = FitPlane(pts)

    pts is a n x 3 matrix containing the East (column 1), North (column 2),
    and Up (column 3) coordinates of n points on the plane

    strike and dip are returned in radians

    FitPlane uses functions Pole and CartToSph 
    
    Python function translated from the Matlab function 
    FitPlane in Allmendinger et al. (2012)
    '''
    
    # Compute the centroid of the selected points
    avge = np.mean(pts[:,0])
    avgn = np.mean(pts[:,1])
    avgu = np.mean(pts[:,2])
    
    # Compute the points vectors minus the centroid
    pts[:,0] = pts[:,0] - avge
    pts[:,1] = pts[:,1] - avgn
    pts[:,2] = pts[:,2] - avgu
    
    # Compute the covariance/orientation matrix
    a = np.zeros((3,3))
    for i in range(pts.shape[0]):
        ce = pts[i,0]
        cn = pts[i,1]
        cu = pts[i,2]
        # compute orientation matrix 
        a[0,0] = a[0,0] + ce*ce
        a[0,1] = a[0,1] + ce*cn
        a[0,2] = a[0,2] + ce*cu
        a[1,1] = a[1,1] + cn*cn
        a[1,2] = a[1,2] + cn*cu
        a[2,2] = a[2,2] + cu*cu
    # The orientation matrix is symmetric so the off-diagonal components can
    # be equated
    a[1,0] = a[0,1]
    a[2,0] = a[0,2]
    a[2,1] = a[1,2]
    
    # calculate the eigenvalues and eigenvectors of the orientation matrix
    # use Matlab function eig
    _,V = np.linalg.eigh(a)
    
    # Calculate pole to best-fit plane = lowest eigenvalue vector
    # in E, N, D coordinates
    ce = V[0,0]
    cn = V[1,0]
    cd = -V[2,0]
    
    # Find trend and plunge of pole to best fit plane
    trd, plg = CartToSph(cn,ce,cd)
    
    # Find Best fit plane
    strike, dip = Pole(trd,plg,0)
    
    return strike, dip