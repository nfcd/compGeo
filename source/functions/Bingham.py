import numpy as np
import matplotlib.pyplot as plt

from SphToCart import SphToCart as SphToCart
from CartToSph import CartToSph as CartToSph
from ZeroTwoPi import ZeroTwoPi as ZeroTwoPi
from Stereonet import Stereonet as Stereonet
from GreatCircle import GreatCircle as GreatCircle
from StCoordLine import StCoordLine as StCoordLine


def Bingham(T,P):
    '''
    Bingham calculates and plots a cylindrical best fit to a distribution of
    lines to find fold axes from poles to bedding or the orientation of a 
    plane from lines on it. The statistical routine is based on algorithms in
    Fisher et al. (1988)

        USE: eigVec, confCone, bestFit = Bingham(T,P)
        
        T and P = Vectors of lines trends and plunges respectively
        
        eigVec = 3 x 3 matrix with eigenvalues (column 1), trends (column 2)
        and plunges (column 3) of the eigenvectors. Maximum eigenvalue and 
        corresponding eigenvector are in row 1, intermediate in row 2, 
        and minimum in row 3.

        confCone = 2 x 2 matrix with the maximum (column 1) and minimum 
        (column 2) radius of the 95% elliptical confidence cone around the 
        eigenvector corresponding to the largest (row 1), and lowest (row 2) 
        eigenvalue

        besFit = 1 x 2 vector containing the strike and dip (right hand rule) 
        of the best fit great circle to the distribution of lines

    NOTE: Input/Output trends and plunges, as well as confidence
    cones are in radians. Bingham plots the input lines, eigenvectors and
    best fit great circle in an equal area stereonet.
   
    Bingham uses functions ZeroTwoPi, Pole, StCoordLine SphToCart, 
    CartToSph, Rotate, GreatCircle, SmallCircle, GeogrToView, and Stereonet

    Python function translated from the Matlab function 
    Bingham in Allmendinger et al. (2012)
    '''

    # Some constants
    pi = np.pi
    east = pi/2
    twopi = pi*2

    # Number of lines
    nlines = len(T)

    # Initialize the orientation matrix
    a = np.zeros((3,3))

    # Fill the orientation matrix with the sums of the squares (for the 
    # principal diagonal) and the products of the direction cosines of each
    # line. cn, ce and cd are the north, east and down direction cosines
    for i in range(nlines): # evtl nlines-1
        cn,ce,cd = SphToCart(T[i],P[i],0)
        a[0,0] = a[0,0] + cn*cn
        a[0,1] = a[0,1] + cn*ce
        a[0,2] = a[0,2] + cn*cd
        a[1,1] = a[1,1] + ce*ce
        a[1,2] = a[1,2] + ce*cd
        a[2,2] = a[2,2] + cd*cd

    # The orientation matrix is symmetric so the off-diagonal components can be
    # equated
    a[1,0] = a[0,1]
    a[2,0] = a[0,2]
    a[2,1] = a[1,2]

    # Calculate the eigenvalues and eigenvectors of the orientation matrix using
    # Python function np.linalg.eig. D is a vector of eigenvalues and V is a 
    # full matrix whose columns are the corresponding eigenvectors
    D, V = np.linalg.eig(a)
    
    # Normalize the eigenvalues by the number of lines and convert the
    # corresponding eigenvectors to the lower hemisphere
    for i in range(3): 
        D[i] = D[i]/nlines
        if V[2,i] < 0:
            V[0,i] = -V[0,i]
            V[1,i] = -V[1,i]
            V[2,i] = -V[2,i]

    # Initialize eigVec
    eigVec = np.zeros((3,3))
    #Fill eigVec
    eigVec[0,0] = D[2]    # Maximum eigenvalue
    eigVec[1,0] = D[1]    # Intermediate eigenvalue
    eigVec[2,0] = D[0]    # Minimum eigenvalue
    # Trend and plunge of largest eigenvalue: column 3 of V
    eigVec[0,1], eigVec[0,2] = CartToSph(V[0,2], V[1,2], V[2,2])
    # Trend and plunge of intermediate eigenvalue: column 2 of V
    eigVec[1,1], eigVec[1,2] = CartToSph(V[0,1], V[1,1], V[2,1])
    # Trend and plunge of minimum eigenvalue: column 1 of V
    eigVec[2,1], eigVec[2,2] = CartToSph(V[0,0], V[1,0], V[2,0])
    
    # Initialize confCone
    confCone = np.zeros((2,2))
    # If there are more than 25 lines, calculate confidence cones at the 95%
    # confidence level. The algorithm is explained in Fisher et al. (1987)
    if nlines >= 25:
        e11 = 0
        e22 = 0
        e12 = 0
        d11 = 0
        d22 = 0
        d12 = 0
        en11 = 1/(nlines*(eigVec[2,0] - eigVec[0,0])**2)
        en22 = 1/(nlines*(eigVec[1,0] - eigVec[0,0])**2)
        en12 = 1/(nlines*(eigVec[2,0] - eigVec[0,0])*(eigVec[1,0] - eigVec[0,0]))
        dn11 = en11
        dn22 = 1/(nlines*(eigVec[2,0] - eigVec[1,0])**2);
        dn12 = 1/(nlines*(eigVec[2,0] - eigVec[1,0])*(eigVec[2,0] - eigVec[0,0]))
        vec = np.zeros((3,3))
        for i in range(3):
            vec[i,0] = np.sin(eigVec[i,2] + east)*np.cos(twopi - eigVec[i,1])
            vec[i,1] = np.sin(eigVec[i,2] + east)*np.sin(twopi - eigVec[i,1])
            vec[i,2] = np.sin(eigVec[i,2] + east)
        for i in range(nlines):
            c1 = np.sin(P[i]+east)*np.cos(twopi-T[i])
            c2 = np.sin(P[i]+east)*np.sin(twopi-T[i])
            c3 = np.cos(P[i]+east)
            u1x = vec[2,0]*c1 + vec[2,1]*c2 + vec[2,2]*c3
            u2x = vec[1,0]*c1 + vec[1,1]*c2 + vec[1,2]*c3
            u3x = vec[0,0]*c1 + vec[0,1]*c2 + vec[0,2]*c3
            e11 = u1x*u1x * u3x*u3x + e11
            e22 = u2x*u2x * u3x*u3x + e22
            e12 = u1x*u2x * u3x*u3x + e12
            d11 = e11
            d22 = u1x*u1x * u2x*u2x + d22
            d12 = u2x*u3x * u1x*u1x + d12
        e22 = en22*e22
        e11 = en11*e11
        e12 = en12*e12
        d22 = dn22*d22
        d11 = dn11*d11
        d12 = dn12*d12
        d = -2*np.log(.05)/nlines
        # initialize f
        f = np.zeros((2,2))
        if abs(e11*e22-e12*e12) >= 0.000001:
            f[0,0] = (1/(e11*e22-e12*e12)) * e22
            f[1,1] = (1/(e11*e22-e12*e12)) * e11
            f[0,1] = -(1/(e11*e22-e12*e12)) * e12
            f[1,0] = f[0,1]
            # Calculate the eigenvalues and eigenvectors of the matrix f using
            # Python function np.linaln.eig. The next lines follow steps 1-4 outlined on 
            # pp. 34-35 of Fisher et al. (1987)
            DD, _ = np.linalg.eig(f)
            if DD[0] > 0 and DD[1] > 0:
                if d/DD[0] <= 1 and d/DD[1] <= 1:
                    confCone[0,1] = np.arcsin(np.sqrt(d/DD[1]))
                    confCone[0,0] = np.arcsin(np.sqrt(d/DD[0]))
        # Repeat the process for the eigenvector corresponding to the smallest
        # eigenvalue
        if abs(d11*d22-d12*d12) >= 0.000001:
            f[0,0] = (1/(d11*d22-d12*d12)) * d22
            f[1,1] = (1/(d11*d22-d12*d12)) * d11
            f[0,1] = -(1/(d11*d22-d12*d12)) * d12
            f[1,0] = f[0,1]
            DD, _ = np.linalg.eig(f)
            if DD[0] > 0.0 and DD[1] > 0.0:
                if d/DD[0] <= 1 and d/DD[1] <= 1:
                    confCone[1,1] = np.arcsin(np.sqrt(d/DD[1]))
                    confCone[1,0] = np.arcsin(np.sqrt(d/DD[0]))

    # Calculate the best fit great circle to the distribution of points
    bestFit = np.zeros(2)
    bestFit[0] = ZeroTwoPi(eigVec[2,1] + east)
    bestFit[1] = east - eigVec[2,2]

    # Plot stereonet
    Stereonet(0, 90*pi/180, 10*pi/180, 1)
    
    # Plot lines
    for i in range(nlines):
        xp,yp = StCoordLine(T[i],P[i],1)
        plt.plot(xp,yp,'k.')
    
    # Plot eigenvectors
    for i in range(3):
        xp,yp = StCoordLine(eigVec[i,1],eigVec[i,2],1)
        plt.plot(xp,yp,'rs')
    
    # Plot best fit great circle
    path = GreatCircle(bestFit[0],bestFit[1],1)
    plt.plot(path[:,0],path[:,1],'r')
    
    # release plot
    plt.show()

    return eigVec, confCone, bestFit