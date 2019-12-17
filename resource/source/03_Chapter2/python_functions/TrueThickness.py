import numpy as np

def TrueThickness(strike, dip, top, base):
    '''
    TrueThickness calculates the thickness (t) of a unit given the
    strike (strike) and dip (dip) of the unit, and points at its top (top)
    and base (base) 

    top and base are 1 x 3 vectors defining the location of top and base
    points in an ENU coordinate system. For each of these vectors, the 
    first entry is the E coordinate, the second entry the N coordinate, 
    and the third entry the U coordinate

    NOTE: strike and dip should be input in radians
    
    Python function translated from the Matlab function TrueThickness in the book:
    Structural Geology Algorithms by Allmendinger, Cardozo, & Fisher, 2012
    '''
    # make the transformation matrix from ENU coordinates to SDP coordinates
    # S = strike direction, D = dip direction, P = pole direction
    a = [[np.sin(strike), np.cos(strike), 0],
    [np.cos(strike)*np.cos(dip), -np.sin(strike)*np.cos(dip), -np.sin(dip)],
    [-np.cos(strike)*np.sin(dip), np.sin(strike)*np.sin(dip), -np.cos(dip)]]
    
    # transform the top and base points from ENU to SDP coordinates 
    topn = np.zeros(3)
    basen = np.zeros(3)
    for i in range(0,3,1):
        for j in range(0,3,1):
            topn[i] = a[i][j]*top[j] + topn[i]
            basen[i] = a[i][j]*base[j] + basen[i]
    
    # compute the thickness of the unit
    t = np.abs(topn[2] - basen[2])
    return t