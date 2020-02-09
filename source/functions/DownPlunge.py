import numpy as np

def DownPlunge(bedseg, trd, plg):
    '''
    DownPlunge constructs the down plunge projection
    of a bed
    
    bedseg is a npoints x 3 array, which holds npoints 
    on the digitized bed, each point defined by
    3 coordinates: X1 = East, X2 = North, X3 = Up
    
    trd and plg are the trend and plunge of the fold axis
    
    NOTE: trd and plg should be entered in radians
    
    Python function translated from the Matlab function 
    DownPlunge in Allmendinger et al. (2012)
    '''
    # Number of points in bed
    nvtex = bedseg.np.shape[0]
    
    # Allocate some arrays
    a=np.zeros((3,3))
    dpbedseg = np.zeros(np.shape(bedseg))
    
    # Calculate the transformation matrix a(i,j)
    a[1,1] = np.sin(trd)*np.sin(plg)
    a[1,2] = np.cos(trd)*np.sin(plg)
    a[1,3] = np.cos(plg)
    a[2,1] = np.cos(trd)
    a[2,2] = -np.sin(trd)
    a[2,3] = 0.0
    a[3,1] = np.sin(trd)*np.cos(plg)
    a[3,2] = np.cos(trd)*np.cos(plg)
    a[3,3] = -np.sin(plg)
    
    # Transform the points' coordinates
    for nv in np.arange(1,nvtex):
        for i in np.arange(1,3):
            dpbedseg[nv,i] = 0.0
            for j in np.arange(1,3):
                dpbedseg[nv,i] = a[i,j]*bedseg[nv,j] + dpbedseg[nv,i]
                
    return dpbedseg