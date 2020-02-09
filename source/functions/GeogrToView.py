import numpy as np
from SphToCart import SphToCart as SphToCart
from CartToSph import CartToSph as CartToSph
from ZeroTwoPi import ZeroTwopi as ZeroTwoPi

def GeogrToView(trd,plg,trdv,plgv):
    '''
	GeogrToView transforms a line from NED to View Direction
	coordinates
	trd = trend of line
   	plg = plunge of line 
   	trdv = trend of view direction 
   	plgv = plunge of view direction 
   	rtrd and rplg are the new trend and plunge of the line 
   	in the view direction.

   	NOTE: Input/Output angles are in radians

   	GeogrToView uses functions ZeroTwoPi, SphToCart and CartToSph

   	Python function translated from the Matlab function 
   	GeogrToView in Allmendinger et al. (2012)
	'''
    #some constants 
    east = np.pi/2.0
    
    #Make transformation matrix between NED and View Direction
    a = np.zeros((3,3))
    a[3,1], a[3,2], a[3,3] = SphToCart(trdv,plgv,0)
    temp1 = trdv + east
    temp2 = 0.0
    a[2,1], a[2,2], a[2,3] = SphToCart(temp1,temp2,0)
    temp1 = trdv
    temp2 = plgv - east
    a[1,1], a[1,2], a[1,3] = SphToCart(temp1,temp2,0)
     
     #Direction cosines of line
    dirCos = np.zeros((1,3))
    dirCos[1], dirCos[2], dirCos[3] = SphToCart(trd,plg,0)
     
     
     # Transform line
    nDirCos = np.zeros((1,3))
    for i in np.arange(1,3):
        nDirCos[i] = a[i,1]*dirCos[1] + a[i,2]*dirCos[2]+ a[i,3]*dirCos[3]
    # Compute line from new direction cosines
    rtrd, rplg = CartToSph(nDirCos[1],nDirCos[2],nDirCos[3])
    
    # Take care of negative plunges
    if rplg < 0.0 :
        rtrd = ZeroTwoPi(rtrd+np.pi)
        rplg = -rplg
    
    return rtrd, rplg