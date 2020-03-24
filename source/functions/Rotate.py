import numpy as np
from SphToCart import SphToCart as SphToCart
from CartToSph import CartToSph as CartToSph

def Rotate(rtrd,rplg,rot,trd,plg,ans0):
    '''
	Rotate rotates a line by performing a coordinate 
	transformation. The algorithm was originally written 
	by Randall A. Marrett

	rtrd = trend of rotation axis
	rplg = plunge of rotation axis
	rot = magnitude of rotation
	trd = trend of the vector to be rotated
	plg = plunge of the vector to be rotated
	ans0 = A character indicating whether the line to be rotated 
	is an axis (ans0 = 'a') or a vector (ans0 = 'v')

	NOTE: All angles are in radians
	
	Rotate uses functions SphToCart and CartToSph
	
	Python function translated from the Matlab function 
	Rotate in Allmendinger et al. (2012)
	'''
    # Allocate some arrays
    a = np.zeros((3,3)) #'Transformation matrix
    pole = np.zeros(3) #'Direction cosines of rotation axis
    plotr = np.zeros(3) #'Direction cosines of rotated vector
    temp = np.zeros(3)  #'Direction cosines of unrotated vector
    
    # Convert rotation axis to direction cosines. Note that the 
    # convention here is X1 = North, X2 = East, X3 = Down
    pole[1] , pole[2], pole[3] = SphToCart(rtrd,rplg,0)
    
    # Calculate the transformation matrix a for the rotation
    # Eq. 5.17
    x = 1.0 - np.cos(rot)
    sinRot = np.sin(rot)
    cosRot = np.cos(rot)
    a[1,1] = cosRot + pole[1]*pole[1]*x
    a[1,2] = -pole[3]*sinRot + pole[1]*pole[2]*x
    a[1,3] = pole[2]*sinRot + pole[1]*pole[3]*x
    a[2,1] = pole[3]*sinRot + pole[2]*pole[1]*x
    a[2,2] = cosRot + pole[2]*pole[2]*x
    a[2,3] = -pole[1]*sinRot + pole[2]*pole[3]*x
    a[3,1] = -pole[2]*sinRot + pole[3]*pole[1]*x
    a[3,2] = pole[1]*sinRot + pole[3]*pole[2]*x
    a[3,3] = cosRot + pole[3]*pole[3]*x
    
    
    # Convert trend and plunge of vector to be rotated into 
    # direction cosines
    temp[1] , temp[2], temp[3] = SphToCart(trd,plg,0)
    
    # Perform the coordinate transformation
    for i in range(0,3,1):
        plotr[i] = 0.0
        for j in range(0,3,1):
            plotr[i] = a[i,j]*temp[j] + plotr[i]
            
    # Convert to lower hemisphere projection if data are
    # axes (ans0 = 'a')
    if plotr[3] < 0.0 and ans0 == 'a' :
        plotr[1] = -plotr[1]
        plotr[2] = -plotr[2]
        plotr[3] = -plotr[3]
        
    # Convert from direction cosines back to trend and plunge
    rtrd , rplg = CartToSph(plotr[1], plotr[2], plotr[3])
    
    return rtrd, rplg