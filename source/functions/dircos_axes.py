import numpy as np

from sph_to_cart import sph_to_cart

def dircos_axes(tx1,px1,tx3):
	"""
	dircos_axes calculates the direction cosines of a right
	handed, orthogonal X1X2X3 cartesian coordinate system
	of any orientation with respect to NED
	
	USE: dc = dircos_axes(tx1,px1,tx3)
	
	tx1 = trend of X1
	px1 = plunge of X1
	tx3 = trend of X3
	dc = 3 x 3 matrix containing the direction cosines
		of X1 (row 1), X2 (row 2), and X3 (row 3)
	
	Note: Input angles should be in radians
	
	Python function translated from the Matlab function
	DirCosAxes in Allmendinger et al. (2012)
	"""
	# some constants
	east = np.pi/2.0
	west = 3.0*east
	# tolerance for near zero values
	tol = 1e-6 
	
	# initialize matrix of direction cosines
	dc = np.zeros((3,3))
	
	# direction cosines of X1
	dc[0,0],dc[0,1],dc[0,2] = sph_to_cart(tx1,px1)
	
	# calculate plunge of axis 3
	# if axis 1 is horizontal
	if abs(px1) < tol:
		dt = abs(tx1-tx3)
		if abs(dt - east) < tol or abs(dt - west) < tol:
			px3 = 0.0
		else:
			px3 = east
	# else
	else:
		# since dot product X1 and X3 = 0
		px3 = np.arctan(-(dc[0,0]*np.cos(tx3)
			+ dc[0,1]*np.sin(tx3))/dc[0,2])
	
	# direction cosines of X3
	dc[2,0],dc[2,1],dc[2,2] = sph_to_cart(tx3,px3)
	
	# X2 is the cross product of X3 and X1
	dc[1,:] = np.cross(dc[2,:],dc[0,:])
	# make it a unit vector
	dc[1,:] = dc[1,:]/np.linalg.norm(dc[1,:])
	
	return dc