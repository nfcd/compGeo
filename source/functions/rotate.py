import numpy as np
from sph_to_cart import sph_to_cart
from cart_to_sph import cart_to_sph

def rotate(rtrd,rplg,rot,trd,plg,ans0):
	"""
	rotate rotates a line by performing a coordinate
	transformation
	
	rtrd = trend of rotation axis
	rplg = plunge of rotation axis
	rot = magnitude of rotation
	trd = trend of the line to be rotated
	plg = plunge of the line to be rotated
	ans0 = A character indicating whether the line 
		to be rotated is an axis (ans0 = "a") or a 
		vector (ans0 = "v")
	trdr and plgr are the trend and plunge of the
		rotated line
	
	NOTE: All angles are in radians
	
	Python function translated from the Matlab function
	Rotate in Allmendinger et al. (2012)
	"""
	# allocate some arrays
	a = np.zeros((3,3)) #Transformation matrix
	raxis = np.zeros(3) #Dir. cosines of rotation axis
	line = np.zeros(3)	#Dir. cosines of line to be rotated
	liner = np.zeros(3) #Dir. cosines of rotated line
	
	# convert rotation axis to direction cosines
	raxis[0] , raxis[1], raxis[2] = sph_to_cart(rtrd,rplg)
	
	# calculate the transformation matrix a for the rotation
	x = 1.0 - np.cos(rot)
	sinrot = np.sin(rot)
	cosrot = np.cos(rot)
	a[0,0] = cosrot + raxis[0]*raxis[0]*x
	a[0,1] = -raxis[2]*sinrot + raxis[0]*raxis[1]*x
	a[0,2] = raxis[1]*sinrot + raxis[0]*raxis[2]*x
	a[1,0] = raxis[2]*sinrot + raxis[1]*raxis[0]*x
	a[1,1] = cosrot + raxis[1]*raxis[1]*x
	a[1,2] = -raxis[0]*sinrot + raxis[1]*raxis[2]*x
	a[2,0] = -raxis[1]*sinrot + raxis[2]*raxis[0]*x
	a[2,1] = raxis[0]*sinrot + raxis[2]*raxis[1]*x
	a[2,2] = cosrot + raxis[2]*raxis[2]*x
	
	# convert trend and plunge of line to be rotated into
	# direction cosines
	line[0] , line[1], line[2] = sph_to_cart(trd,plg)
	
	# perform the coordinate transformation
	for i in range(3):
		for j in range(3):
			liner[i] = a[i,j]*line[j] + liner[i]
	
	# make sure the rotated line is a unit vector
	norm = np.linalg.norm(liner)
	liner = liner/norm
	
	# convert to lower hemisphere projection if axis 
	if liner[2] < 0.0 and ans0 == 'a':
		liner *= -1.0
		
	# convert from direction cosines back to trend and plunge
	trdr , plgr = cart_to_sph(liner[0], liner[1], liner[2])
	
	return trdr, plgr