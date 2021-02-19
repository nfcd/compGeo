import numpy as np
from SphToCart import SphToCart
from CartToSph import CartToSph
from ZeroTwoPi import ZeroTwoPi

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

	GeogrToView uses functions ZeroTwoPi, SphToCart and
	CartToSph

	Python function translated from the Matlab function
	GeogrToView in Allmendinger et al. (2012)
	'''
	#some constants 
	east = np.pi/2.0
	
	#Make transformation matrix between NED and View Direction
	a = np.zeros((3,3))
	a[2,0], a[2,1], a[2,2] = SphToCart(trdv,plgv,0)
	temp1 = trdv + east
	temp2 = 0.0
	a[1,0], a[1,1], a[1,2] = SphToCart(temp1,temp2,0)
	temp1 = trdv
	temp2 = plgv - east
	a[0,0], a[0,1], a[0,2] = SphToCart(temp1,temp2,0)
	
	#Direction cosines of line
	dirCos = np.zeros(3)
	dirCos[0], dirCos[1], dirCos[2] = SphToCart(trd,plg,0)
	
	# Transform line
	nDirCos = np.zeros(3)
	for i in range(0,3):
		nDirCos[i] = a[i,0]*dirCos[0] + a[i,1]*dirCos[1]+ a[i,2]*dirCos[2]
	
	# Compute line from new direction cosines
	rtrd, rplg = CartToSph(nDirCos[0],nDirCos[1],nDirCos[2])
	
	# Take care of negative plunges
	if rplg < 0.0 :
		rtrd = ZeroTwoPi(rtrd+np.pi)
		rplg = -rplg
	
	return rtrd, rplg