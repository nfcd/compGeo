import numpy as np
from sph_to_cart import sph_to_cart
from cart_to_sph import cart_to_sph
from zero_twopi import zero_twopi

def geogr_to_view(trd,plg,trdv,plgv):
	"""
	geogr_to_view transforms a line from NED to View 
	direction coordinates
	trd = trend of line
	plg = plunge of line 
	trdv = trend of view direction 
	plgv = plunge of view direction 
	rtrd and rplg are the new trend and plunge of the line
	in the view direction.
	
	NOTE: Input/Output angles are in radians
	
	Python function translated from the Matlab function
	GeogrToView in Allmendinger et al. (2012)
	"""
	#some constants 
	east = np.pi/2.0
	
	#Make transformation matrix between NED and View Direction
	a = np.zeros((3,3))
	a[2,0], a[2,1], a[2,2] = sph_to_cart(trdv,plgv)
	temp1 = trdv + east
	temp2 = 0.0
	a[1,0], a[1,1], a[1,2] = sph_to_cart(temp1,temp2)
	temp1 = trdv
	temp2 = plgv - east
	a[0,0], a[0,1], a[0,2] = sph_to_cart(temp1,temp2)
	
	#Direction cosines of line
	dc = np.zeros(3)
	dc[0], dc[1], dc[2] = sph_to_cart(trd,plg)
	
	# Transform line
	ndc = np.zeros(3)
	for i in range(3):
		ndc[i] = a[i,0]*dc[0] + a[i,1]*dc[1]+ a[i,2]*dc[2]
	
	# Compute line from new direction cosines
	rtrd, rplg = cart_to_sph(ndc[0],ndc[1],ndc[2])
	
	# Take care of negative plunges
	if rplg < 0.0 :
		rtrd = zero_twopi(rtrd+np.pi)
		rplg = -rplg
	
	return rtrd, rplg