import numpy as np
from sph_to_cart import sph_to_cart
from cart_to_sph import cart_to_sph
from zero_two_pi import zero_two_pi

#some constants
pi_2 = np.pi / 2.0

def geogr_to_view(trd, plg, trdv, plgv):
	'''
	geogr_to_view transforms a line from NED to View Direction
	coordinates
	trd = trend of line
	plg = plunge of line 
	trdv = trend of view direction 
	plgv = plunge of view direction 
	rtrd and rplg are the new trend and plunge of the line
	in the view direction.

	NOTE: Input/Output angles are in radians

	Python function translated from the Matlab function
	GeogrToView in Allmendinger et al. (2012)
	'''
	# If view direction is the default (trdv=0, plgv=90) play as an identity transform
	if trdv == 0.0 and plgv == pi_2:
		return trd, plg

	# Make transformation matrix between NED and View Direction
	T = np.array([
		sph_to_cart(trdv, plgv - pi_2),
		sph_to_cart(trdv + pi_2, 0.0),
		sph_to_cart(trdv, plgv)
	])

	l = np.array(sph_to_cart(trd, plg)) # Direction cosines of line
	t_cn, t_ce, t_cd = np.dot(T, l) # Transform line
	rtrd, rplg = cart_to_sph(t_cn, t_ce, t_cd) # Line from new direction cosines

	# Take care of negative plunges
	if rplg < 0.0:
		rtrd = zero_two_pi(rtrd+np.pi)
		rplg = -rplg

	return rtrd, rplg