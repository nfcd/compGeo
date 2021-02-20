import numpy as np
from sph_to_cart import sph_to_cart
from cart_to_sph import cart_to_sph
from rotation import rotation

def small_circle(axis_trd, axis_plg, cone_angle):
	'''
	small_circle computes the paths of a small circle defined
	by its axis and cone angle

	axis_trd = trend of axis
	axis_plg = plunge of axis
	cone_angle = cone angle

	Python function translated from the Matlab function
	SmallCircle in Allmendinger et al. (2012)
	'''
	v = np.array(sph_to_cart(axis_trd, axis_plg - cone_angle))  # vector to rotate

	# To make the small circle, rotate the starting line
	# 360 degrees in increments of 1 degree
	n = 361 # number of angles from 0 to 360
	NED = np.zeros((3, n))

	for i in range(n):
		rot = np.radians(i)
		R = rotation(axis_trd, axis_plg, rot)
		NED[:, i] = np.dot(R, v)

	# select vectors with non-negative plunge
	mask = NED[2, :] >= 0
	NED = NED[:, mask]

	return cart_to_sph(NED[0, :], NED[1, :], NED[2, :])