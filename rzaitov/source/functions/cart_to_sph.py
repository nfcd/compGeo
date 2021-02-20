import numpy as np
from zero_two_pi import zero_two_pi

def cart_to_sph(cn, ce, cd):
	'''
	cart_to_sph converts from Cartesian to spherical coordinates

	cart_to_sph(cn, ce, cd) returns the trend
	and plunge of a line for input north (cn),
	east (ce), and down (cd) direction cosines

	NOTE: Trend and plunge are returned in radians

	Python function translated from the Matlab function
	CartToSph in Allmendinger et al. (2012)
	'''
	# convert scalar into 1d array https://numpy.org/doc/stable/reference/generated/numpy.atleast_1d.html
	cn, ce, cd = np.atleast_1d(cn, ce, cd)
	if cn.shape != ce.shape or ce.shape != cd.shape:
		raise ValueError("shapes of arguments are incompatible")

	pi = np.pi
	pi_2 = pi / 2.0
	plg = np.arcsin(cd) # Eq. 4.13a
	trd = np.zeros(plg.shape)

	# Trend: If north direction cosine is zero, trend
	# is east or west. Choose which one by the sign of
	# the east direction cosine

	mask = (cn == 0.0) * (ce <  0.0)
	trd[mask] = 3.0 * pi_2 # Eq. 4.14d, trend is west

	mask = (cn == 0.0) * (ce >= 0.0)
	trd[mask] = pi_2      # Eq. 4.14c, trend is east

	mask = (cn != 0.0)
	trd[mask] = np.arctan(ce[mask]/cn[mask]) # Eq. 4.14a

	mask = (cn < 0.0)
	trd[mask] = trd[mask] + pi # Eq. 4.14b

	# Make sure trend is between 0 and 2*pi
	trd = zero_two_pi(trd)

	if trd.shape[0] == 1:
		return trd.item(), plg.item()
	else:
		return trd, plg