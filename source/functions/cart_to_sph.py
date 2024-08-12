import math
from zero_twopi import zero_twopi

def cart_to_sph(cn,ce,cd):
	"""
	cart_to_sph converts from Cartesian to spherical coordinates
	
	cart_to_sph(cn,ce,cd) returns the trend (trd)
	and plunge (plg) of a line with north (cn),
	east (ce), and down (cd) direction cosines
	
	NOTE: Trend and plunge are returned in radians
	
	Python function translated from the Matlab function
	CartToSph in Allmendinger et al. (2012)
	"""
	pi = math.pi
	# plunge 
	plg = math.asin(cd) 
	
	# trend: If north direction cosine is zero, trend
	# is east or west. Choose which one by the sign of
	# the east direction cosine
	if cn == 0.0:
		if ce < 0.0:
			trd = 3.0/2.0 * pi # trend is west
		else:
			trd = pi/2.0 # trend is east
	else:
		trd = math.atan(ce/cn) 
		if cn < 0.0:
			trd += pi 
		# make sure trend is between 0 and 2*pi
		trd = zero_twopi(trd)
	
	return trd, plg