import math
from uncertainties import ufloat, umath # From E. Lebigot
from zero_twopi import zero_twopi

def cart_to_sph_u(cn,ce,cd):
	"""
	cart_to_sph_u converts from Cartesian to spherical coordinates
	
	cart_to_sph_u(cn,ce,cd) returns the trend (trd)
	and plunge (plg) of a line for input north (cn),
	east (ce), and down (cd) direction cosines.
	Input and output values have uncertainties
	
	NOTE: Trend and plunge are returned in radians
	
	Based on Python function cart_to_sph
	"""
	pi = math.pi
	# plunge 
	plg = umath.asin(cd) 
	
	# trend: If north direction cosine is zero, trend
	# is east or west. Choose which one by the sign of
	# the east direction cosine
	if cn.n == 0.0:
		if ce.n < 0.0: # use nominal value n to check
			trd = 3.0/2.0 * pi # trend is west
		else:
			trd = pi/2.0 # trend is east
	# else
	else:
		trd = umath.atan(ce/cn) 
		if cn.n < 0.0: # use nominal value n to check
			# add pi
			trd += pi 
		# make sure trend is between 0 and 2*pi
		trd = zero_twopi_u(trd)
	
	return trd, plg

def zero_twopi_u(a):
    """
    zero_twopi wrapper for UFloat inputs.
    Branching is done on nominal value only.
    """
    return ufloat(zero_twopi(a.n), a.std_dev)