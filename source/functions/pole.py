import math
from zero_twopi import zero_twopi

def pole_from_plane(strike,dip):
	'''
	pole_from_plane returns the trend (trd) and 
	plunge (plg) of a pole, given the strike and 
	dip of the plane
	
	NOTE: Input/Output angles are in radians.
	Input strike and dip is in RHR format
	'''
	# Some constants
	east = math.pi/2
	
	# Pole from plane
	trd = zero_twopi (strike - east)
	plg = east - dip
	
	return trd, plg

def plane_from_pole(trd,plg):
	'''
	plane_from_pole returns the strike and dip
	of a plane, given the trend (trd) and 
	plunge (plg) of its pole
	
	NOTE: Input/Output angles are in radians.
	Output strike and dip is in RHR format
	'''
	# Some constants
	pi = math.pi
	east = pi/2
	
	# Unusual case of pole pointing upwards
	if plg < 0.0:
		trd += pi
		plg *= -1.0
		
	# Calculate plane given its pole
	strike = zero_twopi(trd + east)
	dip = east - plg
		
	return strike, dip