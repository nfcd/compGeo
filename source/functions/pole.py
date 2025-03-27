import math
from zero_twopi import zero_twopi

def pole_from_plane(stk,dip):
	"""
	pole_from_plane returns the trend (trd) and 
	plunge (plg) of a pole, given the strike and 
	dip of the plane
	
	NOTE: Input/Output angles are in radians.
	Input stk and dip is in RHR format
	"""
	# some constants
	east = math.pi/2
	
	# pole from plane
	trd = zero_twopi (stk - east)
	plg = east - dip
	
	return trd, plg

def plane_from_pole(trd,plg):
	"""
	plane_from_pole returns the strike and dip
	of a plane, given the trend (trd) and 
	plunge (plg) of its pole
	
	NOTE: Input/Output angles are in radians.
	Output stk and dip is in RHR format
	"""
	# some constants
	pi = math.pi
	east = pi/2
	
	# unusual case of pole pointing upwards
	if plg < 0.0:
		trd += pi
		plg *= -1.0
		
	# calculate plane given its pole
	stk = zero_twopi(trd + east)
	dip = east - plg
		
	return stk, dip