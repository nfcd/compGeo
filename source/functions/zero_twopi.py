import math

def zero_twopi(a):
	"""
	This function makes sure input azimuth (a)
	is within 0 and 2*pi
	
	NOTE: Azimuth a is input/output in radians
	
	Python function translated from the Matlab function
	ZeroTwoPi in Allmendinger et al. (2012)
	"""
	twopi = 2*math.pi
	if a < 0:
		a += twopi
	elif a >= twopi:
		a -= twopi
	
	return a