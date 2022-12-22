import math

def sph_to_cart(trd,plg):
	"""
	sph_to_cart converts line from spherical (trend 
	and plunge) to Cartesian (direction cosines) 
	coordinates
	
	sph_to_cart(trd,plg) returns the north (cn),
	east (ce), and down (cd) direction cosines of 
	a line with trend = trd and plunge = plg
	
	NOTE: Angles should be entered in radians
	
	Python function based on the Matlab function
	SphToCart in Allmendinger et al. (2012)
	"""
	# Compute direction cosines from trend and plunge
	cn = math.cos(trd) * math.cos(plg)
	ce = math.sin(trd) * math.cos(plg)
	cd = math.sin(plg)
	
	return cn, ce, cd