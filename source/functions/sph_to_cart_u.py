from uncertainties import umath # From E. Lebigot

def sph_to_cart_u(trd,plg):
	"""
	sph_to_cart_u converts line from spherical (trend 
	and plunge) to Cartesian (direction cosines) 
	coordinates
	
	sph_to_cart_u(trd,plg) returns the north (cn),
	east (ce), and down (cd) direction cosines of 
	a line with trend = trd and plunge = plg
	Input and output values have uncertainties
	
	NOTE: Angles should be entered in radians
	
	Based on Python function sph_to_cart
	"""
	# Compute direction cosines from trend and plunge
	cn = umath.cos(trd) * umath.cos(plg)
	ce = umath.sin(trd) * umath.cos(plg)
	cd = umath.sin(plg)
	
	return cn, ce, cd