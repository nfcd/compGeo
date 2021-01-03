import math

def SphToCart(trd,plg,k):
	'''
	SphToCart converts from spherical to Cartesian coordinates

	SphToCart(trd,plg,k) returns the north (cn),
	east (ce), and down (cd) direction cosines of a line.
	
	k: integer to tell whether the trend and plunge of a line
	(k = 0) or strike and dip of a plane in right hand rule
	(k = 1) are being sent in the trd and plg slots. In this
	last case, the direction cosines of the pole to the plane
	are returned

	NOTE: Angles should be entered in radians
	
	Python function translated from the Matlab function
	SphToCart in Allmendinger et al. (2012)
	'''
	# If line (see Table 4.1)
	if k == 0:
		cn = math.cos(trd) * math.cos(plg)
		ce = math.sin(trd) * math.cos(plg)
		cd = math.sin(plg)
	# Else pole to plane (see Table 4.1)
	elif k == 1:
		cn = math.sin(trd) * math.sin(plg)
		ce = -math.cos(trd) * math.sin(plg)
		cd = math.cos(plg)
	
	return cn, ce, cd