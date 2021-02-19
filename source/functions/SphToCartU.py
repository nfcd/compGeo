from uncertainties import umath

def SphToCartU(trd,plg,k):
	'''
	SphToCartU converts from spherical to cartesian coordinates

	SphToCartU(trd,plg,k) returns the north (cn),
	east (ce), and down (cd) direction cosines of a line.
	Notice that these direction cosines have uncertainties

	k is an integer to tell whether the trend and plunge
	of a line (k = 0) or strike and dip of a plane in right
	hand rule (k = 1) are being sent in the trdu and plgu slots.
	In this last case, the direction cosines of the pole to
	the plane are returned
 
	NOTE: trdu and plgu are in radians and they have
	uncertainties in radians
	
	SphToCartU uses the uncertainties package from
	Eric O. Lebigot
	
	Based on Python function SphToCart
	'''
	# If line
	if k == 0:
		cn = umath.cos(trd) * umath.cos(plg)
		ce = umath.sin(trd) * umath.cos(plg)
		cd = umath.sin(plg)
	# Else pole to plane
	elif k == 1:
		cn = umath.sin(trd) * umath.sin(plg)
		ce = -umath.cos(trd) * umath.sin(plg)
		cd = umath.cos(plg)
	
	return cn, ce, cd