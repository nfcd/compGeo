import numpy as np
from uncertainties import ufloat # From E. Lebigot
from uncertainties import umath # From E. Lebigot

def true_thickness_u(strike,dip,top,base):
	"""
	true_thickness_u calculates the thickness (t) of a unit
	given the strike (strike) and dip (dip) of the unit,
	and points at its top (top) and base (base).
	strike and dip, as well as the points have
	uncertainties.
	
	top and base are 1 x 3 arrays defining the location
	of top and base points in an ENU coordinate system.
	For each one of these arrays, the first, second
	and third entries are the E, N and U coordinates.
	These coordinates have uncertainties
	
	NOTE: strike and dip should be input in radians and
		they have uncertainties in radians. The
		returned thickness has also uncertainties
	"""
	# make the transformation matrix from ENU coordinates
	# to SDP coordinates
	sin_str = umath.sin(strike)
	cos_str = umath.cos(strike)
	sin_dip = umath.sin(dip)
	cos_dip = umath.cos(dip)
	a = np.array([[sin_str, cos_str, ufloat(0,0)],
	[cos_str*cos_dip, -sin_str*cos_dip, -sin_dip],
	[-cos_str*sin_dip, sin_str*sin_dip, -cos_dip]])
	
	# transform the top and base points
	# from ENU to SDP coordinates
	topn = np.array([ufloat(0,0), ufloat(0,0), ufloat(0,0)])
	basen = np.array([ufloat(0,0), ufloat(0,0), ufloat(0,0)])
	for i in range(0,3):
		for j in range(0,3):
			topn[i] = a[i,j]*top[j] + topn[i]
			basen[i] = a[i,j]*base[j] + basen[i]
	
	# compute the thickness of the unit
	t = np.abs(basen[2] - topn[2])
	
	return t