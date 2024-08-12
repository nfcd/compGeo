import numpy as np
from uncertainties import umath # From E. Lebigot
from sph_to_cart_u import sph_to_cart_u
from cart_to_sph_u import cart_to_sph_u
from pole import pole_from_plane, plane_from_pole

# these functions are similar to the ones in the module
# angles.py. However, in this case input and output
# values have uncertainties

def angle_bw_lines_u(trd1, plg1, trd2, plg2):
	"""
	angle_bw_lines_u returns the angle between two lines
	of trend and plunge trd1, plg1, trd2, and plg2
	Input and output angles are in radians
	"""
	# convert lines to directions cosines and numpy arrays
	cn1, ce1, cd1 = sph_to_cart_u(trd1, plg1)
	u = np.array([cn1, ce1, cd1])
	cn2, ce2, cd2 = sph_to_cart_u(trd2, plg2)
	v = np.array([cn2, ce2, cd2])
	# angle between lines is arccosine of their dot product
	return umath.acos(np.dot(u, v))

def angle_bw_planes_u(stk1, dip1, stk2, dip2):
	"""
	angle_bw_planes_u returns the angle between two planes
	of strike and dip stk1, dip1, stk2, and dip2
	Input and output angles are in radians
	"""
	# compute poles to lines
	p1_trd, p1_plg = pole_from_plane(stk1, dip1)
	p2_trd, p2_plg = pole_from_plane(stk2, dip2)
	# find angle between poles
	angle = angle_bw_lines_u(p1_trd,p1_plg,p2_trd,p2_plg)
	# angle between planes is the complementary angle
	return (np.pi - angle)

def pole_from_lines_u(trd1, plg1, trd2, plg2):
	"""
	pole_from_lines_u compute the pole to a plane given
	two lines on the plane, with trend and plunge trd1, plg1,
	trd2, and plg2
	Input and output angles are in radians
	"""
	# convert lines to direction cosines and numpy arrays
	cn1, ce1, cd1 = sph_to_cart_u(trd1, plg1)
	u = np.array([cn1, ce1, cd1])
	cn2, ce2, cd2 = sph_to_cart_u(trd2, plg2)
	v = np.array([cn2, ce2, cd2])
	# normal is cross product between vectors
	p = np.cross(u, v)
	# make pole a unit vector by dividing it
	# by its magnitude
	norm = umath.sqrt(np.dot(p, p))
	p = p/norm
	# if pole points upwards, make it point downwards
	if p[2] < 0:
		p *= -1.0
	# return trend and plunge of pole
	return cart_to_sph_u(p[0], p[1], p[2])

def plane_from_app_dips_u(trd1, plg1, trd2, plg2):
	"""
	plane_from_app_dips_u returns the strike and dip of a plane
	from two apparent dips with trend and plunge trd1, plg1,
	trd2, and plg2
	Input and output angles are in radians
	"""
	# compute pole to plane from apparent dips (lines)
	p_trd, p_plg = pole_from_lines_u(trd1, plg1, trd2, plg2)
	# return strike and dip of plane
	return plane_from_pole(p_trd, p_plg)

def int_bw_planes_u(stk1, dip1, stk2, dip2):
	"""
	int_bw_planes_u returns the intersection between two planes
	of strike and dip stk1, dip1, stk2, dip2
	Input and output angles are in radians
	"""
	# compute poles to planes
	p1_trd, p1_plg = pole_from_plane(stk1, dip1)
	p2_trd, p2_plg = pole_from_plane(stk2, dip2)
	# intersection is normal to poles
	return pole_from_lines_u(p1_trd,p1_plg,p2_trd,p2_plg)