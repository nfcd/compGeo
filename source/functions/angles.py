import numpy as np
from sph_to_cart import sph_to_cart
from cart_to_sph import cart_to_sph
from pole import pole_from_plane, plane_from_pole

# Python functions based on the Matlab function
# Angles in Allmendinger et al. (2012)

def angle_bw_lines(trd1, plg1, trd2, plg2):
	"""
	angle_bw_lines returns the angle between two lines
	of trend and plunge trd1, plg1, trd2, and plg2
	Input and output angles are in radians
	"""
	# convert lines to directions cosines and numpy arrays
	cn1, ce1, cd1 = sph_to_cart(trd1, plg1)
	u = np.array([cn1, ce1, cd1])
	cn2, ce2, cd2 = sph_to_cart(trd2, plg2)
	v = np.array([cn2, ce2, cd2])
	# angle between lines is arccosine of their dot product
	return np.arccos(np.dot(u, v))

def angle_bw_planes(str1, dip1, str2, dip2):
	"""
	angle_bw_planes returns the angle between two planes
	of strike and dip str1, dip1, str2, and dip2
	Input and output angles are in radians
	"""
	# compute poles to lines
	pole1_trd, pole1_plg = pole_from_plane(str1, dip1)
	pole2_trd, pole2_plg = pole_from_plane(str2, dip2)
	# find angle between poles
	angle = angle_bw_lines(pole1_trd, pole1_plg, 
		pole2_trd, pole2_plg)
	# angle between planes is the complementary angle
	return (np.pi - angle)

def pole_from_lines(trd1, plg1, trd2, plg2):
	"""
	pole_from_lines compute the pole to a plane given
	two lines on the plane, with trend and plunge trd1, plg1,
	trd2, and plg2
	Input and output angles are in radians
	"""
	# convert lines to direction cosines and numpy arrays
	cn1, ce1, cd1 = sph_to_cart(trd1, plg1)
	u = np.array([cn1, ce1, cd1])
	cn2, ce2, cd2 = sph_to_cart(trd2, plg2)
	v = np.array([cn2, ce2, cd2])
	# normal is cross product between vectors
	pole = np.cross(u, v)
	# make pole a unit vector
	norm = np.linalg.norm(pole)
	pole = pole/norm
	# if pole points upwards, make it point downwards
	if pole[2] < 0:
		pole *= -1.0
	# return trend and plunge of pole
	return cart_to_sph(pole[0], pole[1], pole[2])

def plane_from_app_dips(trd1, plg1, trd2, plg2):
	"""
	plane_from_app_dips returns the strike and dip of a plane
	from two apparent dips with trend and plunge trd1, plg1,
	trd2, and plg2
	Input and output angles are in radians
	"""
	# Compute pole to plane from apparent dips (lines)
	pole_trd, pole_plg = pole_from_lines(trd1,plg1,trd2,plg2)
	# return strike and dip of plane
	return plane_from_pole(pole_trd, pole_plg)

def int_bw_planes(str1, dip1, str2, dip2):
	"""
	int_bw_planes returns the intersection between two planes
	of strike and dip str1, dip1, str2, dip2
	Input and output angles are in radians
	"""
	# compute poles to planes
	pole1_trd, pole1_plg = pole_from_plane(str1, dip1)
	pole2_trd, pole2_plg = pole_from_plane(str2, dip2)
	# intersection is normal to poles
	return pole_from_lines(pole1_trd, pole1_plg, pole2_trd, 
		pole2_plg)