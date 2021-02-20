import math
from sph_to_cart import sph_to_cart
from cart_to_sph import cart_to_sph
from pole_utils import pole_from_plane, plane_from_pole

def angle_from_lines(trd1, plg1, trd2, plg2):
	cn1, ce1, cd1 = sph_to_cart(trd1, plg1)
	cn2, ce2, cd2 = sph_to_cart(trd2, plg2)
	return math.acos(cn1*cn2 + ce1*ce2 + cd1*cd2)

def angle_from_planes(str1, dip1, str2, dip2):
	pole_trend_1, pole_plunge_1 = pole_from_plane(str1, dip1)
	pole_trend_2, pole_plunge_2 = pole_from_plane(str2, dip2)
	return angle_from_lines(pole_trend_1, pole_plunge_1, pole_trend_2, pole_plunge_2)

def plane_from_lines(trd1, plg1, trd2, plg2):
	pole_trd, pole_plg = orthogonal_line(trd1, plg1, trd2, plg2)
	return plane_from_pole(pole_trd, pole_plg)

def line_from_planes(str1, dip1, str2, dip2):
	pole_trend_1, pole_plunge_1 = pole_from_plane(str1, dip1)
	pole_trend_2, pole_plunge_2 = pole_from_plane(str2, dip2)
	return orthogonal_line(pole_trend_1, pole_plunge_1, pole_trend_2, pole_plunge_2)

# find orthogonal line to given two lines
def orthogonal_line(trd1, plg1, trd2, plg2):
	cn1, ce1, cd1 = sph_to_cart(trd1, plg1)
	cn2, ce2, cd2 = sph_to_cart(trd2, plg2)
	# find orthogonal vector by calc cross product by hand
	cn = ce1*cd2 - cd1*ce2
	ce = cd1*cn2 - cn1*cd2
	cd = cn1*ce2 - ce1*cn2
	if cd < 0:
		cn = -cn
		ce = -ce
		cd = -cd
	# cross product of two unit vectors are not necessary unit, so norm it
	r = math.sqrt(cn*cn+ce*ce+cd*cd)
	return cart_to_sph(cn/r, ce/r, cd/r)