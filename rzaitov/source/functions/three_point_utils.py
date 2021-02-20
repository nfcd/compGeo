import numpy as np
from numpy import linalg as la
from cart_to_sph import cart_to_sph
from pole_utils import plane_from_pole

def three_point_utils(p1, p2, p3):
	'''
	ThreePoint calculates the strike (strike) and dip (dip)
	of a plane given the east (E), north (N), and up (U)
	coordinates of three non-collinear points on the plane

	p1, p2 and p3 are 1 x 3 arrays defining the location
	of the points in an ENU coordinate system. For each one
	of these arrays the first entry is the E coordinate,
	the second entry the N coordinate, and the third entry 
	the U coordinate

	NOTE: strike and dip are returned in radians and they
	follow the right-hand rule format
	'''
	# make vectors v (p1 - p3) and u (p2 - p3)
	v = p1 - p2
	u = p2 - p3

	# take the cross product of v and u
	vcu = np.cross(v, u)

	# make this vector a unit vector
	mvcu = la.norm(vcu) # magnitude of the vector
	if mvcu == 0:
		raise ValueError('Error: points are collinear')
	uvcu = vcu/mvcu # unit vector

	# make the pole vector in NED coordinates
	p = np.array([ uvcu[1], uvcu[0], -uvcu[2] ])

	# Make pole point downwards
	if p[2] < 0:
		p = -p

	# find the trend and plunge of the pole
	trd, plg = cart_to_sph(p[0], p[1], p[2])

	# find strike and dip of plane
	strike, dip = plane_from_pole(trd, plg)

	return strike, dip