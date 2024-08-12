import numpy as np
from cart_to_sph import cart_to_sph
from pole import plane_from_pole

def three_points(p1,p2,p3):
	"""
	three_points calculates the strike (stk) and dip (dip)
	of a plane given the east (E), north (N), and up (U)
	coordinates of three non-collinear points on the plane
	
	p1, p2 and p3 are 1 x 3 arrays defining the location
	of the points in an ENU coordinate system. For each one
	of these arrays the first, second and third entries are 
	the E, N and U coordinates, respectively
	
	NOTE: stk and dip are returned in radians and they
	follow the right-hand rule format
	"""
	# make vectors v (p1 - p3) and u (p2 - p3)
	v = p1 - p2
	u = p2 - p3
	# take the cross product of v and u
	vcu = np.cross(v,u)
	
	# make this vector a unit vector
	mvcu = np.linalg.norm(vcu) # magnitude of the vector
	if mvcu == 0: # If points collinear
		raise ValueError("Error: points are collinear")
	
	vcu = vcu/mvcu # unit vector
	
	# make the pole vector in NED coordinates
	p = np.array([vcu[1], vcu[0], -vcu[2]])
	
	# make pole point downwards
	if p[2] < 0:
		p *= -1.0
	
	# find the trend and plunge of the pole
	trd, plg = cart_to_sph(p[0],p[1],p[2])
	
	# find stk and dip of plane
	stk, dip = plane_from_pole(trd, plg)
	
	return stk, dip