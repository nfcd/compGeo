import numpy as np
from numpy import linalg as la
from CartToSph import CartToSph
from Pole import Pole

def ThreePoint(p1,p2,p3):
	'''
	ThreePoint calculates the strike (strike) and dip (dip)
	of a plane given the east (E), north (N), and up (U)
	coordinates of three non-collinear points on the plane
	
	p1, p2 and p3 are 1 x 3 arrays defining the location
	of the points in an ENU coordinate system. For each one
	of these arrays the first, second and third entries are 
	the E, N and U coordinates, respectively

	NOTE: strike and dip are returned in radians and they
	follow the right-hand rule format

	ThreePoint uses functions CartToSph and Pole
	'''
	# make vectors v (p1 - p3) and u (p2 - p3)
	v = p1 - p2
	u = p2 - p3
	# take the cross product of v and u
	vcu = np.cross(v,u)
	
	# make this vector a unit vector
	mvcu = la.norm(vcu) # magnitude of the vector
	if mvcu == 0: # If points collinear
		raise ValueError('Error: points are collinear')
	
	uvcu = vcu/mvcu # unit vector
	
	# make the pole vector in NED coordinates
	p = [uvcu[1], uvcu[0], -uvcu[2]]

	# Make pole point downwards
	if p[2] < 0:
		p = [-elem for elem in p]
	
	# find the trend and plunge of the pole
	trd, plg = CartToSph(p[0],p[1],p[2])
	
	# find strike and dip of plane
	strike, dip = Pole(trd, plg, 0)
	
	return strike, dip