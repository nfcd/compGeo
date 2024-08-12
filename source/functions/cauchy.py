import numpy as np

from dircos_axes import dircos_axes
from sph_to_cart import sph_to_cart
from pole import pole_from_plane

def cauchy(stress,tx1,px1,tx3,stk,dip):
	"""
	Given the stress tensor in a X1X2X3 coordinate system,
	cauchy computes the X1X2X3 tractions on an arbitrarily
	oriented plane 
	
	USE: t,pt = cauchy(stress,tx1,px1,tx3,stk,dip)
	
	stress = 3 x 3 stress tensor
	tx1 = trend of X1
	px1 = plunge of X1
	tx3 = trend of X3
	stk = strike of plane
	dip = dip of plane
	t = 1 x 3 vector with X1, X2 and X3 tractions
	pt = 1 x 3 vector with direction cosines of pole
		to plane with respect to X1X2X3
	
	NOTE = Plane orientation follows the right hand rule
	Input/Output angles are in radians
	
	Python function translated from the Matlab function
	Cauchy in Allmendinger et al. (2012)
	"""
	# compute direction cosines of X1X2X3 with respect
	# to NED
	dc = dircos_axes(tx1,px1,tx3)
	
	# calculate direction cosines of pole to plane
	trd, plg = pole_from_plane(stk,dip)
	p = np.zeros(3)
	p[0],p[1],p[2] = sph_to_cart(trd,plg)
	
	# transform pole to plane to stress coordinates X1X2X3
	# The transformation matrix is the direction cosines of
	# X1X2X3
	pt = np.zeros(3)
	for i in range(3):
		for j in range(3):
			pt[i] = dc[i,j]*p[j] + pt[i]
			
	# calculate the tractions in stress coordinates X1X2X3
	t = np.zeros(3)
	# compute tractions using Cauchy"s law
	for i in range(3):
		for j in range(3):
			t[i] = stress[i,j]*pt[j] + t[i]
	
	return t, pt