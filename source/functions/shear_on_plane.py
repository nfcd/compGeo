import numpy as np
from principal_stress import principal_stress
from sph_to_cart import sph_to_cart
from cart_to_sph import cart_to_sph
from pole import pole_from_plane

def shear_on_plane(stress,tx1,px1,tx3,stk,dip):
	"""
	shear_on_plane calculates the direction and magnitudes of
	the normal and maximum shear tractions on an arbitrarily
	oriented plane
	
	tt,dctt,srat=shear_on_plane(stress,tx1,px1,tx3,stk,dip)
	
	stress = 3 x 3 stress tensor
	tx1 = trend of X1
	px1 = plunge of X1
	tx3 = trend of X3
	stk = strike of plane
	dip = dip of plane
	tt = 3 x 3 matrix with the magnitude (column 1),
		trend (column 2) and plunge (column 3) of the:
		normal traction on the plane (row 1), zero shear
		traction (row 2), and max. shear traction (row 3)
	dctt = 3 x 3 matrix with the direction cosines of unit
		vectors parallel to: normal traction on the plane
		(row 1), zero shear traction (row 2), and maximum
		shear traction (row 3) with respect to NED
	srat = stress ratio
	
	NOTE = Input stress tensor does not need to be along
		principal stress directions
		Plane orientation follows the right hand rule
		Input/Output angles are in radians
	
	Python function translated from the Matlab function
	shear_on_plane in Allmendinger et al. (2012)
	"""
	# compute principal stresses and their orientations
	pstress, dcp = principal_stress(stress,tx1,px1,tx3)
	
	# update stress tensor to principal stress directions
	stress = np.zeros((3,3))
	for i in range(3):
		stress[i,i] = pstress[i,0]
	
	# calculate stress ratio
	srat = (stress[1,1]-stress[0,0])/(stress[2,2]-stress[0,0])
	
	# calculate direction cosines of pole to plane
	trd, plg = pole_from_plane(stk,dip)
	p = np.zeros(3)
	p[0], p[1], p[2] = sph_to_cart(trd,plg)
	
	# transform pole to plane to  principal stress coordinates
	pt = np.zeros(3)
	for i in range(3):
		for j in range(3):
			pt[i] = dcp[i,j]*p[j] + pt[i]
	
	# calculate the tractions in principal stress coordinates
	t = np.zeros(3)
	# compute tractions using Cauchy's law
	for i in range(3):
		for j in range(3):
			t[i] = stress[i,j]*pt[j] + t[i]
	
	# find the b axis by the cross product of t and pt
	b = np.cross(t,pt)
	
	# find the max. shear traction orientation
	# by the cross product of pt and b
	s = np.cross(pt,b)
	
	# convert b and s to unit vectors
	b = b/np.linalg.norm(b)
	s = s/np.linalg.norm(s)
	
	# now we can write the transformation matrix from
	# principal stress coordinates to plane coordinates
	a = np.zeros((3,3))
	a[0,:] = pt
	a[1,:] = b
	a[2,:] = s
	
	# transform stress from principal to plane coordinates
	# do it only for the first row since we are just
	# interested in the plane: sigma'11 = normal traction,
	# sigma'12 = zero, and sigma'13 = max. shear traction
	tt = np.zeros((3,3))
	for i in range(3):
		tt[i,0] = stress[0,0]*a[0,0]*a[i,0] \
		+ stress[1,1]*a[0,1]*a[i,1] + stress[2,2]*a[0,2]*a[i,2]
	
	# transform normal traction, zero shear
	# and max. shear traction to NED coords
	dctt = np.zeros((3,3))
	for i in range(3):
		for j in range(3):
			dctt[0,i] = dcp[j,i]*pt[j] + dctt[0,i]
			dctt[1,i] = dcp[j,i]*b[j] + dctt[1,i]
			dctt[2,i] = dcp[j,i]*s[j] + dctt[2,i]
	
	# compute trend and plunge of normal traction,
	# zero shear, and max. shear traction on plane
	for i in range(3):
		tt[i,1],tt[i,2] = cart_to_sph(dctt[i,0],
			dctt[i,1],dctt[i,2])
	
	return tt, dctt, srat