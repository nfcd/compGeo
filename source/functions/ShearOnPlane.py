import numpy as np
from PrincipalStress import PrincipalStress as PrincipalStress
from SphToCart import SphToCart as SphToCart
from CartToSph import CartToSph as CartToSph

def ShearOnPlane(stress,tX1,pX1,tX3,strike,dip):
	'''
	ShearOnPlane calculates the direction and magnitudes of
	the normal and maximum shear tractions on an arbitrarily
	oriented plane

	USE: TT,dCTT,R = ShearOnPlane(stress,tX1,pX1,tX3,strike,dip)

	stress = 3 x 3 stress tensor
	tX1 = trend of X1
	pX1 = plunge of X1
	tX3 = trend of X3
	strike = strike of plane
	dip = dip of plane
	TT = 3 x 3 matrix with the magnitude (column 1),
		trend (column 2) and plunge (column 3) of the:
		normal traction on the plane (row 1), zero shear
		traction (row 2), and max. shear traction (row 3)
	dCTT = 3 x 3 matrix with the direction cosines of unit
		vectors parallel to: normal traction on the plane
		(row 1), zero shear traction (row 2), and maximum
		shear traction (row 3) with respect to NED
		R = Stress ratio
	
	NOTE = Input stress tensor does not need to be along
		principal stress directions
		Plane orientation follows the right hand rule
		Input/Output angles are in radians

	ShearOnPlane uses functions PrincipalStress, Cauchy and
	CartToSph

	Python function translated from the Matlab function
	ShearOnPlane in Allmendinger et al. (2012)
	'''
	# Compute principal stresses and their orientations
	pstress, dCp = PrincipalStress(stress,tX1,pX1,tX3)
	
	# Update stress tensor to principal stress directions
	stress = np.zeros((3,3))
	for i in range(3):
		stress[i,i] = pstress[i,0]
	
	# Calculate stress ratio
	R = (stress[1,1] - stress[0,0])/(stress[2,2]-stress[0,0])
	
	# Calculate direction cosines of pole to plane
	p = np.zeros(3)
	p[0], p[1], p[2] = SphToCart(strike,dip,1);
	
	# Transform pole to plane to  principal stress coordinates
	pT = np.zeros(3);
	for i in range(3):
		for j in range(3):
			pT[i] = dCp[i,j]*p[j] + pT[i]
	
	# Calculate the tractions in principal stress coordinates
	T = np.zeros(3)
	# Compute tractions using Cauchy's law
	for i in range(3):
		for j in range(3):
			T[i] = stress[i,j]*pT[j] + T[i]
	
	# Find the B axis by the cross product of T and pT
	B = np.cross(T,pT)
	
	# Find the max. shear traction orientation
	# by the cross product of pT and B
	S = np.cross(pT,B)
	
	# Convert B and S to unit vectors
	B = B/np.linalg.norm(B)
	S = S/np.linalg.norm(S)
	
	# Now we can write the transformation matrix from
	# principal stress coordinates to plane coordinates
	a = np.zeros((3,3))
	a[0,:] = pT
	a[1,:] = B
	a[2,:] = S
	
	# Transform stress from principal to plane coordinates
	# Do it only for the first row since we are just
	# interested in the plane: sigma'11 = normal traction,
	# sigma'12 = zero, and sigma'13 = max. shear traction
	TT = np.zeros((3,3))
	for i in range(3):
		TT[i,0] = stress[0,0]*a[0,0]*a[i,0] + stress[1,1]*a[0,1]*a[i,1] +stress[2,2]*a[0,2]*a[i,2]
	
	# Transform normal traction, zero shear
	# and max. shear traction to NED coords
	dCTT = np.zeros((3,3))
	for i in range(3):
		for j in range(3):
			dCTT[0,i] = dCp[j,i]*pT[j] + dCTT[0,i]
			dCTT[1,i] = dCp[j,i]*B[j] + dCTT[1,i]
			dCTT[2,i] = dCp[j,i]*S[j] + dCTT[2,i]
	
	# Compute trend and plunge of normal traction,
	# zero shear, and max. shear traction on plane
	for i in range(3):
		TT[i,1],TT[i,2] = CartToSph(dCTT[i,0],dCTT[i,1],dCTT[i,2])
	
	return TT, dCTT, R