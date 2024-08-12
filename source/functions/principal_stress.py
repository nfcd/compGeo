import numpy as np

from dircos_axes import dircos_axes
from cart_to_sph import cart_to_sph

def principal_stress(stress,tx1,px1,tx3):
	"""
	Given the stress tensor in a X1X2X3 coordinate system
	principal_stress calculates the principal stresses and
	their orientations (trend and plunge) 
	
	USE: pstress,dcp = principal_stress(stress,tx1,px1,tx3)
	
	stress = Symmetric 3 x 3 stress tensor
	tx1 = trend of X1
	px1 = plunge of X1
	tx3 = trend of X3
	pstress = 3 x 3 matrix containing the magnitude
		(column 1), trend (column 2), and plunge
		(column 3) of the	 maximum (row 1),
		intermediate (row 2), and minimum (row 3)
		principal stresses
	dcp = 3 x 3 matrix with direction cosines of the
		principal stress directions: Max. (row 1),
		Int. (row 2), and Min. (row 3) with respect
		to North-East-Down
	
	NOTE: Input/Output angles are in radians
	
	Python function translated from the Matlab function
	PrincipalStress in Allmendinger et al. (2012)
	"""
	# compute direction cosines of X1X2X3
	dc = dircos_axes(tx1,px1,tx3)
	
	# initialize pstress
	pstress = np.zeros((3,3))
	
	# calculate the eigenvalues and eigenvectors
	# of the stress tensor
	D, V = np.linalg.eigh(stress)
	
	# fill principal stress magnitudes
	pstress[0,0] = D[2] # Maximum principal stress
	pstress[1,0] = D[1] # Interm. principal stress
	pstress[2,0] = D[0] # Minimum principal stress
	
	# the direction cosines of the principal stresses are
	# with respect to the X1X2X3 stress coordinate system,
	# so they need to be transformed to the NED coordinate
	# system
	tv = np.zeros((3,3))
	for i in range(3):
		for j in range(3):
			for k in range(3):
				tv[j,i] = dc[k,j]*V[k,i] + tv[j,i]
				
	
	# initialize dcp
	dcp = np.zeros((3,3))
	
	# direction cosines of principal stresses
	for i in range(3):
		for j in range(3):
			dcp[i,j] = tv[j,2-i]
		# avoid precision issues
		# make sure the principal axes are unit vectors
		dcp[i,:] = dcp[i,:]/np.linalg.norm(dcp[i,:])
	
	# trend and plunge of principal stresses
	for i in range(3):
		pstress[i,1],pstress[i,2] = cart_to_sph(dcp[i,0],
			dcp[i,1],dcp[i,2])
	
	return pstress, dcp