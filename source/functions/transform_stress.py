import numpy as np
from dircos_axes import dircos_axes

def transform_stress(stress,tx1,px1,tx3,ntx1,npx1,ntx3):
	"""
	transform_stress transforms a stress tensor from
	old X1X2X3 to new X1'X2'X3' coordinates
	
	nstress=transform_stress(stress,tx1,px1,tx3,ntx1,npx1,ntx3)
	
	stress = 3 x 3 stress tensor
	tx1 = trend of X1 
	px1 = plunge of X1 
	tx3 = trend of X3
	ntx1 = trend of X1'
	npx1 = plunge of X1'
	ntx3 = trend of X3'
	nstress = 3 x 3 stress tensor in new coordinate system
	
	NOTE: All input angles should be in radians
	
	Python function translated from the Matlab function
	TransformStress in Allmendinger et al. (2012)
	"""
	# direction cosines of axes of old coordinate system
	odc = dircos_axes(tx1,px1,tx3)
	
	# direction cosines of axes of new coordinate system
	ndc = dircos_axes(ntx1,npx1,ntx3)
	
	# transformation matrix between old and new
	# coordinate systems
	a = np.zeros((3,3))
	for i in range(3):
		for j in range(3):
			# Use dot product
			a[i,j] = np.dot(ndc[i,:],odc[j,:])
	
	# transform stress
	nstress = np.zeros((3,3))
	for i in range(3):
		for j in range(3):
			for k in range(3):
				for l in range(3):
					nstress[i,j] = a[i,k] * a[j,l] * stress[k,l] \
						+ nstress[i,j]
	
	return nstress