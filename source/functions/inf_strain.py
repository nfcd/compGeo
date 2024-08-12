import numpy as np

from cart_to_sph import cart_to_sph
from zero_twopi import zero_twopi

def inf_strain(e):
	"""
	inf_strain computes infinitesimal strain from an input
	displacement gradient tensor
	
	USE: eps,ome,pstrain,rotc,rot = inf_strain(e)
	
	e = 3 x 3 displacement gradient tensor
	eps = 3 x 3 strain tensor
	ome = 3 x 3 rotation tensor
	pstrain = 3 x 3 matrix with magnitude (column 1), trend
			(column 2) and plunge (column 3) of maximum
			(row 1), intermediate (row 2),and minimum 
			(row 3) principal strains
	rotc = 1 x 3 vector with rotation components
	rot = 1 x 3 vector with rotation magnitude and trend
		and plunge of rotation axis
	
	NOTE: Output trends and plunges of principal strains
		and rotation axes are in radians
	
	Python function translated from the Matlab function
	InfStrain in Allmendinger et al. (2012)
	"""
	# initialize variables
	eps = np.zeros((3,3))
	ome = np.zeros((3,3))
	pstrain = np.zeros((3,3))
	rotc = np.zeros(3)
	rot = np.zeros(3)
	
	# compute strain and rotation tensors
	for i in range(3):
		for j in range(3):
			eps[i,j]=0.5*(e[i,j]+e[j,i])
			ome[i,j]=0.5*(e[i,j]-e[j,i])
	
	# compute principal strains and orientations.
	# Here we use the function eigh. D is a vector
	# of eigenvalues (i.e. principal strains), and V is a
	# full matrix whose columns are the corresponding
	# eigenvectors (i.e. principal strain directions)
	D,V = np.linalg.eigh(eps)
	
	# maximum principal strain
	pstrain[0,0] = D[2]
	pstrain[0,1],pstrain[0,2] = cart_to_sph(V[0,2],V[1,2],V[2,2])
	# intermediate principal strain
	pstrain[1,0] = D[1] 
	pstrain[1,1],pstrain[1,2] = cart_to_sph(V[0,1],V[1,1],V[2,1])
	# minimum principal strain
	pstrain[2,0] = D[0] 
	pstrain[2,1],pstrain[2,2] = cart_to_sph(V[0,0],V[1,0],V[2,0])
	
	# calculate rotation components
	rotc[0]=(ome[1,2]-ome[2,1])*-0.5
	rotc[1]=(-ome[0,2]+ome[2,0])*-0.5
	rotc[2]=(ome[0,1]-ome[1,0])*-0.5
	
	# compute rotation magnitude
	rot[0] = np.sqrt(rotc[0]**2+rotc[1]**2+rotc[2]**2)
	# compute trend and plunge of rotation axis
	rot[1],rot[2] = cart_to_sph(rotc[0]/rot[0],rotc[1]/rot[0],rotc[2]/rot[0])
	# if plunge is negative
	if rot[2] < 0:
		rot[1] = zero_twopi(rot[1]+np.pi)
		rot[2] *= -1
		rot[0] *= -1
	
	return eps, ome, pstrain, rotc, rot