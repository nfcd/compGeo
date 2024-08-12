import numpy as np

from cart_to_sph import cart_to_sph

def fin_strain(e,frame):
	"""
	fin_strain computes finite strain from an input
	displacement gradient tensor
	
	USE: eps,pstrain,dilat,maxsh = fin_strain(e,frame)
	
	e = 3 x 3 Lagrangian or Eulerian displacement gradient
		tensor
	frame = Reference frame. 0 = undeformed (Lagrangian)
		state, 1 = deformed (Eulerian) state
	eps = 3 x 3 Lagrangian or Eulerian strain tensor
	pstrain = 3 x 3 matrix with magnitude (column 1), trend
		(column 2) and plunge (column 3) of maximum
		(row 1), intermediate (row 2), and minimum
		(row 3) elongations
	dilat = dilatation
	maxsh = 1 x 2 vector with max. shear strain and
		orientation with respect to maximum principal
		strain direction. Only valid in 2D
	
	NOTE: Output angles are in radians
	
	Python function translated from the Matlab function
	FinStrain in Allmendinger et al. (2012)
	"""
	# initialize variables
	eps = np.zeros((3,3))
	pstrain = np.zeros((3,3))
	maxsh = np.zeros(2)
	
	# compute strain tensor
	for i in range(3):
		for j in range(3):
			eps[i,j] = 0.5*(e[i][j]+e[j][i])
			for k in range(3):
				# if undeformed reference frame: 
				# Lagrangian strain tensor
				if frame == 0:
					eps[i,j] = eps[i,j] + 0.5*(e[k][i]*e[k][j])
				# if deformed reference frame: 
				# Eulerian strain tensor
				elif frame == 1:
					eps[i,j] = eps[i,j] - 0.5*(e[k][i]*e[k][j])
	
	# compute principal elongations and orientations
	D, V = np.linalg.eigh(eps)
	
	# Principal elongations
	for i in range(3):
		ind = 2-i
		# magnitude
		# if undeformed reference frame: 
		# Lagrangian strain tensor
		if frame == 0:
			pstrain[i,0] = np.sqrt(1.0+2.0*D[ind])-1.0
		# if deformed reference frame:
		# Eulerian strain tensor
		elif frame == 1:
			pstrain[i,0] = np.sqrt(1.0/(1.0-2.0*D[ind]))-1.0
		# orientations
		pstrain[i,1],pstrain[i,2] = cart_to_sph(V[0,ind],
			V[1,ind],V[2,ind])
	
	# dilatation
	dilat = (1.0+pstrain[0,0])*(1.0+pstrain[1,0])* \
		(1.0+pstrain[2,0]) - 1.0
	
	# maximum shear strain: This only works if plane strain
	lmax = (1.0+pstrain[0,0])**2 # maximum quad. elongation
	lmin = (1.0+pstrain[2,0])**2 # minimum quad. elongation
	# maximum shear strain: Ramsay (1967) Eq. 3.46
	maxsh[0] = (lmax-lmin)/(2.0*np.sqrt(lmax*lmin))
	# angle of maximum shear strain with respect to maximum
	# principal strain. Ramsay (1967) Eq. 3.45
	# if undeformed reference frame
	if frame == 0:
		maxsh[1] = np.pi/4.0
	# if deformed reference frame
	elif frame == 1:
		maxsh[1] = np.arctan(np.sqrt(lmin/lmax))
	
	return eps, pstrain, dilat, maxsh