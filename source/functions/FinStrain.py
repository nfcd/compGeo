import numpy as np

from CartToSph import CartToSph

def FinStrain(e,frame):
	'''
	FinStrain computes finite strain from an input
	displacement gradient tensor
	
	USE: eps,pstrain,dilat,maxsh = FinStrain(e,frame)
	
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
	
	FinStrain uses function CartToSph
	
	Python function translated from the Matlab function
	FinStrain in Allmendinger et al. (2012)
	'''
	# Initialize variables
	eps = np.zeros((3,3))
	pstrain = np.zeros((3,3))
	maxsh = np.zeros(2)

	# Compute strain tensor
	for i in range(0,3):
		for j in range(0,3):
			eps[i,j] = 0.5*(e[i][j]+e[j][i])
			for k in range(0,3):
				# If undeformed reference frame: 
				# Lagrangian strain tensor, Eq. 8.29
				if frame == 0:
					eps[i,j] = eps[i,j] + 0.5*(e[k][i]*e[k][j])
				# If deformed reference frame: 
				# Eulerian strain tensor, Eq. 8.30
				elif frame == 1:
					eps[i,j] = eps[i,j] - 0.5*(e[k][i]*e[k][j])

	# Compute principal elongations and orientations
	D, V = np.linalg.eigh(eps)

	# Principal elongations
	for i in range(0,3):
		ind = 2-i
		# Magnitude
		# If undeformed reference frame: 
		# Lagrangian strain tensor, Eq. 8.33
		if frame == 0:
			pstrain[i,0] = np.sqrt(1.0+2.0*D[ind])-1.0
		# If deformed reference frame:
		# Eulerian strain tensor, Eq. 8.34
		elif frame == 1:
			pstrain[i,0] = np.sqrt(1.0/(1.0-2.0*D[ind]))-1.0
		# Orientations
		pstrain[i,1],pstrain[i,2] = CartToSph(V[0,ind],V[1,ind],V[2,ind])

	# Dilatation
	dilat = (1.0+pstrain[0,0])*(1.0+pstrain[1,0])*(1.0+pstrain[2,0]) - 1.0

	# Maximum shear strain: This only works if plane strain
	lmax = (1.0+pstrain[0,0])**2 # Maximum quad. elongation
	lmin = (1.0+pstrain[2,0])**2 # Minimum quad. elongation
	# Maximum shear strain: Ramsay (1967) Eq. 3.46
	maxsh[0] = (lmax-lmin)/(2.0*np.sqrt(lmax*lmin))
	# Angle of maximum shear strain with respect to maximum
	# principal strain. Ramsay (1967) Eq. 3.45
	# If undeformed reference frame
	if frame == 0:
		maxsh[1] = np.pi/4.0
	# If deformed reference frame
	elif frame == 1:
		maxsh[1] = np.arctan(np.sqrt(lmin/lmax))

	return eps, pstrain, dilat, maxsh