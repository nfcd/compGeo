import numpy as np

from DirCosAxes import DirCosAxes as DirCosAxes

def TransformStress(stress,tX1,pX1,tX3,ntX1,npX1,ntX3):
    '''
	TransformStress transforms a stress tensor from 
	old X1X2X3 to new X1'X2'X3' coordinates

	USE: nstress = TransformStress(stress,tX1,pX1,tX3,ntX1,npX1,ntX3)

	stress = 3 x 3 stress tensor
	tX1 = trend of X1 
	pX1 = plunge of X1 
	tX3 = trend of X3
	ntX1 = trend of X1'
	npX1 = plunge of X1'
	ntX3 = trend of X3'
	nstress = 3 x 3 stress tensor in new coordinate system
	
	NOTE: All input angles should be in radians
	
	TransformStress uses function DirCosAxes
	
	Python function translated from the Matlab function
	TransformStress in Allmendinger et al. (2012)
	'''
    
    # Direction cosines of axes of old coordinate system
    odC = DirCosAxes(tX1,pX1,tX3)

	# Direction cosines of axes of new coordinate system
    ndC = DirCosAxes(ntX1,npX1,ntX3)
	
	# Transformation matrix between old and new 
	# coordinate systems
    a = np.zeros((3,3))
    for i in range(3):
        for j in range(3):
        	# Use dot product
            a[i,j] = np.dot(ndC[i,:],odC[j,:])
	
	# Transform stress: Eq. 7.5
    nstress = np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                for l in range(3):
                    nstress[i,j] = a[i,k] * a[j,l] * stress[k,l] + nstress[i,j]
	
    return nstress