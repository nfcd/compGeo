import numpy as np

from DirCosAxes import DirCosAxes as DirCosAxes
from SphToCart import SphToCart as SphToCart

def Cauchy(stress,tX1,pX1,tX3,strike,dip):
	'''
	Given the stress tensor in a X1X2X3 coordinate system,
	Cauchy computes the X1X2X3 tractions on an arbitrarily
	oriented plane 

	USE: T,pT = Cauchy(stress,tX1,pX1,tX3,strike,dip)

	stress = 3 x 3 stress tensor
	tX1 = trend of X1
	pX1 = plunge of X1
	tX3 = trend of X3
	strike = strike of plane
	dip = dip of plane
	T = 1 x 3 vector with X1, X2 and X3 tractions
	pT = 1 x 3 vector with direction cosines of pole
		to plane with respect to X1X2X3

	NOTE = Plane orientation follows the right hand rule
	Input/Output angles are in radians

	Cauchy uses functions DirCosAxes and SphToCart

	Python function translated from the Matlab function
	Cauchy in Allmendinger et al. (2012)
	'''
	# Compute direction cosines of X1X2X3 with respect
	# to NED
	dC = DirCosAxes(tX1,pX1,tX3)

	# Calculate direction cosines of pole to plane
	p = np.zeros(3)
	p[0],p[1],p[2] = SphToCart(strike,dip,1)
	
	# Transform pole to plane to stress coordinates X1X2X3
	# The transformation matrix is the direction cosines of
	# X1X2X3
	pT = np.zeros(3)
	for i in range(3):
		for j in range(3):
			pT[i] = dC[i,j]*p[j] + pT[i]
			
	# Calculate the tractions in stress coordinates X1X2X3
	T = np.zeros(3)
	# Compute tractions using Cauchy's law: Eq. 7.4
	for i in range(3):
		for j in range(3):
			T[i] = stress[i][j]*pT[j] + T[i]
	
	return T, pT