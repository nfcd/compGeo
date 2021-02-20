import numpy as np
from numpy import linalg as la
from sph_to_cart import sph_to_cart
from cart_to_sph import cart_to_sph

def calc_mv(T, P):
	'''
	calc_mv calculates the mean vector for a group of lines

	calc_mv(T,P) calculates the trend (trd) and plunge (plg)
	of the mean vector, its mean resultant length (Rave) for a series of lines
	whose trends and plunges are stored in the arrays T and P

	NOTE: Input/Output trends and plunges

	Python function translated from the Matlab function
	CalcMV in Allmendinger et al. (2012)
	'''
	# Number of lines
	nlines = len(T)

	# Now add up all the individual vectors. Eq. 4.15a
	CN, CE, CD = sph_to_cart(T, P)

	# Initialize the 3 direction cosines which contain the
	# sums of the individual vectors 
	NED_SUM = np.array([ np.sum(CN), np.sum(CE), np.sum(CD) ])

	# R is the length of the resultant vector and
	# Rave is the mean resultant length. Eqs. 4.15b and d
	R = la.norm(NED_SUM)
	Rave = R / nlines
	# If Rave is lower than 0.1, the mean vector is
	# insignificant, return error
	if Rave < 0.1:
		raise ValueError('Mean vector is insignificant')
	else:
		# Divide the resultant vector by its length to get
		# the direction cosines of the unit vector
		NED_SUM = NED_SUM / R
		# Convert the mean vector to the lower hemisphere
		if NED_SUM[2] < 0.0:
			NED_SUM = -NED_SUM
		# Convert the mean vector to trend and plunge
		trd, plg = cart_to_sph(NED_SUM[0], NED_SUM[1], NED_SUM[2])

	return trd, plg, Rave