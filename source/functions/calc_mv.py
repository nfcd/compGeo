import math
from sph_to_cart import sph_to_cart
from cart_to_sph import cart_to_sph

def calc_mv(T,P):
	"""
	calc_mv calculates the mean vector for a group of lines.
	It calculates the trend (trd) and plunge (plg) of the
	mean vector, its mean resultant length (rave), and
	Fisher statistics (concentration factor (conc), 99 (d99)
	and 95 (d95) % uncertainty cones); for a series of lines
	whose trends and plunges are stored in the arrays T and P
	
	NOTE: Input/Output values are in radians
	
	Python function translated from the Matlab function
	CalcMV in Allmendinger et al. (2012)
	"""
	# number of lines
	nlines = len(T)

	# initialize the 3 direction cosines which contain the
	# sums of the individual vectors 
	cn_sum = ce_sum = cd_sum = 0.0
	
	# now add up all the individual vectors
	for i in range(nlines):
		cn,ce,cd = sph_to_cart(T[i],P[i])
		cn_sum += cn
		ce_sum += ce
		cd_sum += cd
	
	# r is the length of the resultant vector and
	# rave is the mean resultant length
	r = math.sqrt(cn_sum*cn_sum+ce_sum*ce_sum+cd_sum*cd_sum)
	rave = r/nlines
	# if rave is lower than 0.1, the mean vector is
	# insignificant, return error
	if rave < 0.1:
		raise ValueError("Mean vector is insignificant")
	#Else 
	else:
		# divide the resultant vector by its length to get
		# the direction cosines of the unit vector
		cn_sum /= r
		ce_sum /= r
		cd_sum /= r
		# convert the mean vector to the lower hemisphere
		if cd_sum < 0.0:
			cn_sum *= -1.0
			ce_sum *= -1.0
			cd_sum *= -1.0
		# convert the mean vector to trend and plunge
		trd, plg = cart_to_sph(cn_sum,ce_sum,cd_sum)
		# if there are enough measurements calculate the
		# Fisher statistics (Fisher et al., 1987)
		conc = d99 = d95 = 0.0
		if r < nlines:
			if nlines < 16:
				afact = 1.0-(1.0/nlines)
				conc = (nlines/(nlines-r))*afact**2
			else:
				conc = (nlines-1.0)/(nlines-r)
		if rave >= 0.65 and rave < 1.0:
			afact = 1.0/0.01
			bfact = 1.0/(nlines-1.0)
			d99 = math.acos(1.0-((nlines-r)/r)*(afact**bfact-1.0))
			afact = 1.0/0.05
			d95 = math.acos(1.0-((nlines-r)/r)*(afact**bfact-1.0))
	
	return trd, plg, rave, conc, d99, d95