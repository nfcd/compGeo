import math
from ZeroTwoPi import ZeroTwoPi as ZeroTwoPi

def Pole(trd,plg,k):
	'''
	Pole returns the pole to a plane or the plane from a pole

	If k = 0, Pole returns the strike (trd1) and dip (plg1)
	of a plane, given the trend (trd) and plunge (plg)
	of its pole.

	If k = 1, Pole returns the trend (trd1) and plunge (plg1)
	of a pole, given the strike (trd) and dip (plg)
	of its plane.

	NOTE: Input/Output angles are in radians.
	Input/Output strike and dip follow the RHR format

	Pole uses function ZeroTwoPi
	'''
	# Some constants
	east = math.pi/2

	# Eq. 3.2
	# Calculate plane given its pole
	if k == 0:
		if plg >= 0:
			plg1 = east - plg
			trd1 = ZeroTwoPi(trd + east)
		else: # Unusual case of pole pointing upwards
			plg1 = east + plg
			trd1 = ZeroTwoPi(trd - east)
	# Else calculate pole given its plane
	elif k == 1:
		plg1 = east - plg;
		trd1 = ZeroTwoPi (trd - east)
		
	return trd1, plg1