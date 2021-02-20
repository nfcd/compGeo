import numpy as np
from zero_two_pi import zero_two_pi

# Some constants
pi_4 = np.pi/4
sqrt_2 = np.sqrt(2)

def eq_angle_stereonet(trd, plg):
	trd, plg = lower_hemisphere(trd, plg)
	plg_2 = plg / 2

	# Equal angle stereonet, Eq. 3.6
	xp = np.tan(pi_4 - plg_2)*np.sin(trd)
	yp = np.tan(pi_4 - plg_2)*np.cos(trd)

	return xp, yp

def eq_area_stereonet(trd, plg):
	trd, plg = lower_hemisphere(trd, plg)
	plg_2 = plg / 2

	# Equal area stereonet, Eq. 3.7
	xp = sqrt_2*np.sin(pi_4 - plg_2)*np.sin(trd)
	yp = sqrt_2*np.sin(pi_4 - plg_2)*np.cos(trd)

	return xp, yp

# Take care of negative plunges
def lower_hemisphere(trd, plg):
	condition = plg < 0
	trd = np.where(condition, zero_two_pi(trd + np.pi), trd)
	plg = np.where(condition, -plg, plg)

	return trd, plg