import math
from zero_twopi import zero_twopi

def st_coord_line(trd,plg,stype):
	"""
	st_coord_line computes the coordinates of a line
	in an equal angle or equal area stereonet of unit radius
	
	trd = trend of line
	plg = plunge of line
	stype = Stereonet type: 0 = equal angle, 1 = equal area
	xp and yp: Coordinates of the line in the stereonet
	
	NOTE: trend and plunge should be entered in radians
	
	Python function translated from the Matlab function
	StCoordLine in Allmendinger et al. (2012)
	"""
	# Take care of negative plunges
	if plg < 0:
		trd = zero_twopi(trd+math.pi)
		plg = -plg
	
	# Some constants
	pis4 = math.pi/4
	s2 = math.sqrt(2)
	plgs2 = plg/2
	
	# Equal angle stereonet
	if stype == 0:
		xp = math.tan(pis4 - plgs2)*math.sin(trd)
		yp = math.tan(pis4 - plgs2)*math.cos(trd)
	# Equal area stereonet
	elif stype == 1:
		xp = s2*math.sin(pis4 - plgs2)*math.sin(trd)
		yp = s2*math.sin(pis4 - plgs2)*math.cos(trd)
	
	return xp, yp