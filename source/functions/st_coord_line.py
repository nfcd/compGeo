import math
from zero_twopi import zero_twopi

def st_coord_line(trd,plg,stype):
	"""
	st_coord_line computes the coordinates of a line
	on an equal angle or equal area stereonet of unit radius
	
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
		plg *= -1.0
	
	# Equal angle stereonet
	if stype == 0:
		x = math.tan(math.pi/4 - plg/2)
	# Equal area stereonet
	elif stype == 1:
		x = math.sqrt(2) * math.sin(math.pi/4 - plg/2)
	
	# Compute coordinates
	xp = x * math.sin(trd)
	yp = x * math.cos(trd)
	
	return xp, yp