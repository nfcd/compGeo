import numpy as np
from pole import pole_from_plane
from rotate import rotate
from st_coord_line import st_coord_line

def great_circle(strike,dip,stype):
	"""
	great_circle computes the great circle path of a plane
	in an equal angle or equal area stereonet of unit radius
	
	strike = strike of plane
	dip = dip of plane
	stype = Stereonet type: 0 = equal angle, 1 = equal area
	path = x and y coordinates of points in great circle path
	
	NOTE: strike and dip should be entered in radians.
	
	Python function translated from the Matlab function
	great_circle in Allmendinger et al. (2012)
	"""
	pi = np.pi
	# Compute the pole to the plane. This will be the axis of
	# rotation to make the great circle
	trda, plga = pole_from_plane(strike,dip)
	
	# Now pick the strike line at the intersection of the
	# great circle with the primitive of the stereonet
	trd, plg = strike, 0.0
	
	# To make the great circle, rotate the line 180 degrees
	# in increments of 1 degree
	rot = np.radians(np.arange(0,181,1))
	path = np.zeros((rot.shape[0],2))
	
	for i in range(rot.size):
	# Avoid joining ends of path
		if rot[i] == pi:
			rot[i] = rot[i] * 0.9999
		# Rotate line
		rtrd, rplg = rotate(trda,plga,rot[i],trd,plg,"a")
		# Calculate stereonet coordinates of rotated line
		# and add to great circle path
		path[i,0], path[i,1] = st_coord_line(rtrd,rplg,stype)
	
	return path