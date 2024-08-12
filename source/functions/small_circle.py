import numpy as np
from rotate import rotate
from st_coord_line import st_coord_line
from zero_twopi import zero_twopi

def small_circle(trda,plga,cangle,stype):
	"""
	small_circle computes the paths of a small circle defined
	by its axis and cone angle, for an equal angle or equal
	area stereonet of unit radius
	
	trda = trend of axis
	plga = plunge of axis
	cangle = cone angle
	stype = Stereonet type: 0 = equal angle, 1 = equal area
	path1 and path2 = vectors with the x and y coordinates
		of the points in the small circle paths
	np1 and np2 = Number of points in path1 and path2
	
	NOTE: All angles should be in radians
	
	Python function translated from the Matlab function
	SmallCircle in Allmendinger et al. (2012)
	"""
	pi = np.pi
	# find where to start the small circle
	if (plga - cangle) >= 0.0:
		trd = trda
		plg = plga - cangle
	else:
		if plga == pi/2.0:
			plga *= 0.9999
		angle = np.arccos(np.cos(cangle)/np.cos(plga))
		trd = zero_twopi(trda+angle)
		plg = 0.0
	
	# to make the small circle, rotate the starting line
	# 360 degrees in increments of 1 degree
	rot = np.radians(np.arange(0,361,1))
	path1 = np.zeros((rot.shape[0],2))
	path2 = np.zeros((rot.shape[0],2))
	np1 = np2 = 0
	for i in range(rot.shape[0]):
		# rotate line: The line is considered as a vector
		rtrd , rplg = rotate(trda,plga,rot[i],trd,plg,"v")
		# calculate stereonet coordinates and add to path
		# if rotated plunge is positive add to 1st path
		if rplg >= 0.0:
			path1[np1,0] , path1[np1,1] = st_coord_line(rtrd,rplg,
				stype)
			np1 += 1
		# else add to 2nd path
		else:
			path2[np2,0] , path2[np2,1] = st_coord_line(rtrd,rplg,
				stype)
			np2 += 1
	
	return path1, path2, np1, np2