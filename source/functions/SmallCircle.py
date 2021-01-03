import numpy as np
from Rotate import Rotate as Rotate
from StCoordLine import StCoordLine as StCoordLine
from ZeroTwoPi import ZeroTwoPi as ZeroTwoPi

def SmallCircle(trda,plga,coneAngle,sttype):
	'''
	SmallCircle computes the paths of a small circle defined
	by its axis and cone angle, for an equal angle or equal
	area stereonet of unit radius
	
	trda = trend of axis
	plga = plunge of axis
	coneAngle = cone angle
	sttype = type of stereonet. 0 = equal angle, 1 = equal area
	path1 and path2 = vectors with the x and y coordinates of
				the points in the small circle paths
	np1 and np2 = Number of points in path1 and path2
	
	NOTE: All angles should be in radians

	SmallCircle uses functions ZeroTwoPi, StCoordLine and
	Rotate

	Python function translated from the Matlab function
	SmallCircle in Allmendinger et al. (2012)
	'''
	pi = np.pi
	# Find where to start the small circle
	if (plga - coneAngle) >= 0.0:
		trd = trda
		plg = plga - coneAngle
	else:
		if plga == pi/2.0:
			plga = plga*0.9999
		angle = np.arccos(np.cos(coneAngle)/np.cos(plga))
		trd = ZeroTwoPi(trda+angle)
		plg = 0.0
	
	# To make the small circle, rotate the starting line
	# 360 degrees in increments of 1 degree
	rot = np.arange(0,361,1)*pi/180
	path1 = np.zeros((rot.shape[0],2))
	path2 = np.zeros((rot.shape[0],2))
	np1 = 0
	np2 = 0
	for i in range(rot.shape[0]):
		# Rotate line: Notice that the line is considered as
		# a vector
		rtrd , rplg = Rotate(trda,plga,rot[i],trd,plg,'v')
		# Calculate stereonet coordinates and add to the 
		# right path
		# If plunge of rotated line is positive add to
		# first path
		if rplg >= 0.0:
			path1[np1,0] , path1[np1,1] = StCoordLine(rtrd,rplg,sttype)
			np1 = np1 +1
		# Else add to the second path
		else:
			path2[np2,0] , path2[np2,1] = StCoordLine(rtrd,rplg,sttype)
			np2 = np2 +1
	
	return path1, path2, np1, np2
	

