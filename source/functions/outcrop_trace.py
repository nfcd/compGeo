import numpy as np

def outcrop_trace(strike,dip,p1,XG,YG,ZG):
	"""
	outcrop_trace estimates the outcrop trace of a plane,
	given the strike (strike) and dip (dip) of the plane,
	the ENU coordinates of a point (p1) where the plane
	outcrops, and a DEM of the terrain as a regular grid
	of points with E (XG), N (YG) and U (ZG) coordinates.
	
	After using this function, to draw the outcrop trace
	of the plane, you just need to draw the contour 0 on
	the grid XG,YG,DG
	
	NOTE: strike and dip should be input in radians
		p1 must be an array
		XG and YG arrays should be constructed using
		the numpy function meshgrid
	"""
	# make the transformation matrix from ENU coordinates to
	# SDP coordinates. We just need the third row of this
	# matrix
	a = np.zeros((3,3))
	a[2,0] = -np.cos(strike)*np.sin(dip) 
	a[2,1] = np.sin(strike)*np.sin(dip) 
	a[2,2] = -np.cos(dip);
	
	# Initialize DG
	n, m = XG.shape
	DG = np.zeros((n,m))
	
	# Estimate the P coordinate of the outcrop point p1
	P1 = a[2,0]*p1[0] + a[2,1]*p1[1] + a[2,2]*p1[2]
	
	# Estimate the P coordinate at each point of the DEM
	# grid and subtract P1
	for i in range(0,n):
		for j in range(0,m):
			DG[i,j] = P1 - (a[2,0]*XG[i,j] + a[2,1]*YG[i,j] 
				+ a[2,2]*ZG[i,j])
	
	return DG