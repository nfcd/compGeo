import numpy as np

def down_plunge(bs,trd,plg):
	"""
	down_plunge constructs the down plunge projection of a bed
	
	bs is a npoints x 3 array, which holds npoints
	on the digitized bed, each point defined by
	3 coordinates: X1 = East, X2 = North, X3 = Up
	
	trd and plg are the trend and plunge of the fold axis
	and they should be entered in radians
	
	dpbs are the bed's transformed coordinates
	
	Python function translated from the Matlab function
	DownPlunge in Allmendinger et al. (2012)
	"""
	# Number of points in bed
	nvtex = bs.shape[0]
	
	# Allocate some arrays
	a=np.zeros((3,3))
	dpbs = np.zeros((np.shape(bs)))
	
	# Calculate the transformation matrix a(i,j)
	a[0,0] = np.sin(trd)*np.sin(plg)
	a[0,1] = np.cos(trd)*np.sin(plg)
	a[0,2] = np.cos(plg)
	a[1,0] = np.cos(trd)
	a[1,1] = -np.sin(trd)
	a[2,0] = np.sin(trd)*np.cos(plg)
	a[2,1] = np.cos(trd)*np.cos(plg)
	a[2,2] = -np.sin(plg)
	
	# Perform transformation
	for nv in range(0,nvtex):
		for i in range(0,3):
			dpbs[nv,i] = 0.0
			for j in range(0,3):
					dpbs[nv,i] = a[i,j]*bs[nv,j] + dpbs[nv,i]
	
	return dpbs