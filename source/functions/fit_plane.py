import numpy as np

from pole import plane_from_pole
from cart_to_sph import cart_to_sph

def fit_plane(pts):
	"""
	fit_plane computes the best-fit plane for a group of
	points (position vectors) on the plane
	
	USE: stk, dip, stdev = fit_plane(pts)
	
	pts is a n x 3 matrix containing the East (column 1),
	North (column 2), and Up (column 3) coordinates
	of n points on the plane
	
	stk and dip are returned in radians
	
	stdev is the standard deviation of the distance of
	each point from the best-fit plane
	"""
	# compute the centroid of the selected points
	avge = np.mean(pts[:,0])
	avgn = np.mean(pts[:,1])
	avgu = np.mean(pts[:,2])
	
	# compute the points vectors minus the centroid
	pts[:,0] = pts[:,0] - avge
	pts[:,1] = pts[:,1] - avgn
	pts[:,2] = pts[:,2] - avgu
	
	# compute the covariance/orientation matrix
	a = np.zeros((3,3))
	for i in range(pts.shape[0]):
		ce = pts[i,0]
		cn = pts[i,1]
		cu = pts[i,2]
		# compute orientation matrix
		a[0,0] = a[0,0] + ce*ce
		a[0,1] = a[0,1] + ce*cn
		a[0,2] = a[0,2] + ce*cu
		a[1,1] = a[1,1] + cn*cn
		a[1,2] = a[1,2] + cn*cu
		a[2,2] = a[2,2] + cu*cu
	# the orientation matrix is symmetric so the 
	# off-diagonal components are equal
	a[1,0] = a[0,1]
	a[2,0] = a[0,2]
	a[2,1] = a[1,2]
	
	# calculate the eigenvalues and eigenvectors of the
	# orientation matrix: use function eigh
	D, V = np.linalg.eigh(a)
	
	# calculate pole to best-fit plane = lowest eigenvalue
	# vector in N, E, D coordinates
	cn = V[1,0]
	ce = V[0,0]
	cd = -V[2,0]
	
	# find trend and plunge of pole to best fit plane
	trd, plg =cart_to_sph(cn,ce,cd)
	
	# find Best fit plane
	stk, dip = plane_from_pole(trd,plg)
	
	# calculate standard deviation = square root of
	# minimum eigenvalue
	stdev = np.sqrt(D[0])
	
	return stk, dip, stdev