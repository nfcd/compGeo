import numpy as np
import matplotlib.pyplot as plt

from sph_to_cart import sph_to_cart
from cart_to_sph import cart_to_sph
from zero_twopi import zero_twopi
from stereonet import stereonet
from great_circle import great_circle
from st_coord_line import st_coord_line


def bingham(T,P,stype,ax):
	"""
	bingham calculates and plots a cylindrical best fit to
	a distribution of poles to bedding. The statistical
	routine is based on algorithms in Fisher et al. (1988)
	
	USE: ev, conf, bf = bingham(T,P,stype,ax)
		
	T and P = Vectors of lines trends and plunges
			respectively
	stype = Stereonet type: 0 = equal angle, 1 = equal area
	ax = axis handle for the plot
		
	ev = 3 x 3 matrix with eigenvalues (column 1), trends
	(column 2) and plunges (column 3) of the eigenvectors.
	Maximum eigenvalue and corresponding eigenvector are
	in row 1, intermediate in row 2, and minimum in row 3.
	
	conf = 2 x 2 matrix with the maximum (column 1) and
	minimum (column 2) radius of the 95% elliptical
	confidence cone around the eigenvector corresponding
	to the largest (row 1), and lowest (row 2) eigenvalue
	
	bf = 1 x 2 vector containing the strike and dip
	(right hand rule) of the best fit great circle to
	the distribution of lines
	
	NOTE: Input/Output trends and plunges, as well as
	confidence cones are in radians. bingham plots the
	input poles, eigenvectors and best fit great circle
	in an equal area stereonet.
	
	Python function translated from the Matlab function
	Bingham in Allmendinger et al. (2012)
	"""
	# some constants
	pi = np.pi
	east = pi/2
	twopi = pi*2
	
	# number of lines
	nlines = len(T)
	
	# initialize the orientation matrix
	a = np.zeros((3,3))
	
	# fill the orientation matrix with the sums of the
	# squares (for the principal diagonal) and the products
	# of the direction cosines of each line. cn, ce and cd
	# are the north, east and down direction cosines
	for i in range(nlines): 
		cn,ce,cd = sph_to_cart(T[i],P[i])
		a[0,0] = a[0,0] + cn*cn
		a[0,1] = a[0,1] + cn*ce
		a[0,2] = a[0,2] + cn*cd
		a[1,1] = a[1,1] + ce*ce
		a[1,2] = a[1,2] + ce*cd
		a[2,2] = a[2,2] + cd*cd
	
	# the orientation matrix is symmetric so the off-diagonal
	# components can be equated
	a[1,0] = a[0,1]
	a[2,0] = a[0,2]
	a[2,1] = a[1,2]
	
	# calculate the eigenvalues and eigenvectors of the
	# orientation matrix using function eigh.
	# D is a vector of eigenvalues and V is a full matrix
	# whose columns are the corresponding eigenvectors
	D, V = np.linalg.eigh(a)
	
	# normalize the eigenvalues by the number of lines and
	# convert the corresponding eigenvectors to the lower
	# hemisphere
	for i in range(3): 
		D[i] = D[i]/nlines
		if V[2,i] < 0:
			V[0,i] = -V[0,i]
			V[1,i] = -V[1,i]
			V[2,i] = -V[2,i]
	
	# initialize ev
	ev = np.zeros((3,3))
	# fill ev
	ev[0,0] = D[2]	  # Maximum eigenvalue
	ev[1,0] = D[1]	  # Intermediate eigenvalue
	ev[2,0] = D[0]	  # Minimum eigenvalue
	# trend and plunge of largest eigenvalue: column 3 of V
	ev[0,1], ev[0,2] = cart_to_sph(V[0,2], V[1,2], 
		V[2,2])
	# trend and plunge of interm. eigenvalue: column 2 of V
	ev[1,1], ev[1,2] = cart_to_sph(V[0,1], V[1,1], 
		V[2,1])
	# trend and plunge of minimum eigenvalue: column 1 of V
	ev[2,1], ev[2,2] = cart_to_sph(V[0,0], V[1,0], 
		V[2,0])
	
	# initialize conf
	conf = np.zeros((2,2))
	# if there are more than 25 lines, calculate confidence
	# cones at the 95% confidence level. The algorithm is
	# explained in Fisher et al. (1987)
	if nlines >= 25:
		e11 = e22 = e12 = d11 = d22 = d12 = 0
		en11 = 1/(nlines*(ev[2,0]-ev[0,0])**2)
		en22 = 1/(nlines*(ev[1,0]-ev[0,0])**2)
		en12 = 1/(nlines*(ev[2,0]-ev[0,0])*(ev[1,0]-
			ev[0,0]))
		dn11 = en11
		dn22 = 1/(nlines*(ev[2,0]-ev[1,0])**2)
		dn12 = 1/(nlines*(ev[2,0]-ev[1,0])*(ev[2,0]-
			ev[0,0]))
		vec = np.zeros((3,3))
		for i in range(3):
			vec[i,0] = np.sin(ev[i,2]+east)*np.cos(twopi-
				ev[i,1])
			vec[i,1] = np.sin(ev[i,2]+east)*np.sin(twopi-
				ev[i,1])
			vec[i,2] = np.cos(ev[i,2]+east)
		for i in range(nlines):
			c1 = np.sin(P[i]+east)*np.cos(twopi-T[i])
			c2 = np.sin(P[i]+east)*np.sin(twopi-T[i])
			c3 = np.cos(P[i]+east)
			u1x = vec[2,0]*c1 + vec[2,1]*c2 + vec[2,2]*c3
			u2x = vec[1,0]*c1 + vec[1,1]*c2 + vec[1,2]*c3
			u3x = vec[0,0]*c1 + vec[0,1]*c2 + vec[0,2]*c3
			e11 = u1x*u1x * u3x*u3x + e11
			e22 = u2x*u2x * u3x*u3x + e22
			e12 = u1x*u2x * u3x*u3x + e12
			d11 = e11
			d22 = u1x*u1x * u2x*u2x + d22
			d12 = u2x*u3x * u1x*u1x + d12
		e22 = en22*e22
		e11 = en11*e11
		e12 = en12*e12
		d22 = dn22*d22
		d11 = dn11*d11
		d12 = dn12*d12
		d = -2*np.log(.05)/nlines
		# initialize f
		f = np.zeros((2,2))
		if abs(e11*e22-e12*e12) >= 0.000001:
			f[0,0] = (1/(e11*e22-e12*e12)) * e22
			f[1,1] = (1/(e11*e22-e12*e12)) * e11
			f[0,1] = -(1/(e11*e22-e12*e12)) * e12
			f[1,0] = f[0,1]
			# calculate the eigenvalues and eigenvectors
			# of the matrix f using function eigh
			# The next lines follow steps 1-4 outlined
			# on pp. 34-35 of Fisher et al. (1987)
			DD, _ = np.linalg.eigh(f)
			if DD[0] > 0 and DD[1] > 0:
				if d/DD[0] <= 1 and d/DD[1] <= 1:
					conf[0,1] = np.arcsin(np.sqrt(d/DD[1]))
					conf[0,0] = np.arcsin(np.sqrt(d/DD[0]))
		# repeat the process for the eigenvector
		# corresponding to the smallest eigenvalue
		if abs(d11*d22-d12*d12) >= 0.000001:
			f[0,0] = (1/(d11*d22-d12*d12)) * d22
			f[1,1] = (1/(d11*d22-d12*d12)) * d11
			f[0,1] = -(1/(d11*d22-d12*d12)) * d12
			f[1,0] = f[0,1]
			DD, _ = np.linalg.eigh(f)
			if DD[0] > 0.0 and DD[1] > 0.0:
				if d/DD[0] <= 1 and d/DD[1] <= 1:
					conf[1,1] = np.arcsin(np.sqrt(d/DD[1]))
					conf[1,0] = np.arcsin(np.sqrt(d/DD[0]))
	
	# calculate the best fit great circle
	# to the distribution of points
	bf = np.zeros(2)
	bf[0] = zero_twopi(ev[2,1] + east)
	bf[1] = east - ev[2,2]
	
	# Plot stereonet
	stereonet(0, 90*pi/180, 10*pi/180, stype, ax)
	
	# Plot lines
	for i in range(nlines):
		xp,yp = st_coord_line(T[i],P[i],stype)
		ax.plot(xp,yp,"k.")
	
	# Plot eigenvectors
	for i in range(3):
		xp,yp = st_coord_line(ev[i,1],ev[i,2],stype)
		ax.plot(xp,yp,"rs")
		ax.text(xp-0.03,yp+0.03,str(i+1),c="r")
	
	# Plot best fit great circle
	path = great_circle(bf[0],bf[1],stype)
	ax.plot(path[:,0],path[:,1],"r")
	
	return ev, conf, bf