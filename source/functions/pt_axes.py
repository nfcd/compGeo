import numpy as np
import matplotlib.pyplot as plt

from cart_to_sph import cart_to_sph
from zero_twopi import zero_twopi
from sph_to_cart import sph_to_cart
from stereonet import stereonet
from great_circle import great_circle
from st_coord_line import st_coord_line
from pole import pole_from_plane

def pt_axes(fault,slip,sense, fpsv,ax):
	"""
	pt_axes computes the P and T axes from the orientation
	of several fault planes and their slip vectors. Results
	are plotted in an equal area stereonet
	
	USE: P,T,senseC = pt_axes(fault,slip,sense,ax)
	
	fault = nfaults x 2 vector with strikes and dips of
		faults
	slip = nfaults x 2 vector with trend and plunge of
		slip vectors
	sense = nfaults x 1 vector with sense of faults
	ax = axis handle for the plot
	fpsv = A flag to tell wether the fault plane and
		slip vector are plotted (1) or not
	P = nfaults x 2 vector with trend and plunge of P axes
	T = nfaults x 2 vector with trend and plunge of T axes
	senseC = nfaults x 1 vector with corrected sense of slip
	
	NOTE: Input/Output angles are in radians
	
	Python function based on the Matlab function
	PTAxes in Allmendinger et al. (2012)
	"""
	pi = np.pi
	
	# initialize some vectors
	p = np.zeros(3)
	u = np.zeros(3)
	eps = np.zeros((3,3))
	P = np.zeros(fault.shape)
	T = np.zeros(fault.shape)
	senseC = sense
	
	# for all faults
	for i in range(fault.shape[0]):
		# Direction cosines of pole to fault and slip vector
		trd, plg = pole_from_plane(fault[i,0],fault[i,1])
		p[0],p[1],p[2] = sph_to_cart(trd, plg)
		u[0],u[1],u[2] = sph_to_cart(slip[i,0],slip[i,1])
		# compute u(i)*p(j) + u(j)*p(i)
		for j in range(3):
			for k in range(3):
				eps[j,k]=u[j]*p[k]+u[k]*p[j]
		# compute orientations of principal axes of strain
		# here we use the function eigh
		_,V = np.linalg.eigh(eps)
		# P orientation
		P[i,0],P[i,1] = cart_to_sph(V[0,2],V[1,2],V[2,2])
		if P[i,1] < 0:
			P[i,0] = zero_twopi(P[i,0]+pi)
			P[i,1] *= -1
		# T orientation
		T[i,0],T[i,1] = cart_to_sph(V[0,0],V[1,0],V[2,0]) 
		if T[i,1] < 0.0:
			T[i,0] = zero_twopi(T[i,0]+pi)
			T[i,1] *= -1
		# determine 3rd component of pole cross product slip
		cross = p[0] * u[1] - p[1] * u[0]
		# use cross and first character in sense to
		# determine if kinematic axes should be switched
		s2 = "p"
		if sense[i][0] == "T" or sense[i][0] == "t": 
			s2 = "Y"
		if (sense[i][0]=="R" or sense[i][0]=="r") and cross>0.0:
			s2 = "Y"
		if (sense[i][0]=="L" or sense[i][0]=="l") and cross<0.0: 
			s2 = "Y"
		if s2 == "Y":
			temp1 = P[i,0]
			temp2 = P[i,1]
			P[i,0] = T[i,0]
			P[i,1] = T[i,1]
			T[i,0] = temp1
			T[i,1] = temp2
			if cross < 0.0: 
				senseC[i] = "TL"
			if cross > 0.0:
				senseC[i] = "TR"
		else:
			if cross < 0.0:
				senseC[i] = "NR"
			if cross > 0.0:
				senseC[i] = "NL"
	
	# plot in equal area stereonet
	stereonet(0,90*pi/180,10*pi/180,1,ax)
	# plot P and T axes
	for i in range(fault.shape[0]):
		if fpsv == 1:
			# plot fault
			path = great_circle(fault[i,0],fault[i,1],1)
			ax.plot(path[:,0],path[:,1],"k")
			# plot slip vector (black circle)
			xp,yp = st_coord_line(slip[i,0],slip[i,1],1)
			ax.plot(xp,yp,"ko","MarkerFaceColor","k")
		# plot P axis (blue circle)
		xp,yp = st_coord_line(P[i,0],P[i,1],1)
		ax.plot(xp,yp,"bo","MarkerFaceColor","b")
		# plot T axis (red circle)
		xp,yp = st_coord_line(T[i,0],T[i,1],1)
		ax.plot(xp,yp,"ro","MarkerFaceColor","r")
	
	return P, T, senseC