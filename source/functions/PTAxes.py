import numpy as np
import matplotlib.pyplot as plt

from CartToSph import CartToSph as CartToSph
from ZeroTwoPi import ZeroTwoPi as ZeroTwoPi
from SphToCart import SphToCart as SphToCart
from Stereonet import Stereonet as Stereonet
from GreatCircle import GreatCircle as GreatCircle
from StCoordLine import StCoordLine as StCoordLine

def PTAxes(fault,slip,sense, fpsv):
	'''
	PTAxes computes the P and T axes from the orientation
	of several fault planes and their slip vectors. Results
	are plotted in an equal area stereonet
	
	USE: P,T,senseC = PTAxes(fault,slip,sense)
	
	fault = nfaults x 2 vector with strikes and dips of
		faults
	slip = nfaults x 2 vector with trend and plunge of
		slip vectors
	sense = nfaults x 1 vector with sense of faults
	fpsv = A flag to tell wether the fault plane and
		slip vector are plotted (1) or not
	P = nfaults x 2 vector with trend and plunge of P axes
	T = nfaults x 2 vector with trend and plunge of T axes
	senseC = nfaults x 1 vector with corrected sense of slip
	
	NOTE: Input/Output angles are in radians
	
	PTAxes uses functions SphToCart, ZeroTwoPi, CartToSph,
	Stereonet, GreatCircle and StCoordLine
	
	Python function based on the Matlab function
	PTAxes in Allmendinger et al. (2012)
	'''
	pi = np.pi
	# Initialize some vectors
	n = np.zeros(3)
	u = np.zeros(3)
	eps = np.zeros((3,3))
	P = np.zeros((np.size(fault,0),2))
	T = np.zeros((np.size(fault,0),2))
	senseC = sense

	# For all faults
	for i in range(np.size(fault,0)):
		# Direction cosines of pole to fault and slip vector
		n[0],n[1],n[2] = SphToCart(fault[i,0],fault[i,1],1)
		u[0],u[1],u[2] = SphToCart(slip[i,0],slip[i,1],0)
		# Compute u(i)*n(j) + u(j)*n(i)
		for j in range(3):
			for k in range(3):
				eps[j,k]=u[j]*n[k]+u[k]*n[j]
		# Compute orientations of principal axes of strain
		# Here we use the function eigh
		_,V = np.linalg.eigh(eps)
		# P orientation
		P[i,0],P[i,1] = CartToSph(V[0,2],V[1,2],V[2,2])
		if P[i,1] < 0:
			P[i,0] = ZeroTwoPi(P[i,0]+pi)
			P[i,1] = P[i,1]*-1.0
		# T orientation
		T[i,0],T[i,1] = CartToSph(V[0,0],V[1,0],V[2,0]) 
		if T[i,1] < 0.0:
			T[i,0] = ZeroTwoPi(T[i,0]+pi)
			T[i,1] = T[i,1]*-1.0
		# Determine 3rd component of pole cross product slip
		cross = n[0] * u[1] - n[1] * u[0]
		# Use cross and first character in sense to
		# determine if kinematic axes should be switched
		s2 = 'N'
		if sense[i][0] == 'T' or sense[i][0] == 't': 
			s2 = 'Y'
		if (sense[i][0]=='R' or sense[i][0]=='r') and cross>0.0:
			s2 = 'Y'
		if (sense[i][0]=='L' or sense[i][0]=='l') and cross<0.0: 
			s2 = 'Y'
		if s2 == 'Y':
			temp1 = P[i,0]
			temp2 = P[i,1]
			P[i,0] = T[i,0]
			P[i,1] = T[i,1]
			T[i,0] = temp1
			T[i,1] = temp2
			if cross < 0.0: 
				senseC[i] = 'TL'
			if cross > 0.0:
				senseC[i] = 'TR'
		else:
			if cross < 0.0:
				senseC[i] = 'NR'
			if cross > 0.0:
				senseC[i] = 'NL'

	# Plot in equal area stereonet
	Stereonet(0,90*pi/180,10*pi/180,1)
	# Plot P and T axes
	for i in range(np.size(fault,0)):
		if fpsv == 1:
			# Plot fault
			path = GreatCircle(fault[i,0],fault[i,1],1)
			plt.plot(path[:,0],path[:,1],'k')
			# Plot slip vector (black circle)
			xp,yp = StCoordLine(slip[i,0],slip[i,1],1)
			plt.plot(xp,yp,'ko','MarkerFaceColor','k')
		# Plot P axis (blue circle)
		xp,yp = StCoordLine(P[i,0],P[i,1],1)
		plt.plot(xp,yp,'bo','MarkerFaceColor','b')
		# Plot T axis (red circle)
		xp,yp = StCoordLine(T[i,0],T[i,1],1)
		plt.plot(xp,yp,'ro','MarkerFaceColor','r')

	# Show plot
	plt.show()

	return P,T,senseC
