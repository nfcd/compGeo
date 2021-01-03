import numpy as np
import matplotlib.pyplot as plt
from PrincipalStress import PrincipalStress as PrincipalStress
from SphToCart import SphToCart as SphToCart
from CartToSph import CartToSph as CartToSph

def StressMohrCircle(stress,tX1,pX1,tX3,planes):
	'''
	Given the stress tensor in a X1X2X3 coordinate system,
	and a group of n planes, StressMohrCircle draws the Mohr
	Circle for stress (including the planes). It
	also returns the normal and max. shear tractions on the
	planes and their orientations
	
	USE: ns, ons = StressMohrCircle(stress,tX1,pX1,tX3,planes)
	
	stress = 3 x 3 stress tensor
	tX1 = trend of X1
	pX1 = plunge of X1
	tX3 = trend of X3
	planes =  n x 2 vector with strike and dip of planes
	ns = n x 2 vector with the normal and max. shear
		tractions on the planes
	ons = n x 4 vector with the trend and plunge of the
		normal traction (columns 1 and 2), and the
		max. shear traction (columns 3 and 4)

	NOTE = Planes orientation follows the right hand rule
		Input and output angles are in radians

	StressMohrCircle uses functions PrincipalStress
	'''
	# tolerance for near zero values
	tol = 1e-6
	
	# Compute principal stresses and their orientations
	pstress, dCp = PrincipalStress(stress,tX1,pX1,tX3)
	
	# Update stress tensor to principal stresses
	stress = np.zeros((3,3))
	for i in range(3):
		stress[i,i] = pstress[i,0]
	
	# Draw sigma1-sigma3 circle
	c = (stress[0,0] + stress[2,2])/2.0
	r = (stress[0,0] - stress[2,2])/2.0
	th =np.arange(0,2*np.pi,np.pi/50)
	x = r * np.cos(th)+c
	y = r*np.sin(th)
	fig, ax = plt.subplots()
	plt.plot(x,y,'k-')
	
	# Draw sigma1-sigma2 circle
	c = (stress[0,0] + stress[1,1])/2.0
	r = (stress[0,0] - stress[1,1])/2.0
	th =np.arange(0,2*np.pi,np.pi/50)
	x = r * np.cos(th)+c
	y = r*np.sin(th)
	plt.plot(x,y,'k-')
	
	# Draw sigma2-sigma3 circle
	c = (stress[1,1] + stress[2,2])/2.0
	r = (stress[1,1] - stress[2,2])/2.0
	th =np.arange(0,2*np.pi,np.pi/50)
	x = r * np.cos(th)+c
	y = r*np.sin(th)
	plt.plot(x,y,'k-')
	
	
	# Initialize pole to plane
	p = np.zeros(3)
	
	# Initialize vectors with normal and
	# max. shear tractions
	ns = np.zeros((np.size(planes,0),2))
	ons = np.zeros((np.size(planes,0),4))
	
	# Compute normal and max. shear tractions
	for i in range(np.size(planes,0)):
		
		# Calculate direction cosines of pole to plane
		p[0],p[1],p[2] = SphToCart(planes[i,0],planes[i,1],1)
		
		# Trend and plunge of pole = dir. normal traction
		ons[i,0],ons[i,1] = CartToSph(p[0],p[1],p[2])
		
		# Transform pole to	 principal stress coordinates
		pT = np.zeros(3)
		for j in range(3):
			for k in range(3):
				pT[j] = dCp[j,k]*p[k] + pT[j]
		
		# Calculate the tractions in principal stress
		# coordinates
		T = np.zeros(3)
		for j in range(3):
			for k in range(3):
				T[j] = stress[j,k]*pT[k] + T[j]
		
		# Find the B and S axes
		B = np.cross(T,pT);
		S = np.cross(pT,B);
		B = B/np.linalg.norm(B);
		S = S/np.linalg.norm(S);
	
		# Transformation matrix from principal
		# stress coordinates to plane coordinates
		a = np.zeros((3,3))
		a[0,:] = pT
		a[1,:] = B
		a[2,:] = S
		
		# normal and max. shear tractions
		ns[i,0] = stress[0,0]*a[0,0]*a[0,0] + stress[1,1]*a[0,1]*a[0,1]+ stress[2,2]*a[0,2]*a[0,2]
		ns[i,1] = stress[0,0]*a[0,0]*a[2,0] + stress[1,1]*a[0,1]*a[2,1]+ stress[2,2]*a[0,2]*a[2,2]
		
		# Calculate direction cosines of max.
		# shear traction with respect to NED
		ds = np.zeros(3)
		for j in range(3):
			for k in range(3):
				ds[j] = dCp[k,j]*S[k] + ds[j]
		
		# Trend and plunge of max. shear traction
		ons[i,2],ons[i,3] = CartToSph(ds[0],ds[1],ds[2])
		
		# Cross product of pole and max. shear traction
		ps = np.cross(p,ds)
		
		# Make clockwise shear traction negative
		if np.abs(ps[2]) < tol: # Dip slip
			if ds[2] > 0.0: # Normal slip
				if pT[0]*pT[2] < 0.0:
					ns[i,1] *= -1.0
			else:	# Reverse slip
				if pT[0]*pT[2] >= 0.0:
					ns[i,1] *= -1.0
		else:	# Oblique slip
			if ps[2] < 0.0:
				ns[i,1] *= -1.0
				
	# Plot planes
	plt.plot(ns[:,0],ns[:,1],'ks')
	
	# Make axes equal and plot grid
	plt.axis ('equal')
	plt.grid()
	
	# Move x-axis to center and y-axis to origin
	ax.spines['bottom'].set_position('center')
	ax.spines['left'].set_position('zero')
	
	# Eliminate upper and right axes
	ax.spines['right'].set_color('none')
	ax.spines['top'].set_color('none')
	
	# Show ticks in the left and lower axes only
	ax.xaxis.set_ticks_position('bottom')
	ax.yaxis.set_ticks_position('left')
	
	# Add labels at end of axes
	ax.set_xlabel('σ',x=1)
	ax.set_ylabel('τ',y=1,rotation=0)
	
	# show plot
	plt.show()
	
	return ns, ons