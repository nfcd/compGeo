import numpy as np
import matplotlib.pyplot as plt
from principal_stress import principal_stress
from sph_to_cart import sph_to_cart
from cart_to_sph import cart_to_sph
from pole import pole_from_plane

def mohr_circle_stress(stress,tx1,px1,tx3,planes,ax):
	"""
	Given the stress tensor in a X1X2X3 coordinate system,
	and a group of n planes, mohr_circle_stress draws the Mohr
	Circle for stress (including the planes). It
	also returns the normal and max. shear tractions on the
	planes and their orientations
	
	ns,ons = mohr_circle_stress(stress,tx1,px1,tx3,planes,ax)
	
	stress = 3 x 3 stress tensor
	tx1 = trend of X1
	px1 = plunge of X1
	tx3 = trend of X3
	planes =  n x 2 vector with strike and dip of planes
	ax = handle to axes where the Mohr Circle will be drawn
	ns = n x 2 vector with the normal and max. shear
		tractions on the planes
	ons = n x 4 vector with the trend and plunge of the
		normal traction (columns 1 and 2), and the
		max. shear traction (columns 3 and 4)
	
	NOTE = Planes orientation follows the right hand rule
		Input and output angles are in radians
	"""
	# tolerance for near zero values
	tol = 1e-6
	
	# compute principal stresses and their orientations
	pstress, dcp = principal_stress(stress,tx1,px1,tx3)
	
	# update stress tensor to principal stresses
	stress = np.zeros((3,3))
	for i in range(3):
		stress[i,i] = pstress[i,0]
	
	# draw sigma1-sigma3 circle
	c = (stress[0,0] + stress[2,2])/2.0
	r = (stress[0,0] - stress[2,2])/2.0
	th =np.arange(0,2*np.pi,np.pi/50)
	costh = np.cos(th)
	sinth = np.sin(th)
	x = r * costh + c
	y = r * sinth
	ax.plot(x,y,"k-")
	
	# draw sigma1-sigma2 circle
	c = (stress[0,0] + stress[1,1])/2.0
	r = (stress[0,0] - stress[1,1])/2.0
	x = r * costh + c
	y = r * sinth
	ax.plot(x,y,"k-")
	
	# draw sigma2-sigma3 circle
	c = (stress[1,1] + stress[2,2])/2.0
	r = (stress[1,1] - stress[2,2])/2.0
	x = r * costh + c
	y = r * sinth
	ax.plot(x,y,"k-")
	
	
	# initialize pole to plane
	p = np.zeros(3)
	
	# initialize vectors with normal and
	# max. shear tractions
	ns = np.zeros((planes.shape[0],2))
	ons = np.zeros((planes.shape[0],4))
	
	# compute normal and max. shear tractions
	for i in range(planes.shape[0]):
		
		# calculate direction cosines of pole to plane
		trd, plg = pole_from_plane(planes[i,0],planes[i,1])
		p[0],p[1],p[2] = sph_to_cart(trd, plg)
		
		# trend and plunge of pole = dir. normal traction
		ons[i,0],ons[i,1] = trd, plg
		
		# transform pole to	 principal stress coordinates
		pt = np.zeros(3)
		for j in range(3):
			for k in range(3):
				pt[j] = dcp[j,k]*p[k] + pt[j]
		
		# calculate the tractions in principal stress
		# coordinates
		t = np.zeros(3)
		for j in range(3):
			for k in range(3):
				t[j] = stress[j,k]*pt[k] + t[j]
		
		# find the b and s axes
		b = np.cross(t,pt)
		s = np.cross(pt,b)
		b = b/np.linalg.norm(b)
		s = s/np.linalg.norm(s)
	
		# transformation matrix from principal
		# stress coordinates to plane coordinates
		a = np.zeros((3,3))
		a[0,:] = pt
		a[1,:] = b
		a[2,:] = s
		
		# normal and max. shear tractions
		ns[i,0] = stress[0,0]*a[0,0]*a[0,0] + stress[1,1]\
			*a[0,1]*a[0,1]+ stress[2,2]*a[0,2]*a[0,2]
		ns[i,1] = stress[0,0]*a[0,0]*a[2,0] + stress[1,1]\
			*a[0,1]*a[2,1]+ stress[2,2]*a[0,2]*a[2,2]
		
		# calculate direction cosines of max.
		# shear traction with respect to NED
		ds = np.zeros(3)
		for j in range(3):
			for k in range(3):
				ds[j] = dcp[k,j]*s[k] + ds[j]
		
		# trend and plunge of max. shear traction
		ons[i,2],ons[i,3] = cart_to_sph(ds[0],ds[1],ds[2])
		
		# cross product of pole and max. shear traction
		ps = np.cross(p,ds)
		
		# make clockwise shear traction negative
		if np.abs(ps[2]) < tol: # Dip slip
			if ds[2] > 0.0: # Normal slip
				if pt[0]*pt[2] < 0.0:
					ns[i,1] *= -1.0
			else:	# Reverse slip
				if pt[0]*pt[2] >= 0.0:
					ns[i,1] *= -1.0
		else:	# Oblique slip
			if ps[2] < 0.0:
				ns[i,1] *= -1.0
	
	# plot planes
	ax.plot(ns[:,0],ns[:,1],"ks")
	
	# make axes equal and plot grid
	ax.axis ("equal")
	ax.grid()
	
	# move x-axis to center and y-axis to origin
	ax.spines["bottom"].set_position("center")
	ax.spines["left"].set_position("zero")
	
	# eliminate upper and right axes
	ax.spines["right"].set_color("none")
	ax.spines["top"].set_color("none")
	
	# show ticks in the left and lower axes only
	ax.xaxis.set_ticks_position("bottom")
	ax.yaxis.set_ticks_position("left")
	
	# add labels at end of axes
	ax.set_xlabel(r"$\sigma$",x=1)
	ax.set_ylabel(r"$\tau$",y=1,rotation=0)
	
	return ns, ons