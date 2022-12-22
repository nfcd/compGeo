import numpy as np
import matplotlib.pyplot as plt

def hoop(geom,stress):
	"""
	hoop computes the hoop and radial stresses 
	around a circular hole. assuming that 
	the circle is on a principal plane,
	smax (s1) is N-S (theta = 90 or 270), 
	and smin (s3) is E-W (theta = 0 or 180)
	Based on Jaeger et al. (2007)
	
	USE: shm, srm, fig, ax = hoop(geom,stress)
	
	geom: A 1 x 2 vector with the number of points 
		along the radius, and the number of points 
		around the circle. These values are used
		to construct the grid around the circle
	stress: A 1 x 3 vector with the value of s1, s3,
		and fluid pressure (pf), all in MPa
	shm, srm: maximum hoop and radial stresses and
		their theta orientations
		
	fig and ax are handles to the figure and axes
	"""
	pi = np.pi # pi 
	
	# Geometry
	M = geom[0] # points along radius
	N = geom[1] # points around circle
	R1 = 1.0 # radius of hole = 1.0 
	R2 = R1 * 5 # outer radius = 5.0
	nR = np.linspace(R1,R2,M)
	nT = np.linspace(0,2*pi,N)
	R, T = np.meshgrid(nR,nT)
	# Convert grid to cartesian coordintes
	X = R*np.cos(T)
	Y = R*np.sin(T)
	m,n = X.shape
	
	# Principal stresses and pore pressure (MPa)
	s1 = stress[0] 
	s3 = stress[1]
	pf = stress[2]
	
	# Initialize hoop and radial stresses
	sh = np.zeros(X.shape)
	sr = np.zeros(X.shape)
	
	# Initialize maximum hoop and radial stresses
	shm = np.zeros(2)
	srm = np.zeros(2)
	
	# Compute hoop and radial stresses
	for i in range(m):
		for j in range(n):
			# Hoop stress
			sh[i,j] = (s1+s3)/2*(1+(R1/R[i,j])**2) \
				-(s3-s1)/2*(1+3*(R1/R[i,j])**4) \
				*np.cos(2*T[i,j]) - pf*(R1/R[i,j])**2
			# Radial stress
			sr[i,j] = (s1+s3)/2*(1-(R1/R[i,j])**2) \
				+(s3-s1)/2*(1-4*(R1/R[i,j])**2+3*(R1/R[i,j])**4) \
				*np.cos(2*T[i,j]) + pf*(R1/R[i,j])**2
			# maximum hoop stress
			if sh[i,j] > shm[0]:
				shm[0] = sh[i,j]
				shm[1] = T[i,j]*180/pi
			# maximum radial stress
			if sr[i,j] > srm[0]:
				srm[0] = sr[i,j]
				srm[1] = T[i,j]*180/pi
	
	# Create figure
	fig, ax = plt.subplots(2,2,figsize=(10,10))
	
	# Plot hoop stress
	ax[0,0].axis("equal")
	ax[0,0].set_frame_on(False)
	ax[0,0].get_xaxis().set_visible(False)
	ax[0,0].get_yaxis().set_visible(False)	
	ax[0,0].plot(X[:,0],Y[:,0],"k",linewidth=1.5)
	ax[0,0].plot(X[:,n-1],Y[:,n-1],"k",linewidth=1.5)
	cbar = ax[0,0].contourf(X, Y, sh, cmap="jet")	 
	cstr = ax[0,0].contour(X,Y,sh, colors = "k")
	fig.colorbar(cbar, ax=ax[0,0], label="MPa")
	ax[0,0].set_title("Hoop stress",fontweight="bold")
	
	# Plot radial stress
	ax[0,1].axis("equal")
	ax[0,1].set_frame_on(False)
	ax[0,1].get_xaxis().set_visible(False)
	ax[0,1].get_yaxis().set_visible(False)	
	ax[0,1].plot(X[:,0],Y[:,0],"k",linewidth=1.5)
	ax[0,1].plot(X[:,n-1],Y[:,n-1],"k",linewidth=1.5)
	cbar = ax[0,1].contourf(X, Y, sr, cmap="jet")
	cstr = ax[0,1].contour(X,Y,sr, colors = "k")
	fig.colorbar(cbar, ax=ax[0,1], label="MPa")
	ax[0,1].set_title("Radial stress",fontweight="bold")
	
	# Plot variation of hoop and radial stress along s3
	ax[1,0].plot(R[0,:]/R1,sh[0,:],"r.-", label="Hoop")
	ax[1,0].plot(R[0,:]/R1,sr[0,:],"b.-", label="Radial")
	ax[1,0].grid(b=True)
	ax[1,0].set_xlabel("Normalized radial distance")
	ax[1,0].set_ylabel("Stress (MPa)")
	ax[1,0].legend(loc="upper right")
	ax[1,0].set_title("Stress variation along s3",
								fontweight="bold")
	
	# Plot variation of hoop and radial stress around circle
	ax[1,1].plot(T[:,0]*180/pi,sh[:,0],"r.-", label="Hoop")
	ax[1,1].plot(T[:,0]*180/pi,sr[:,0],"b.-", label="Radial")
	ax[1,1].grid(b=True)
	ax[1,1].set_xlim([0, 360])
	ax[1,1].set_xticks([0, 90, 180, 270, 360])
	ax[1,1].set_xlabel("Angle around the hole (deg)")
	ax[1,1].set_ylabel("Stress (MPa)")
	ax[1,1].legend(loc="upper right")
	ax[1,1].set_title("Stress variation around circle",
								fontweight="bold")
	
	return shm, srm, fig, ax