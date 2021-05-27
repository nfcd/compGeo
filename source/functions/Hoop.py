import numpy as np
import matplotlib.pyplot as plt

def Hoop(geom,stress):
	'''
	Hoop computes the hoop and radial stresses 
	around a circular hole. assuming that 
	the circle is on a principal plane,
	smax (s1) is N-S (theta = 90 or 270), 
	and smin (s3) is E-W (theta = 0 or 180)
	Based on Jaeger et al. (2007)

	USE: [shm, srm] = Hoop(geom,stress)

	geom: A 1 x 2 vector with the number of points 
		along the radius, and the number of points 
		around the circle. These values are used
		to construct the grid around the circle
	stress: A 1 x 3 vector with the value of s1, s3,
		and fluid pressure (pf), all in MPa
	shm, srm: maximum hoop and radial stresses and
		their theta orientations
	'''

	# Geometry
	M = geom[0] # points along radius
	N = geom[1] # points around circle
	R1 = 1.0 # radius of hole = 1.0 
	R2 = R1 * 5 # outer radius = 5.0
	nR = np.linspace(R1,R2,M)
	nT = np.linspace(0,2*np.pi,N)
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
	for i in range(0,m):
		for j in range(0,n):
			# Hoop stress. Eq. 9.10
			sh[i,j] = (s1+s3)/2*(1+(R1/R[i,j])**2) \
				-(s3-s1)/2*(1+3*(R1/R[i,j])**4) \
				*np.cos(2*T[i,j]) - pf*(R1/R[i,j])**2
			# Radial stress. Eq. 9.10
			sr[i,j] = (s1+s3)/2*(1-(R1/R[i,j])**2) \
				+(s3-s1)/2*(1-4*(R1/R[i,j])**2+3*(R1/R[i,j])**4) \
				*np.cos(2*T[i,j]) + pf*(R1/R[i,j])**2
			# maximum hoop stress
			if sh[i,j] > shm[0]:
				shm[0] = sh[i,j]
				shm[1] = T[i,j]*180/np.pi
			# maximum radial stress
			if sr[i,j] > srm[0]:
				srm[0] = sr[i,j]
				srm[1] = T[i,j]*180/np.pi

	# Create figure
	fig = plt.figure(figsize=(10,10))
	ax1 = fig.add_subplot(2, 2, 1)
	ax2 = fig.add_subplot(2, 2, 2)
	ax3 = fig.add_subplot(2, 2, 3)
	ax4 = fig.add_subplot(2, 2, 4)
	
	# Plot hoop stress
	ax1.axis('equal')
	ax1.set_frame_on(False)
	ax1.get_xaxis().set_visible(False)
	ax1.get_yaxis().set_visible(False)	
	ax1.plot(X[:,0],Y[:,0],'k',linewidth=1.5)
	ax1.plot(X[:,n-1],Y[:,n-1],'k',linewidth=1.5)
	cbar = ax1.contourf(X, Y, sh, cmap='jet')	 
	cstr = ax1.contour(X,Y,sh, colors = 'k')
	fig.colorbar(cbar, ax=ax1, label='MPa')
	ax1.set_title('Hoop stress',fontweight="bold")
	
	# Plot radial stress
	ax2.axis('equal')
	ax2.set_frame_on(False)
	ax2.get_xaxis().set_visible(False)
	ax2.get_yaxis().set_visible(False)	
	ax2.plot(X[:,0],Y[:,0],'k',linewidth=1.5)
	ax2.plot(X[:,n-1],Y[:,n-1],'k',linewidth=1.5)
	cbar = ax2.contourf(X, Y, sr, cmap='jet')
	cstr = ax2.contour(X,Y,sr, colors = 'k')
	fig.colorbar(cbar, ax=ax2, label='MPa')
	ax2.set_title('Radial stress',fontweight="bold")
	
	# Plot variation of hoop and radial stress along s3
	ax3.plot(R[0,:]/R1,sh[0,:],'r.-', label='Hoop')
	ax3.plot(R[0,:]/R1,sr[0,:],'b.-', label='Radial')
	ax3.grid(b=True)
	ax3.set_xlabel('Normalized radial distance')
	ax3.set_ylabel('Stress (MPa)')
	ax3.legend(loc='upper right')
	ax3.set_title('Stress variation along s3',fontweight="bold")
	
	# Plot variation of hoop and radial stress around circle
	ax4.plot(T[:,0]*180/np.pi,sh[:,0],'r.-', label='Hoop')
	ax4.plot(T[:,0]*180/np.pi,sr[:,0],'b.-', label='Radial')
	ax4.grid(b=True)
	ax4.set_xlim([0, 360])
	ax4.set_xticks([0, 90, 180, 270, 360])
	ax4.set_xlabel('Angle around the hole (deg)')
	ax4.set_ylabel('Stress (MPa)')
	ax4.legend(loc='upper right')
	ax4.set_title('Stress variation around circle',fontweight="bold")

	return shm, srm