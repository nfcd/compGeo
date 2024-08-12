import numpy as np
import matplotlib.pyplot as plt

def general_shear(pts,st1,gamma,kk,ninc,ax):
	"""
	general_shear computes displacement paths, kinematic
	vorticity numbers and progressive finite strain 
	history, for general shear with a pure shear stretch,
	no area change, and a single shear strain
	
	USE: paths,wk,pfs = 
		general_shear(pts,st1,gamma,kk,ninc,ax)
	
	pts = npoints x 2 matrix with X1 and X3 coord. of points
	st1 = Pure shear stretch parallel to shear zone
	gamma = Engineering shear strain
	kk = An integer that indicates whether the maximum 
		finite stretch is parallel (kk = 0), or 
		perpendicular (kk=1) to the shear direction
	ninc = number of strain increments
	ax = an array of two axis handles for the plots
	paths = displacement paths of points
	wk = Kinematic vorticity number
	pfs = progressive finite strain history. column 1 =
		orientation of maximum stretch with respect to 
		X1, column 2 = maximum stretch magnitude
	
	NOTE: Intermediate principal stretch is 1.0 (Plane
		strain). Output orientations are in radians
		
	Python function translated from the Matlab function
	GeneralShear in Allmendinger et al. (2012)
	"""
	# compute minimum principal stretch and incr. stretches
	st1inc =st1**(1.0/ninc)
	st3 =1.0/st1
	st3inc =st3**(1.0/ninc)
	
	# incremental engineering shear strain
	gammainc = gamma/ninc
	
	# initialize displacement paths
	npts = pts.shape[0] # number of points
	paths = np.zeros((ninc+1,npts,2))
	paths[0,:,:] = pts # initial points of paths
	
	# calculate incremental deformation gradient tensor
	# if max. stretch parallel to shear direction Eq. 8.45
	if kk == 0:
		F=np.zeros((2,2))
		F[0,]=[st1inc, (gammainc*(st1inc-st3inc))/
			(2.0*np.log(st1inc))]
		F[1,]=[0.0, st3inc]
	# if max. stretch perpendicular to shear direction Eq. 8.46
	elif kk == 1:
		F=np.zeros((2,2))
		F[0,]= [st3inc, (gammainc*(st3inc-st1inc))/
							(2.0*np.log(st3inc))]
		F[1,]= [0.0, st1inc]
	
	# compute displacement paths
	for i in range(npts): # for all points
		for j in range(ninc+1): # for all strain increments
			for k in range(2):
				for L in range(2):
					paths[j,i,k] = F[k,L]*paths[j-1,i,L] + paths[j,i,k]
		# plot displacement path of point
		xx = paths[:,i,0]
		yy = paths[:,i,1]
		ax[0].plot(xx,yy,"k.-")
	
	# plot initial and final polygons
	inpol = np.zeros((npts+1,2))
	inpol[0:npts,]=paths[0,0:npts,:]
	inpol[npts,] = inpol[0,]
	ax[0].plot(inpol[:,0],inpol[:,1],"b-")
	finpol = np.zeros((npts+1,2))
	finpol[0:npts,]=paths[ninc,0:npts,:]
	finpol[npts,] = finpol[0,]
	ax[0].plot(finpol[:,0],finpol[:,1],"r-")
	
	# set axes
	ax[0].set_xlabel(r"$\mathbf{X_1}$")
	ax[0].set_ylabel(r"$\mathbf{X_3}$")
	ax[0].grid()
	ax[0].axis("equal")
	
	# determine the eigenvectors of the flow (apophyses)
	# since F is not symmetrical, use function eig
	_,V = np.linalg.eig(F)
	theta2 = np.arctan(V[1,1]/V[0,1])
	wk = np.cos(theta2)
	
	# initalize progressive finite strain history. 
	# we are not including the initial state
	pfs = np.zeros((ninc,ninc))
	
	# calculate progressive finite strain history
	for i in range(1,ninc+1):
		# determine the finite deformation gradient tensor
		finF = np.linalg.matrix_power(F, i)
		# determine Green deformation tensor
		G = np.dot(finF,finF.conj().transpose())
		# stretch magnitude and orientation: Maximum 
		# eigenvalue and their corresponding eigenvectors
		# of Green deformation tensor
		D, V = np.linalg.eigh(G)
		pfs[i-1,0] = np.arctan(V[1,1]/V[0,1])
		pfs[i-1,1] = np.sqrt(D[1])
	
	# plot progressive finite strain history
	ax[1].plot(pfs[:,0]*180/np.pi,pfs[:,1],"k.-")
	ax[1].set_xlabel(r"$\Theta\;(\circ)$")
	ax[1].set_ylabel("Maximum finite stretch")
	ax[1].set_xlim(-90,90)
	ax[1].set_ylim(1,max(pfs[:,1])+0.5)
	ax[1].grid()
	
	return paths, wk, pfs