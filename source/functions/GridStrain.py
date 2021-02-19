import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from scipy.spatial import Delaunay

from lscov import lscov as lscov
from InfStrain import InfStrain

def GridStrain(pos,disp,k,par,plotpar,plotst):
	'''
	GridStrain computes the infinitesimal strain of a network
	of stations with displacements in x (east) and y (north).
	Strain in z is assumed to be zero (plane strain)

	USE: cent,eps,ome,pstrain,rotc = GridStrain(pos,disp,k,par,plotpar,plotst)

	pos = nstations x 2 matrix with x (east) and y (north)
		positions of stations in meters
	disp = nstations x 2 matrix with x (east) and y (north)
		displacements of stations in meters
	k = Type of computation: Delaunay (k = 0), nearest
		neighbor (k = 1), or distance weighted (k = 2)
	par = Parameters for nearest neighbor or distance
		weighted computation. If Delaunay (k = 0), enter
		a scalar corresponding to the minimum internal
		angle of a triangle valid for computation.
		If nearest neighbor (k = 1), input a 1 x 3 vector
		with grid spacing, number of nearest neighbors,
		and maximum distance to neighbors. If distance
		weighted (k = 2), input a 1 x 2 vector with grid
		spacing and distance weighting factor alpha
	plotpar = Parameter to color the cells: Max elongation
		(plotpar = 0), minimum elongation
		(plotpar = 1), rotation (plotpar = 2),
		or dilatation (plotpar = 3)
	plotst = A flag to plot the stations (1) or not (0)
	cent = ncells x 2 matrix with x and y positions of cells
		centroids
	eps = 3 x 3 x ncells array	with strain tensors of
		the cells
	ome = 3 x 3 x ncells array with rotation tensors of
		the cells
	pstrain = 3 x 3 x ncells array with magnitude and
		orientation of principal strains of
		the cells
	rotc = ncells x 3 matrix with rotation components
		of cells

	NOTE: Input/Output angles are in radians. Output
		azimuths are given with respect to North
		pos, disp, grid spacing, max. distance to
		neighbors, and alpha should be in meters
	
	GridStrain uses functions lscov and InfStrain

	Python function translated from the Matlab function
	GridStrain in Allmendinger et al. (2012)
	'''
	pi = np.pi
	# If Delaunay
	if k == 0:
		# Indexes of triangles vertices
		# Use function Delaunay
		tri = Delaunay(pos)
		inds = tri.simplices
		# Number of cells
		ncells = np.size(inds,0)
		# Number of stations per cell = 3
		nstat = 3
		# Centers of cells
		cent = np.zeros((ncells,2))
		for i in range(0,ncells):
			# Triangle vertices
			v1x=pos[inds[i,0],0]
			v2x=pos[inds[i,1],0]
			v3x=pos[inds[i,2],0]
			v1y=pos[inds[i,0],1]
			v2y=pos[inds[i,1],1]
			v3y=pos[inds[i,2],1]
			# Center of cell
			cent[i,0]=(v1x + v2x + v3x)/3.0
			cent[i,1]=(v1y + v2y + v3y)/3.0
			# Triangle internal angles
			s1 = np.sqrt((v3x-v2x)**2 + (v3y-v2y)**2)
			s2 = np.sqrt((v1x-v3x)**2 + (v1y-v3y)**2)
			s3 = np.sqrt((v2x-v1x)**2 + (v2y-v1y)**2)
			a1 = np.arccos((v2x-v1x)*(v3x-v1x)/(s3*s2)+\
						(v2y-v1y)*(v3y-v1y)/(s3*s2))
			a2 = np.arccos((v3x-v2x)*(v1x-v2x)/(s1*s3)+\
						(v3y-v2y)*(v1y-v2y)/(s1*s3))
			a3 = np.arccos((v2x-v3x)*(v1x-v3x)/(s1*s2)+\
						(v2y-v3y)*(v1y-v3y)/(s1*s2))
			# If any of the internal angles is less than
			# specified minimum, invalidate triangle
			if a1 < par or a2 < par or a3 < par:
				inds[i,:] = np.zeros(3)
	# If nearest neighbor or distance weighted
	else:
		# Construct grid
		xmin = min(pos[:,0]); xmax = max(pos[:,0])
		ymin = min(pos[:,1]); ymax = max(pos[:,1])
		cellsx = int(np.ceil((xmax-xmin)/par[0]))
		cellsy = int(np.ceil((ymax-ymin)/par[0]))
		xgrid = np.arange(xmin,(xmin+(cellsx+1)*par[0]),par[0])
		ygrid = np.arange(ymin,(ymin+(cellsy+1)*par[0]),par[0])
		XX,YY = np.meshgrid(xgrid,ygrid)
		# Number of cells
		ncells = cellsx * cellsy
		# Number of stations per cell (nstat) and
		# other parameters
		# If nearest neighbor
		if k == 1:
			nstat = par[1] # max neighbors
			sqmd = par[2]**2 # max squared distance
		# If distance weighted
		elif k == 2:
			nstat = np.size(pos,0) # all stations
			dalpha = 2.0*par[1]*par[1] # 2*alpha*alpha
		# Cells' centers
		cent = np.zeros((ncells,2))
		count = 0
		for i in range(0,cellsy):
			for j in range(0,cellsx):
				cent[count,0] = (XX[i,j]+XX[i,j+1])/2.0
				cent[count,1] = (YY[i,j]+YY[i+1,j])/2.0
				count += 1
		# Initialize stations indexes for cells to -1
		inds = np.ones((ncells,nstat), dtype=int)*-1
		# Initialize weight matrix for distance weighted
		wv = np.zeros((ncells,nstat*2))
		# For all cells set stations indexes
		for i in range(0,ncells):
			# Initialize sq distances to -1.0
			sds = np.ones(nstat)*-1.0
			# For all stations
			for j in range(0,np.size(pos,0)):
				# Sq distance from cell center to station
				dx = cent[i,0] - pos[j,0]
				dy = cent[i,1] - pos[j,1]
				sd = dx**2 + dy**2
				# If nearest neighbor
				if k == 1:
					# If within the max sq distance
					if sd <= sqmd:
						minsd = min(sds)
						mini = np.argmin(sds)
						# If less than max neighbors
						if minsd == -1.0:
							sds[mini] = sd
							inds[i,mini] = j
						# If max neighbors
						else:
							# If sq distance is less 
							# than neighbors max sq distance
							maxsd = max(sds)
							maxi = np.argmax(sds)
							if sd < maxsd:
								sds[maxi] = sd
								inds[i,maxi] = j
				# If distance weighted
				elif k == 2:
					# All stations indexes
					inds[i,:] = np.arange(nstat)
					# Eq. 8.27: Weight factor
					weight = np.exp(-sd/dalpha)
					wv[i,j*2] = weight
					wv[i,j*2+1] = weight
	
	# Initialize arrays
	y = np.zeros(nstat*2)
	M = np.zeros((nstat*2,6))
	e = np.zeros((3,3)) 
	eps = np.zeros((3,3,ncells))
	ome = np.zeros((3,3,ncells)) 
	pstrain = np.zeros((3,3,ncells))
	rotc = np.zeros((ncells,3)) 
		
	# For each cell
	for i in range(0,ncells):
		# If required minimum number of stations
		if min(inds[i,:]) >= 0:
			# Eq. 8.24: Displacements column vector y
			# and design matrix M. X1 = North, X2 = East
			for j in range(0,nstat):
				ic = inds[i,j]
				y[j*2] = disp[ic,1]
				y[j*2+1] = disp[ic,0]
				M[j*2,:] = [1.,0.,pos[ic,1],pos[ic,0],0.,0.]
				M[j*2+1,:] = [0.,1.,0.,0.,pos[ic,1],pos[ic,0]]
			# Eqs. 8.25-8.26: Find x using function lscov
			# If Delaunay or nearest neighbor
			if k == 0 or k == 1:
				x = lscov(M,y)
			# If distance weighted
			elif k == 2:
				x = lscov(M,y,wv[i,:])
			# Displacement gradient tensor
			for j in range(0,2):
				e[j,0] = x[j*2+2]
				e[j,1] = x[j*2+3]
			# Compute strain
			eps[:,:,i],ome[:,:,i],pstrain[:,:,i],\
				rotc[i,:],_ = InfStrain(e)

	# Variable to plot
	# If maximum principal strain
	if plotpar == 0:
		vp = pstrain[0,0,:]
		lcb = "emax"
	# If minimum principal strain
	elif plotpar == 1:
		vp = pstrain[2,0,:]
		lcb = "emin"
	# If rotation: 
	# For plane strain, rotation = rotc(3)
	elif plotpar == 2:
		vp = rotc[:,2]*180/pi
		lcb = "Rotation (deg)"
	# If dilatation
	elif plotpar == 3:
		vp = pstrain[0,0,:]+pstrain[1,0,:]+pstrain[2,0,:]
		lcb = "dilatation"

	# Make a figure
	fig, ax = plt.subplots()
	fig.set_size_inches(15.0, 7.5)
	
	# Patches and colors for cells
	patches = []
	colors = []
	
	# Fill cells patches and colors
	# If Delaunay
	if k == 0:
		for i in range(0,ncells):
			# If minimum number of stations
				if min(inds[i,:]) >= 0:
					xpyp = [[pos[inds[i,0],0],pos[inds[i,0],1]],\
							[pos[inds[i,1],0],pos[inds[i,1],1]],\
							[pos[inds[i,2],0],pos[inds[i,2],1]]]
					# length in km
					xpyp = np.divide(xpyp,1e3)
					polygon = Polygon(xpyp, True)
					patches.append(polygon)
					colors.append(vp[i])
	# If nearest neighbor or distance weighted
	if k == 1 or k == 2:
		count = 0
		for i in range(0,cellsy):
			for j in range(0,cellsx):
				# If minimum number of stations
				if min(inds[count,:]) >= 0:
					xpyp = [[XX[i,j],YY[i,j]],[XX[i,j+1],YY[i,j+1]],\
						[XX[i+1,j+1],YY[i+1,j+1]],[XX[i+1,j],YY[i+1,j]]]
					# length in km
					xpyp = np.divide(xpyp,1e3)
					polygon = Polygon(xpyp, True)
					patches.append(polygon)
					colors.append(vp[count])
				count += 1

	# Collect cells patches
	pcoll = PatchCollection(patches)
	# Cells colors
	pcoll.set_array(np.array(colors))
	# Color map is blue to red
	pcoll.set_cmap('bwr')
	# Positive values are red, negative are
	# blue and zero is white
	vmin = min(vp) 
	vmax = max(vp)
	norm=mcolors.TwoSlopeNorm(vmin=vmin,vcenter=0.0,vmax=vmax)
	pcoll.set_norm(norm)
	
	# Draw cells
	ax.add_collection(pcoll)
	
	# Plot stations
	if plotst == 1:
		plt.plot(pos[:,0]*1e-3,pos[:,1]*1e-3,'k.',markersize=2) 
	
	# Axes
	plt.axis('equal')
	plt.xlabel('x (km)')
	plt.ylabel('y (km)')
	
	# Color bar with nice ticks
	intv = (vmax-vmin)*0.25
	ticks=[vmin,vmin+intv,vmin+2*intv,vmin+3*intv,vmax]
	lticks = ['{:.2e}'.format(ticks[0]),\
			'{:.2e}'.format(ticks[1]),'{:.2e}'.format(ticks[2]),\
			'{:.2e}'.format(ticks[3]),'{:.2e}'.format(ticks[4])]
	cbar = fig.colorbar(pcoll, label=lcb, ticks=ticks)
	cbar.ax.set_yticklabels(lticks)
	
	# Show plot
	plt.show()
	
	return cent,eps,ome,pstrain,rotc