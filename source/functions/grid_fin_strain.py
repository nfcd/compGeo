import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from matplotlib.patches import Polygon
from matplotlib.collections import PatchCollection
from scipy.spatial import Delaunay
from lscov import lscov
from fin_strain import fin_strain

def grid_fin_strain(pos,disp,frame,k,par,plotpar,plotst,fig,ax):
	"""
	grid_fin_strain computes the finite strain of a group
	of points with displacements in x and y.
	Strain in z is assumed to be zero (plane strain)
	
	USE: cent,eps,pstrain,dilat,maxsh = grid_fin_strain(pos,
			disp,frame,k,par,plotpar,plotst,fig,ax)
	
	pos = npoints x 2 matrix with x and y position
		of points
	disp = nstations x 2 matrix with x and y
		displacements of points
	frame = Reference frame. 0 = undeformed (Lagrangian)
		state, 1 = deformed (Eulerian) state
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
		(plotpar = 1), dilatation (plotpar = 2),
		or max. shear strain (plotpar = 3)
	plotst = A flag to plot the stations (1) or not (0)
	fig = figure handle for the plot
	ax = axis handle for the plot
	cent = ncells x 2 matrix with x and y positions of the
		cells centroids
	eps = 3 x 3 x ncells array	with strain tensors of
		the cells
	pstrain = 3 x 3 x ncells array with magnitude and
		orientation of principal strains of the cells
	dilat = ncells x 1 vector with dilatation of the cells
	maxsh = ncells x 2 matrix with max. shear strain and
		orientation with respect to maximum principal
		strain direction, of the cells.
		Only valid for plane strain
	
	NOTE: Input/Output angles are in radians. Output
		azimuths are given with respect to y
		pos, disp, grid spacing, max. distance to
		neighbors, and alpha should be in the same
		length units
	"""
	# if Delaunay
	if k == 0:
		# indexes of triangles vertices
		# use function Delaunay
		tri = Delaunay(pos)
		inds = tri.simplices
		# number of cells
		ncells = inds.shape[0]
		# number of stations per cell = 3
		nstat = 3
		# centers of cells
		cent = np.zeros((ncells,2))
		for i in range(ncells):
			# triangle vertices
			v1x=pos[inds[i,0],0]
			v2x=pos[inds[i,1],0]
			v3x=pos[inds[i,2],0]
			v1y=pos[inds[i,0],1]
			v2y=pos[inds[i,1],1]
			v3y=pos[inds[i,2],1]
			# center of cell
			cent[i,0]=(v1x + v2x + v3x)/3.0
			cent[i,1]=(v1y + v2y + v3y)/3.0
			# triangle internal angles
			s1 = np.sqrt((v3x-v2x)**2 + (v3y-v2y)**2)
			s2 = np.sqrt((v1x-v3x)**2 + (v1y-v3y)**2)
			s3 = np.sqrt((v2x-v1x)**2 + (v2y-v1y)**2)
			a1 = np.arccos((v2x-v1x)*(v3x-v1x)/(s3*s2)+\
				(v2y-v1y)*(v3y-v1y)/(s3*s2))
			a2 = np.arccos((v3x-v2x)*(v1x-v2x)/(s1*s3)+\
				(v3y-v2y)*(v1y-v2y)/(s1*s3))
			a3 = np.arccos((v2x-v3x)*(v1x-v3x)/(s1*s2)+\
				(v2y-v3y)*(v1y-v3y)/(s1*s2))
			# if any of the internal angles is less than
			# specified minimum, invalidate triangle
			if a1 < par or a2 < par or a3 < par:
				inds[i,:] = np.zeros(3)
	# if nearest neighbor or distance weighted
	else:
		# construct grid
		xmin = min(pos[:,0]); xmax = max(pos[:,0])
		ymin = min(pos[:,1]); ymax = max(pos[:,1])
		cellsx = int(np.ceil((xmax-xmin)/par[0]))
		cellsy = int(np.ceil((ymax-ymin)/par[0]))
		xgrid = np.arange(xmin,(xmin+(cellsx+1)*par[0]),par[0])
		ygrid = np.arange(ymin,(ymin+(cellsy+1)*par[0]),par[0])
		XX,YY = np.meshgrid(xgrid,ygrid)
		# number of cells
		ncells = cellsx * cellsy
		# number of stations per cell (nstat) and
		# other parameters
		# if nearest neighbor
		if k == 1:
			nstat = par[1] # max neighbors
			sqmd = par[2]**2 # max squared distance
		# if distance weighted
		elif k == 2:
			nstat = pos.shape[0] # all stations
			dalpha = 2.0*par[1]*par[1] # 2*alpha*alpha
		# cells' centers
		cent = np.zeros((ncells,2))
		count = 0
		for i in range(cellsy):
			for j in range(cellsx):
				cent[count,0] = (XX[i,j]+XX[i,j+1])/2.0
				cent[count,1] = (YY[i,j]+YY[i+1,j])/2.0
				count += 1
		# initialize stations indexes for cells to -1
		inds = np.ones((ncells,nstat), dtype=int)*-1
		# initialize weight matrix for distance weighted
		wv = np.zeros((ncells,nstat*2))
		# for all cells set stations indexes
		for i in range(ncells):
			# initialize sq distances to -1.0
			sds = np.ones(nstat)*-1.0
			# for all stations
			for j in range(pos.shape[0]):
				# square distance from cell center to station
				dx = cent[i,0] - pos[j,0]
				dy = cent[i,1] - pos[j,1]
				sd = dx**2 + dy**2
				# if nearest neighbor
				if k == 1:
					# if within the max sq distance
					if sd <= sqmd:
						minsd = min(sds)
						mini = np.argmin(sds)
						# if less than max neighbors
						if minsd == -1.0:
							sds[mini] = sd
							inds[i,mini] = j
						# if max neighbors
						else:
							# ff sq distance is less
							# than neighbors max sq distance
							maxsd = max(sds)
							maxi = np.argmax(sds)
							if sd < maxsd:
								sds[maxi] = sd
								inds[i,maxi] = j
				# if distance weighted
				elif k == 2:
					# all stations indexes
					inds[i,:] = np.arange(nstat)
					# Eq. 8.27: Weight factor
					weight = np.exp(-sd/dalpha)
					wv[i,j*2] = weight
					wv[i,j*2+1] = weight
	
	# initialize arrays
	y = np.zeros(nstat*2)
	M = np.zeros((nstat*2,6))
	e = np.zeros((3,3)) 
	eps = np.zeros((3,3,ncells)) 
	pstrain = np.zeros((3,3,ncells))
	dilat = np.zeros((ncells,1))
	maxsh = np.zeros((ncells,2)) 
		
	# for each cell
	for i in range(ncells):
		# if required minimum number of stations
		if min(inds[i,:]) >= 0:
			# Eq. 8.24: Displacements column vector y
			# and design matrix M. X1 = y, X2 = x
			for j in range(nstat):
				ic = inds[i,j]
				y[j*2] = disp[ic,1]
				y[j*2+1] = disp[ic,0]
				M[j*2,:] = [1.,0.,pos[ic,1],pos[ic,0],0.,0.]
				M[j*2+1,:] = [0.,1.,0.,0.,pos[ic,1],pos[ic,0]]
			# Eqs. 8.25-8.26: Find x using function lscov
			# if Delaunay or nearest neighbor
			if k == 0 or k == 1:
				x = lscov(M,y)
			# if distance weighted
			elif k == 2:
				x = lscov(M,y,wv[i,:])
			# displacement gradient tensor
			for j in range(2):
				e[j,0] = x[j*2+2]
				e[j,1] = x[j*2+3]
			# compute strain
			eps[:,:,i],pstrain[:,:,i],dilat[i,:],\
				maxsh[i,:] = fin_strain(e,frame)

	# variable to plot
	# if maximum principal strain
	if plotpar == 0:
		vp = pstrain[0,0,:]
		vmin = 0.0
		vmax = 2.0
		lcb = "emax"
	# if minimum principal strain
	elif plotpar == 1:
		vp = pstrain[2,0,:]
		vmin = -2.0
		vmax = 0.0
		lcb = "emin"
	# if dilatation:
	elif plotpar == 2:
		vp = dilat[:]
		vmin = -1.0
		vmax = 1.0
		lcb = "dilatation"
	# if max. shear strain
	elif plotpar == 3:
		vp = maxsh[:,0]
		vmin = 0.0
		vmax = 2.0
		lcb = "max. shear strain"

	# Patches and colors for cells
	patches = []
	colors = []
	
	# fill cells patches and colors
	# if Delaunay
	if k == 0:
		for i in range(ncells):
			# If minimum number of stations
				if min(inds[i,:]) >= 0:
					xpyp = [[pos[inds[i,0],0],pos[inds[i,0],1]],\
							[pos[inds[i,1],0],pos[inds[i,1],1]],\
							[pos[inds[i,2],0],pos[inds[i,2],1]]]
					polygon = Polygon(xpyp, True)
					patches.append(polygon)
					colors.append(vp[i])
	# if nearest neighbor or distance weighted
	if k == 1 or k == 2:
		count = 0
		for i in range(cellsy):
			for j in range(cellsx):
				# if minimum number of stations
				if min(inds[count,:]) >= 0:
					xpyp = [[XX[i,j],YY[i,j]],[XX[i,j+1],YY[i,j+1]],\
						[XX[i+1,j+1],YY[i+1,j+1]],[XX[i+1,j],YY[i+1,j]]]
					polygon = Polygon(xpyp, True)
					patches.append(polygon)
					colors.append(vp[count])
				count += 1

	# collect cells patches
	pcoll = PatchCollection(patches)
	# cells colors
	pcoll.set_array(np.array(colors))
	# color map is blue to red
	pcoll.set_cmap("jet")
	norm = mcolors.Normalize(vmin=vmin,vmax=vmax)
	pcoll.set_norm(norm)
	
	# draw cells
	ax.add_collection(pcoll)
	
	# plot stations
	if plotst == 1:
		ax.plot(pos[:,0],pos[:,1],"k.",markersize=2) 
	
	# axes
	ax.axis("equal")
	ax.set_xlabel("x")
	ax.set_ylabel("y")
	
	# color bar with nice ticks
	intv = (vmax-vmin)*0.25
	ticks=[vmin,vmin+intv,vmin+2*intv,vmin+3*intv,vmax]
	lticks = ["{:.2}".format(ticks[0]),\
			"{:.2}".format(ticks[1]),"{:.2}".format(ticks[2]),\
			"{:.2}".format(ticks[3]),"{:.2}".format(ticks[4])]
	cbar = fig.colorbar(pcoll, label=lcb, ticks=ticks)
	cbar.ax.set_yticklabels(lticks)
	
	return cent, eps, pstrain, dilat, maxsh