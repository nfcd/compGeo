import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

def flex2d(geom,elas,loads,fig,ax):
	"""
	flex2d computes the deflection profile produced by 
	a group of load columns on an elastic lithosphere 
	resting on a fluid-like foundation 
	(i.e. astenosphere). The program is based on 
	Hetenyi (1946) solution for an infinite elastic 
	beam on a fluid-like foundation
	
	USE: w, wp = flex2d(geom,elas,loads,ax)
	
	geom: A 1 x 2 vector with the lateral extent of the 
		domain in meters, and the distance between
		points (x interval) in meters where the
		deflection will be computed
	elas: A 1 x 4 vector with the Young Modulus (in Pa), 
		Poisson ratio, Elastic thickness in meters, 
		and density of the foundation in kg/m^3
	loads: A nloads x 4 vector with the left x coordinate,
		right x coordinate, height, and density of each
		load column. Lengths should be in meters and 
		density in kg/m^3
	fig: a figure handle for the plots
	ax: an array of two axis handles for the plots
	w: The deflection of the lithosphere in meters at the 
		points specified by the extent and x interval
	wp: 3 x 2 matrix with key deflection parameters: 
		1st row is the maximum deflection (maximum basin
		depth) and x at this location
		2nd row is the zero deflection and x at this 
		location (basin width)
		3rd row is the minimum deflection (forebulge) 
		and x at this location
	"""
	# geometry
	extent = geom[0] # extent in x
	xint = geom[1] # interval in x
	x = np.arange(0,extent+1,xint) # points in x
	w = np.zeros(len(x)) # initialize displacement
	
	# elastic and flexural parameters
	E = elas[0] # Young Modulus
	v = elas[1] # Poisson ratio
	h = elas[2] # elastic thickness
	# flexural rigidity
	rigid = (E*h*h*h)/(12*(1.0-v*v)) 
	densup = elas[3] # density of foundation
	g = 9.81 # gravity
	k = densup*g # support of foundation
	# flexural parameter
	alpha = ((4*rigid)/(k))**(1/4) 
	
	# loads
	lxmin = loads[:,0]
	lxmax = loads[:,1]
	lh = loads[:,2]
	ldens = loads[:,3]
	
	# compute deflection profile
	# for all the loads columns
	for i in range(len(lxmin)):
		q = lh[i] * ldens[i] * 9.81
		# for all points in x
		for j in range(len(x)):
			tolf = abs(x[j] - lxmin[i])
			tort = abs(x[j] - lxmax[i])
			dA = np.exp(-tolf/alpha)*np.cos(tolf/alpha)
			dB = np.exp(-tort/alpha)*np.cos(tort/alpha)
			# if below the load
			if x[j] >= lxmin[i] and x[j] <= lxmax[i]:
				w[j] = w[j] + (q/(2.0*k))*(2.0-dA-dB)
			# if to the left of the load
			elif x[j] < lxmin[i]:
				w[j] = w[j] + (q/(2.0*k))*(dA-dB)
			# if to the right of the load
			elif x[j] > lxmax[i]:
				w[j] = w[j] - (q/(2.0*k))*(dA-dB)
	
	# key deflection parameters
	wp = np.zeros((3,2))
	# maximum basin depth
	v = max(w)
	ind = w.argmax()
	wp[0,0] = v
	wp[0,1] = x[ind]
	# basin width
	ind = np.where(w <=0)[0]
	wp[1,0] = 0.0
	wp[1,1] = x[ind[0]]
	# forebulge
	v = min(w)
	ind = w.argmin()
	wp[2,0] = v
	wp[2,1] = x[ind]
	
	# plot
	# input loads
	for i in range(len(lxmin)):
		xcol = [lxmin[i], lxmin[i], lxmax[i], lxmax[i]]
		ycol = [0, lh[i], lh[i], 0]
		ax[0].plot(xcol,ycol,"k-")
	maxlh = max(lh)*1.25
	ax[0].axis([0, geom[0], -maxlh, maxlh])
	ax[0].grid()
	ax[0].set_xlabel("m")
	ax[0].set_ylabel("m")
	ax[0].xaxis.set_major_formatter(ticker.FormatStrFormatter
		("%0.0e"))
	ax[0].set_title("Input loads",fontweight="bold")
	
	# deflected loads plus deflection
	wl = np.interp(lxmin,x,w)
	wr = np.interp(lxmax,x,w)
	for i in range(len(lxmin)):
		xcol = [lxmin[i], lxmin[i], lxmax[i], lxmax[i]]
		ycol = [-wl[i], lh[i]-wl[i], lh[i]-wr[i], -wr[i]]
		ax[1].plot(xcol,ycol,"k-")
	ax[1].plot(x,-w,"b-")
	ax[1].axis([0, geom[0], -maxlh, maxlh])
	ax[1].grid()
	ax[1].set_xlabel("m")
	ax[1].set_ylabel("m")
	ax[1].xaxis.set_major_formatter(ticker.FormatStrFormatter
		("%0.0e"))
	ax[1].set_title("Deflected loads",fontweight="bold")
	fig.subplots_adjust(hspace=0.6)
	
	return w, wp