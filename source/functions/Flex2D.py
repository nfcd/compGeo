import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

def Flex2D(geom,elas,loads):
	'''
	Flex2D computes the deflection profile produced by 
	a group of load columns on an elastic lithosphere 
	resting on a fluid-like foundation 
	(i.e. astenosphere). The program is based on 
	Hetenyi (1946) solution for an infinite elastic 
	beam on a fluid-like foundation

	USE: w, wp = Flex2D(geom,elas,loads)

	geom: A 1 x 2 vector with the lateral extent of the 
		domain in meters, and the distance between
		points (x interval) in meters where the
		deflection will be computed
	elas: A 1 x 4 vector with the Young Modulus (in Pa), 
		Poisson ration, Elastic thickness in meters, 
		and density of the foundation in kg/m^3
	loads: A nloads x 4 vector with the left x coordinate,
		right x coordinate, height, and density of each
		load column. Lengths should be in meters and 
		density in kg/m^3
	w: The deflection of the lithosphere in meters at the 
		points specified by the extent and x interval
	wp: 3 x 2 matrix with key deflection parameters: 
		1st row is the maximum deflection (maximum basin
		depth) and x at this location
		2nd row is the zero deflection and x at this 
		location (basin width)
		3rd row is the minimum deflection (forebulge) 
		and x at this location
	'''
	# Geometry
	extent = geom[0] # Extent in x
	xint = geom[1] # Interval in x
	x = np.arange(0,extent+1,xint) # points in x
	w = np.zeros(len(x)) # initialize displacement

	# Elastic and flexural parameters
	E = elas[0] # Young Modulus
	v = elas[1] # Poisson Ratio
	h = elas[2] # Elastic thickness
	# Flexural rigidity, Eq. 9.11
	rigid = (E*h*h*h)/(12*(1.0-v*v)) 
	densup = elas[3] # Density of foundation
	g = 9.81 # Gravity
	k = densup*g # Support of foundation
	# Flexural parameter, Eq. 9.13
	alpha = ((4*rigid)/(k))**(1/4) 
	
	# Loads
	lxmin = loads[:,0]
	lxmax = loads[:,1]
	lh = loads[:,2]
	ldens = loads[:,3]

	# Compute deflection profile
	# for all the loads columns
	for i in range(0,len(lxmin)):
		# Eqs. 9.14-9.16
		q = lh[i] * ldens[i] * 9.81
		# for all points in x
		for j in range(0,len(x)):
			tolf = abs(x[j] - lxmin[i])
			tort = abs(x[j] - lxmax[i])
			dA = np.exp(-tolf/alpha)*np.cos(tolf/alpha)
			dB = np.exp(-tort/alpha)*np.cos(tort/alpha)
			# If below the load
			if x[j] >= lxmin[i] and x[j] <= lxmax[i]:
				w[j] = w[j] + (q/(2.0*k))*(2.0-dA-dB)
			# If to the left of the load
			elif x[j] < lxmin[i]:
				w[j] = w[j] + (q/(2.0*k))*(dA-dB)
			# If to the right of the load
			elif x[j] > lxmax[i]:
				w[j] = w[j] - (q/(2.0*k))*(dA-dB)
	
	# Key deflection parameters
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

	# Plot
	# Input loads profile
	fig = plt.figure(figsize=(10,8))
	ax1 = fig.add_subplot(2, 1, 1)
	ax2 = fig.add_subplot(2, 1, 2)
	for i in range(0,len(lxmin)):
		xcol = [lxmin[i], lxmin[i], lxmax[i], lxmax[i]]
		ycol = [0, lh[i], lh[i], 0]
		ax1.plot(xcol,ycol,'k-')
	maxlh = max(lh)*1.25
	ax1.axis([0, geom[0], -maxlh, maxlh])
	ax1.grid()
	ax1.set_xlabel('m')
	ax1.set_ylabel('m')
	ax1.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0e'))
	ax1.set_title('Input loads',fontweight="bold")

	# Deflected loads plus deflection
	wl = np.interp(lxmin,x,w)
	wr = np.interp(lxmax,x,w)
	for i in range(0,len(lxmin)):
		xcol = [lxmin[i], lxmin[i], lxmax[i], lxmax[i]]
		ycol = [-wl[i], lh[i]-wl[i], lh[i]-wr[i], -wr[i]]
		ax2.plot(xcol,ycol,'k-')
	ax2.plot(x,-w,'b-')
	ax2.axis([0, geom[0], -maxlh, maxlh])
	ax2.grid()
	ax2.set_xlabel('m')
	ax2.set_ylabel('m')
	ax2.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0e'))
	ax2.set_title('Deflected loads',fontweight="bold")
	fig.subplots_adjust(hspace=0.6)
	return w, wp
