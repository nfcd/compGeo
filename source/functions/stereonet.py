import numpy as np
import matplotlib.pyplot as plt
from pole import plane_from_pole
from geogr_to_view import geogr_to_view
from small_circle import small_circle
from great_circle import great_circle

def stereonet(trdv,plgv,intrad,stype):
	"""
	stereonet plots an equal angle or equal area stereonet
	of unit radius in any view direction
	
	USE: fig, ax = stereonet(trdv,plgv,intrad,stype)
	
	trdv = trend of view direction
	plgv = plunge of view direction
	intrad = interval in radians between great or small circles
	stype = Stereonet type: 0 = equal angle, 1 = equal area
	
	fig and ax are handles to the figure and axes
	
	NOTE: All angles should be entered in radians
	
	Python function translated from the Matlab function
	Stereonet in Allmendinger et al. (2012)
	"""
	pi = np.pi
	# some constants
	east = pi/2.0
	west = 3.0*east
	
	# Make a figure
	fig, ax = plt.subplots(figsize=(15,7.5))
	
	# Stereonet reference circle
	r = 1.0 # radius of stereonet
	th = np.radians(np.arange(0,361,1))
	x = r * np.cos(th)
	y = r * np.sin(th)
	# Plot stereonet reference circle
	ax.plot(x,y, "k")
	ax.axis ([-1, 1, -1, 1])
	ax.axis ("equal")
	ax.axis("off")
	
	# Number of small circles
	ncircles = int(pi/(intrad*2.0))
	
	# small circles, start at North
	trd = plg = 0.0
	
	# If view direction is not the default
	# transform line to view direction
	if trdv != 0.0 or plgv != east:
		trd, plg = geogr_to_view(trd,plg,trdv,plgv)
	
	# Plot small circles
	for i in range(1,ncircles+1):
		cangle = i * intrad
		path1, path2, np1, np2 = small_circle(trd,plg,cangle,
			stype)
		ax.plot(path1[np.arange(0,np1),0], 
			path1[np.arange(0,np1),1], color="gray",
			linewidth=0.5)
		if np2 > 0:
			ax.plot(path2[np.arange(0,np2),0], 
				path2[np.arange(0,np2),1], color="gray", 
				linewidth=0.5)
	
	# Great circles
	for i in range(ncircles*2+1):
		# Western half
		if i <= ncircles:
			# Pole of great circle
			trd = west
			plg = i * intrad
		# Eastern half
		else:
			# Pole of great circle
			trd = east
			plg = (i-ncircles) * intrad
		# If pole is vertical shift it a little bit
		if plg == east:
			plg = plg * 0.9999
		# If view direction is not the default 
		# transform line to view direction
		if trdv != 0.0 or plgv != east:
			trd, plg = geogr_to_view(trd,plg,trdv,plgv)
		# Compute plane from pole
		strike, dip = plane_from_pole(trd,plg)
		# Plot great circle
		path = great_circle(strike,dip,stype)
		ax.plot(path[:,0],path[:,1],color="gray",linewidth=0.5)
		
	# return figure and axes handles
	return fig, ax