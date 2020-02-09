import numpy as np
from Pole import Pole as Pole
from Rotate import Rotate as Rotate
from StCoordLine import StCoordLine as StCoordLine

def GreatCircle(strike,dip,sttype):
    '''
	GreatCircle computes the great circle path of a plane 
	in an equal angle or equal area stereonet of unit radius
	
	strike = strike of plane
	dip = dip of plane
	sttype = type of stereonet: 0 = equal angle, 1 = equal area
	path = x and y coordinates of points in great circle path
	
	NOTE: strike and dip should be entered in radians. 
	
	GreatCircle uses functions StCoordLine, Pole and Rotate
	
	Python function translated from the Matlab function
	GreatCircle in Allmendinger et al. (2012)
	'''
    pi = np.pi
	# Compute the pole to the plane. This will be the axis of 
	# rotation to make the great circle
    trda, plga = Pole(strike,dip,1)
    
    
	# Now pick a line at the intersection of the great circle 
	# with the primitive of the stereonet
    trd = strike
    plg = 0.0
    
	# To make the great circle, rotate the line 180 degrees 
	# in increments of 1 degree
    rot = np.arange(0,180,1)*pi/180
    path = np.zeros((rot.np.shape(1),2))
    for i in range(rot.size):
	# Avoid joining ends of path
        if rot[i] == pi:
            rot[i] = rot[i]*0.9999
		# Rotate line
        rtrd, rplg = Rotate(trda,plga,rot(i),trd,plg,'a')
		# Calculate stereonet coordinates of rotated line 
		# and add to great circle path
        path[i,1], path[i,2] = StCoordLine(rtrd,rplg,sttype)
    
    return path
    