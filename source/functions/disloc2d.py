from numpy import arctan,arctan2,sin,cos,sign,sqrt,pi,mod,log

def disloc2d(tip,base,slip,nu,obsx,obsy):
	"""
	This function calculates displacements on a 2D planar 
	fault of finite extent, modeled by two edge dislocations 
	in a homogeneous, isotropic elastic half space
	
	Arguments:
		tip = tuple of (x,y) coordinates of the fault tip
		base = tuple of (x,y) coordinates of the fault base
		slip = slip on the fault, + for reverse, - for normal
		nu = Poisson's ratio (scalar)
		obsx = x coordinates of observation points 
		obsy = y coordinates of observation points
	Returns:
		ux = x components of displ. vectors at obs. points 
		uy = y components of displ. vectors at obs. points
		
	disloc2d uses function displacement
	
	Written by David Oakley (david.o.oakley@uis.no)
	"""
	dip = arctan2(tip[1]-base[1],tip[0]-base[0])
	s1 = slip*cos(dip)
	s2 = slip*sin(dip)
	[ux_part1,uy_part1] = displacement(tip[0],tip[1],s1,s2,
		nu,obsx,obsy)
	[ux_part2,uy_part2] = displacement(base[0],base[1],-s1,-s2,
		nu,obsx,obsy)
	ux = ux_part1+ux_part2
	uy = uy_part1+uy_part2
	return ux,uy
	
def displacement(xi1,xi2,s1,s2,nu,x1,x2):
	"""
	Calculate displacements (u1,u2) in 2D in a half space at
	points (x1,x2) due to an edge dislocation at (xi1,xi2) with 
	slip vector (s1,s2). This uses Eqs from Segall (2010), 
	as corrected in the Errata to that book. Indices 1 and 2 
	correspond to the x and y directions respectively.
	Arguments:
		xi1 = x coordinate of the dislocation tip 
		xi2 = y coordinate of the dislocation tip 
		x1 = x coordinates of observation points 
		x2 = y coordinates of observation points 
		s1 = x component of slip vector 
		s2 = y component of slip vector 
		nu = Poisson"s ratio 
	Returns:
		ux = x components of displ. vectors at obs. points 
		uy = y components of displ. vectors at obs. points 
	"""
	# This occurs if the fault dips to the left
	if sign(s1)==sign(s2): 
		# Flip to the sign convention that Segall"s equations use
		s1 = -s1
		s2 = -s2
	# These equations are written for xi1 = 0. If that's not 
	# the case, this makes it equivalent to that case
	x1 = x1-xi1 
	r1_sq = x1**2.+(x2-xi2)**2.
	r2_sq = x1**2.+(x2+xi2)**2.
	r1 = sqrt(r1_sq)
	r2 = sqrt(r2_sq)
	log_r2_r1 = log(r2/r1)
	# Calculate the angles relative to the vertical axis
	theta1 = arctan2(x1,(x2-xi2))
	theta2 = arctan2(x1,(x2+xi2))
	dip = arctan(s2/s1)
	# This puts dip in the range [0,pi]
	if dip<0:
		dip=dip+pi 
	# The following puts the atan branch cuts along the fault 
	# by rotating theta1 to point down the fault and theta2 to 
	# point up (above the half space) along it
	# Shift theta1 to be measured up from the fault
	theta1 = theta1+pi/2.+dip
	# Shift theta2 to be measured from a line pointing up 
	# opposite the fault from the image dislocation
	theta2 = theta2+dip-pi/2. 
	theta1 = mod(theta1,2.*pi)
	theta2 = mod(theta2,2.*pi)
	# Make a correction for rounding errors that can occur 
	# when theta1 or theta2 is very close to 0 or 2*pi.
	theta1[2.*pi-theta1<1e-10] = 0.
	theta2[2.*pi-theta2<1e-10] = 0.
	theta_diff = theta1-theta2
	# These are two very long equations....
	u1 = (-s1/(pi*(1.-nu)))*(((1.-nu)/2.)*(-theta_diff)+x1*(x2-xi2)/(4.*r1_sq) - x1*(x2+(3.-4.*nu)*xi2)/(4.*r2_sq)+xi2*x2*x1*(x2+xi2)/r2_sq**2.) + (s2/(pi*(1.-nu)))*(((1.-2.*nu)/4.)*log_r2_r1-(x2-xi2)**2./(4.*r1_sq) + (x2**2.+xi2**2.-4.*(1.-nu)*xi2*(x2+xi2))/(4.*r2_sq)+x2*xi2*(x2+xi2)**2./r2_sq**2.)
	u2 = (-s1/(pi*(1.-nu)))*(((1.-2.*nu)/4.)*log_r2_r1+(x2-xi2)**2./(4.*r1_sq) - ((x2+xi2)**2.-2.*xi2**2.-2.*(1.-2.*nu)*xi2*(x2+xi2))/(4.*r2_sq) + x2*xi2*(x2+xi2)**2./r2_sq**2.)+(s2/(pi*(1.-nu)))*(((1.-nu)/2.)*(theta_diff) + x1*(x2-xi2)/(4.*r1_sq)-x1*(x2+(3.-4.*nu)*xi2)/(4.*r2_sq)-xi2*x2*x1*(x2+xi2)/r2_sq**2.)
	
	return u1, u2