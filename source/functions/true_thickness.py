import numpy as np

def true_thickness(strike,dip,top,base):
	"""
	true_thickness calculates the thickness (t) of a unit
	given the strike (strike) and dip (dip) of the unit,
	and points at its top (top) and base (base)
	
	top and base are 1 x 3 arrays defining the location
	of top and base points in an ENU coordinate system.
	For each one of these arrays, the first, second
	and third entries are the E, N and U coordinates
	
	NOTE: strike and dip should be input in radians
	"""
	# make the transformation matrix a 
	# from ENU coordinates to SDP coordinates
	sin_str = np.sin(strike)
	cos_str = np.cos(strike)
	sin_dip = np.sin(dip)
	cos_dip = np.cos(dip)
	a = np.array([[sin_str, cos_str, 0],
	[cos_str*cos_dip, -sin_str*cos_dip, -sin_dip],
	[-cos_str*sin_dip, sin_str*sin_dip, -cos_dip]])
	
	# transform the top and base points
	# from ENU to SDP coordinates
	topn = np.zeros(3)
	basen = np.zeros(3)
	for i in range(0,3):
		for j in range(0,3):
			topn[i] = a[i,j]*top[j] + topn[i]
			basen[i] = a[i,j]*base[j] + basen[i]
	
	# compute the thickness of the unit
	t = np.abs(basen[2] - topn[2])
	
	return t