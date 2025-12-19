import math

def true_thickness(stk,dip,top,base):
	"""
	Calculate true thickness
	
	Parameters
    ----------
    stk, dip : float
        Strike and dip in radians
    top, base : sequence of length 3 coordinates of the top 
		and base points of the unit in East (E), North (N), 
		Up (U) system. 

    Returns
    -------
    Float
        True thickness
	"""
	# make the transformation matrix from ENU coordinates
	# to SDP coordinates
	sin_str = math.sin(stk)
	cos_str = math.cos(stk)
	sin_dip = math.sin(dip)
	cos_dip = math.cos(dip)

	a = [
        [ sin_str,            cos_str,           0       ],
        [-cos_str*cos_dip,    sin_str*cos_dip,   sin_dip ],
        [-cos_str*sin_dip,    sin_str*sin_dip,  -cos_dip ]
    ]
	
	# transform the top and base points
	# from ENU to SDP coordinates
	topn  = [0, 0, 0]
	basen = [0, 0, 0]
	for i in range(3):
		for j in range(3):
			topn[i] += a[i][j]*top[j]
			basen[i] += a[i][j]*base[j]
	
	# compute the thickness of the unit
	t = basen[2] - topn[2]
	# ensure thickness is positive
	if t < 0: 
		t *= -1.0
	
	return t