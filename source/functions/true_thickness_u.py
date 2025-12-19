from uncertainties import umath # From E. Lebigot

def true_thickness_u(stk,dip,top,base):
	"""
	Calculate true thickness with uncertainty propagation
	
	Parameters
    ----------
    stk, dip : UFloat
        Strike and dip and their uncertainties in radians
    top, base : sequence of length 3 coordinates of the top 
		and base points of the unit in East (E), North (N), 
		Up (U) system. Each coordinate has its uncertainty
		(UFloat).

    Returns
    -------
    UFloat
        True thickness with uncertainty
	"""
	# make the transformation matrix from ENU coordinates
	# to SDP coordinates
	sin_str = umath.sin(stk)
	cos_str = umath.cos(stk)
	sin_dip = umath.sin(dip)
	cos_dip = umath.cos(dip)

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
	if t.n < 0: # use nominal value n to check sign
		t *= -1.0
	
	return t