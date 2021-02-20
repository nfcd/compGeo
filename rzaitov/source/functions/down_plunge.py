import numpy as np

def down_plunge(trd, plg):
	'''
	down_plunge constructs the down-plunge projection matrix from trend and plunge of the fold axis
	NOTE: trd and plg should be entered in radians
	'''

	return np.array([
		[ np.sin(trd)*np.sin(plg), np.cos(trd)*np.sin(plg), np.cos(plg)  ],
		[ np.cos(trd),             -np.sin(trd),            0            ],
		[ np.sin(trd)*np.cos(plg), np.cos(trd)*np.cos(plg), -np.sin(plg) ]
	])