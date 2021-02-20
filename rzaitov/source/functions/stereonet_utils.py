import math
import numpy as np
import matplotlib.pyplot as plt
from geogr_to_view import geogr_to_view
from small_circle import small_circle
from great_circle import great_circle
from pole_utils import plane_from_pole

# some constants
pi = np.pi
pi_2 = pi / 2.0
east = pi_2
west = 3.0 * east

def plot_stereonet(stereonet, intrad, trdv = 0, plgv = pi / 2):
	if intrad <= 0 or intrad > np.pi:
		raise ValueError('intrad should be in half-open range (0, pi]')

	# Make a larger figure
	plt.rcParams['figure.figsize'] = [15, 7.5]
	plt.axis([-1, 1, -1, 1])
	plt.axis('equal')
	plt.axis('off')

	r = 1.0 # radius pf stereonet
	TH = np.radians(np.arange(0, 361, 1))
	X = r * np.cos(TH)
	Y = r * np.sin(TH)
	plt.plot(X,Y, 'k')

	n_circles = pi_2 / intrad # number of circles in the half of a hemisphere
	n_circles = int(n_circles) - 1 if n_circles.is_integer() else int(math.floor(n_circles))

	# Small circles. Start at the North
	n_trd, n_plg = 0.0, 0.0
	trd, plg = geogr_to_view(n_trd, n_plg, trdv, plgv)

	# possible cone angles
	north_cones = [(i+1)*intrad for i in range(n_circles)]
	south_cones = [pi-(i+1)*intrad for i in range(n_circles)]
	equator_cone = [pi_2]
	cones = north_cones + equator_cone + south_cones
	for cone_angle in cones:
		SC_T, SC_P = small_circle(trd, plg, cone_angle)
		X, Y = stereonet(SC_T, SC_P)
		plt.plot(X, Y, color='gray', linewidth=0.5)

	# Great circles
	western = [(west, pi_2-i*intrad) for i in range(n_circles+1)]
	eastern = [(east, pi_2-i*intrad) for i in range(n_circles+1)]
	vertical = [(east, 0.0)]
	poles = western + vertical + eastern
	for trd, plg in poles:
		trd, plg = geogr_to_view(trd, plg, trdv, plgv)
		strike, dip = plane_from_pole(trd, plg)
		GC_T, CG_P = great_circle(strike, dip)
		X, Y = stereonet(GC_T, CG_P)
		plt.plot(X, Y, color='gray', linewidth=0.5)