import numpy as np
from sph_to_cart import sph_to_cart
from cart_to_sph import cart_to_sph

def rotation(rtrd, rplg, rot):
	# Convert rotation axis to direction cosines. The
	# convention here is X1 = North, X2 = East, X3 = Down
	cn, ce, cd = sph_to_cart(rtrd, rplg)

	# Calculate the transformation matrix a for the rotation
	# Eq. 5.17
	x = 1.0 - np.cos(rot)
	sin_rot = np.sin(rot)
	cos_rot = np.cos(rot)

	return np.array([
		[ cos_rot + cn*cn*x, -cd*sin_rot + cn*ce*x, ce*sin_rot + cn*cd*x ],
		[ cd*sin_rot + ce*cn*x, cos_rot + ce*ce*x, -cn*sin_rot + ce*cd*x ],
		[ -ce*sin_rot + cd*cn*x, cn*sin_rot + cd*ce*x, cos_rot + cd*cd*x ]
	])


def rotate_line(T, trd, plg):
	cned = np.array(sph_to_cart(trd, plg))
	cn, ce, cd = np.dot(T, cned)
	return cart_to_sph(cn, ce, cd)


def rotate_axis(T, trd, plg):
	cned = np.array(sph_to_cart(trd, plg))
	rotated_cned = np.dot(T, cned)
	if rotated_cned[2] < 0:
		rotated_cned *= -1
	cn, ce, cd = rotated_cned
	return cart_to_sph(cn, ce, cd)