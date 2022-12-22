import numpy as np

def lscov(A, B, w=None):
	"""Least-squares solution in presence of known covariance
	
	:math:`A \\cdot x = B`, that is, :math:`x` minimizes
	:math:`(B - A \\cdot x)^T \\cdot \\text{diag}(w) \\cdot (B - A \\cdot x)`.
	The matrix :math:`w` typically contains either counts or inverse
	variances.
	
	Parameters
	----------
	A: matrix or 2d ndarray
		input matrix
	B: vector or 1d ndarray
		input vector
	
	Notes
	--------
	https://de.mathworks.com/help/matlab/ref/lscov.html
	Code written by Paul Muller in connection with the ggf package
	"""
	# https://stackoverflow.com/questions/27128688/how-to-use-least-squares-with-weight-matrix-in-python
	# https://de.mathworks.com/help/matlab/ref/lscov.html
	if w is None:
		Aw = A.copy()
		Bw = B.T.copy()
	else:
		W = np.sqrt(np.diag(np.array(w).flatten()))
		Aw = np.dot(W, A)
		Bw = np.dot(B.T, W)
	
	# set rcond=1e-10 to prevent diverging odd indices in x
	# (problem specific to ggf/stress computation)
	x, residuals, rank, s = np.linalg.lstsq(Aw, Bw.T, rcond=1e-10)
	return np.array(x).flatten()