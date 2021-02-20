import math
two_pi = 2 * math.pi

# This function could accept a number or a numpy array
# Here we exploit that integer division in Python rounds down (not a truncation as in other languages such as C or C++)
# https://docs.python.org/3/library/operator.html#operator.floordiv
def zero_two_pi(a):
	return a - (a // two_pi) * two_pi