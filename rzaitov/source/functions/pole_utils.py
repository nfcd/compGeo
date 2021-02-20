import numpy as np
from zero_two_pi import zero_two_pi

pi_2 = np.pi/2

def pole_from_plane(strike, dip):
  plunge = pi_2 - dip;
  trend = zero_two_pi(strike - pi_2)
  return trend, plunge

def plane_from_pole(trend, plunge):
  condition = (plunge >= 0)
  dip = np.where(condition, pi_2 - plunge, pi_2 + plunge)
  strike = np.where(condition, trend + pi_2, trend - pi_2)
  strike = zero_two_pi(strike)
  return strike, dip