import pandas as pd
import numpy as np

from Pole import Pole as Pole
from Bingham import Bingham as Bingham


# load strike and dips
data = pd.read_csv('../data/ch5-6/beasd_sep.txt', header = None)

# strike and dip in radians
S = data[0]*np.pi/180
D = data[1]*np.pi/180

# Initialize poles
T = np.zeros(len(S))
P = np.zeros(len(S))

# Find poles
for i in range(len(S)):
     T[i],P[i] = Pole(S[i],D[i],1) # Pole from plane

# Compute cylindrical best fit
eigVec,confCone,bestFit = Bingham(T,P)

# Write best fit old axis to the console
print(eigVec[2,1]*180/np.pi) # Trend should be 125 deg
print(eigVec[2,2]*180/np.pi) # Plunge should be 26 deg