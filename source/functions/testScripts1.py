import numpy as np

from DirCosAxes import DirCosAxes as DirCosAxes
from Cauchy import Cauchy as Cauchy
from TransformStress import TransformStress as TransformStress
from PrincipalStress import PrincipalStress as PrincipalStress
from ShearOnPlane import ShearOnPlane as ShearOnPlane

# Stress tensor in principal stress coordinate system
stress = [[50, 0, 0],[ 0, 40, 0],[ 0, 0, 30]]

# trend X1, plunge X1, and trend X3
tX1 = 80*np.pi/180
pX1 = 0*np.pi/180
tX3 = 170*np.pi/180

# Test DirCosAxes
dC = DirCosAxes(tX1,pX1,tX3)

# plane
strike = 40*np.pi/180
dip = 65*np.pi/180

# Test Cauchy
T,pT = Cauchy(stress,tX1,pX1,tX3,strike,dip)

# Test TransformStress
ntX1 = 30*np.pi/180
npX1 = 45*np.pi/180
ntX3 = 210*np.pi/180
nstress = TransformStress(stress,tX1,pX1,tX3,ntX1,npX1,ntX3)

# Test Principal Stress
pstress,dCp = PrincipalStress(stress,tX1,pX1,tX3)

pstress1,dCp1 = PrincipalStress(nstress,ntX1,npX1,ntX3)

# Test shear on Plane
TT,dCTT,R = ShearOnPlane(stress,tX1,pX1,tX3,strike,dip)

TT1,dCTT1,R1 = ShearOnPlane(nstress,ntX1,npX1,ntX3,strike,dip)