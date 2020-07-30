import numpy as np

from DirCosAxes import DirCosAxes as DirCosAxes
from CartToSph import CartToSph as CartToSph

def PrincipalStress(stress,tX1,pX1,tX3):
    '''
    Given the stress tensor in a X1X2X3 coordinate system 
    PrincipalStress calculates the principal stresses and 
    their orientations (trend and plunge) 

    USE: pstress,dCp = PrincipalStress(stress,tX1,pX1,tX3)

    stress = Symmetric 3 x 3 stress tensor
    tX1 = trend of X1
    pX1 = plunge of X1
    tX3 = trend of X3
    pstress = 3 x 3 matrix containing the magnitude 
    		  (column 1), trend (column 2), and plunge 
    		  (column 3) of the  maximum (row 1), 
    		  intermediate (row 2), and minimum (row 3) 
    		  principal stresses
    dCp = 3 x 3 matrix with direction cosines of the 
    	principal stress directions: Max. (row 1), 
    	Int. (row 2), and Min. (row 3) with respect 
    	to North-East-Down

    NOTE: Input/Output angles are in radians

    PrincipalStress uses functions DirCosAxes and CartToSph
     
    Python function translated from the Matlab function
    PrincipalStress in Allmendinger et al. (2012)
    '''
     
    # Compute direction cosines of X1X2X3
    dC = DirCosAxes(tX1,pX1,tX3)
    
    # Initialize pstress
    pstress = np.zeros((3,3))
    
    # Calculate the eigenvalues and eigenvectors 
    # of the stress tensor
    D, V = np.linalg.eigh(stress)

    # Fill principal stress magnitudes
    pstress[0,0] = D[2] # Maximum principal stress
    pstress[1,0] = D[1] # Interm. principal stress
    pstress[2,0] = D[0] # Minimum principal stress


    # The direction cosines of the principal stresses are 
    # with respect to the X1X2X3 stress coordinate system, 
    # so they need to be transformed to the NED coordinate 
    # system
    tV = np.zeros((3,3))
    for i in range(3):
        for j in range(3):
            for k in range(3):
                tV[j,i] = dC[k,j]*V[k,i] + tV[j,i]
                
    # Initialize dCp
    dCp = np.zeros((3,3))
    
    # Direction cosines of principal stresses
    for i in range(3):
        for j in range(3):
            dCp[i,j] = tV[j,2-i]
        # Avoid precision issues
        # Make sure the principal axes are unit vectors
        dCp[i,:] = dCp[i,:]/np.linalg.norm(dCp[i,:])
    
    # Trend and plunge of principal stresses
    for i in range(3):
        pstress[i,1],pstress[i,2] = CartToSph(dCp[i,0],dCp[i,1],dCp[i,2])
      
    return pstress,dCp