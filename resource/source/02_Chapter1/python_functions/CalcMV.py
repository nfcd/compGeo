# -*- coding: utf-8 -*-
import math

def CalcMV(T, P):
    '''
    CalcMV calculates the mean vector for a given series of lines

    CalcMV(T,P) calculates the trend (trd) 
    and plunge (plg) of the mean vector, its normalized length (Rave), and
    Fisher statistics (concentration factor (conc), 99 (d99) and 
    95 (d95) % uncertainty cones); for a series of lines whose trends
    and plunges are stored in the vectors T and P

    NOTE: Input/Output trends and plunges, as well as uncertainty
    cones are in radians

    CalcMV uses functions SphToCart and CartToSph
    
    Python function translated from the Matlab function CalcMV
    of Structural Geology Algorithms by Allmendinger, Cardozo, & Fisher, 2011
    '''
    # Number of lines
    nlines = len(T) #???
    
    # Initialize the 3 direction cosines which contain the sums of the
    # individual vectors (i.e. the coordinates of the resultant vector)
    CNsum = 0.0
    CEsum = 0.0
    CDsum = 0.0
    
    #Now add up all the individual vectors
    for i in range(nlines):
        cn,ce,cd = SphToCart(T[i],P[i],0)
        CNsum += cn
        CEsum += ce
        CDsum += cd
    
    # R is the length of the resultant vector and Rave is the lenght of
    # the resultant vector normalized by the number of lines
    R = math.sqrt(CNsum*CNsum + CEsum*CEsum + CDsum*CDsum)
    Rave = R/nlines
    # If Rave is lower than 0.1, the mean vector is insignificant, return error
    if Rave < 0.1:
        print('Error: Mean vector is insignificant.')
    #Else 
    else:
        # Divide the resultant vector by its length to get the average
        # unit vector
        CNsum = CNsum/R
        CEsum = CEsum/R
        CDsum = CDsum/R
        # Use the following 'if' statement if you want to convert the
        # mean vector to the lower hemisphere
        if CDsum < 0.0:
            CNsum = -CNsum
            CEsum = -CEsum
            CDsum = -CDsum
        # Convert the mean vector from direction cosines to trend and plunge
        trd, plg = CartToSph(CNsum,CEsum,CDsum)
        # If there are enough measurements calculate the Fisher Statistics
        # For more information on these statistics see Fisher et al. (1987)
        if R < nlines:
            if nlines < 16:
                afact = 1.0-(1.0/nlines)
                conc = (nlines/(nlines-R))*afact**2
            else:
                conc = (nlines-1.0)/(nlines-R)
        if Rave >= 0.65 and Rave < 1.0:
            afact = 1.0/0.01
            bfact = 1.0/(nlines-1.0)
            print(R, afact, bfact)
            d99 = math.acos(1.0-((nlines-R)/R)*(afact**bfact-1.0))
            afact = 1.0/0.05
            d95 = math.acos(1.0-((nlines-R)/R)*(afact**bfact-1.0))
    return trd, plg, Rave, conc, d99, d95
