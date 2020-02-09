import numpy as np
import matplotlib.pyplot as plt
from Pole import Pole as Pole
from GeogrToView import GeogrToView as GeogrToView
from SmallCircle import SmallCircle as SmallCircle
from GreatCircle import GreatCircle as GreatCircle

def Stereonet(trdv,plgv,intrad,sttype):
    '''
	Stereonet plots an equal angle or equal area stereonet 
	of unit radius in any view direction
	
	trdv = trend of view direction
	plgv = plunge of view direction
	intrad = interval in radians between great or small circles 
	sttype = type of stereonet. 0 = equal angle, 1 = equal area
	
	NOTE: All angles should be entered in radians
	
	Stereonet uses functions Pole, GeogrToView, 
	SmallCircle and GreatCircle
	
	Python function translated from the Matlab function
	Stereonet in Allmendinger et al. (2012)
	'''
    pi = np.pi
    # some constants
    east = pi/2.0
    west = 3.0*east
    
    # Plot stereonet reference circle
    r = 1.0 # radius pf stereonet
    TH = np.arange(0,360,1)*pi/180
    X = r * np.cos(TH)
    Y = r * np.sin(TH)
    plt.plot(X,Y, 'k')
    plt.axis ([-1, 1, -1, 1])
    plt.axis ('equal')
    plt.axis('off')
    plt.hold(True)
    
    # Number of small circles
    nCircles = pi/(intrad*2.0)
    
    # small circles
    # start at the North
    trd = 0.0
    plg = 0.0
    # If view direction is not the default (trdv=0,plgv=90)
    # transform line to view direction
    if trdv == 0.0 and plgv == east:
        trd, plg = GeogrToView(trd,plg,trdv,plgv)
    for i in np.arange(1,nCircles):
        coneAngle = i*intrad
        path1, path2, np1, np2 = SmallCircle(trd,plg,coneAngle,sttype)
        plt.plot(path1[np.arange(1,np1),1], path1[np.arange(1,np1),2], 'b')
        if np2 > 0:
            plt.plot(path2[np.arange(1,np2),1], path2[np.arange(1,np2),2], 'b')
    
    # Great Circle
    for i in np.arange(0,nCircles*2):
        # Western half
        if i <= nCircles:
            # Pole of great circle
            trd = west
            plg = i*intrad
            # Eastern half
        else:
            # Pole of great circle
            trd = east
            plg = (i-nCircles)*intrad
        # If pole is vertical, shift it a little bit
        if plg == east:
            plg = plg * 0.9999
        # If view direction is not the default (trdv=0,plgv=90)
        # transform line to view direction
        if trdv == 0.0 and plgv == east:
            trd, plg = GeogrToView(trd,plg,trdv,plgv)
        # Compute plane from pole
        strike, dip = Pole(trd,plg,0)
        # Plot great circle
        path = GreatCircle(strike,dip,sttype)
        plt.plot(path[:,1], path[:,2], 'b')
    plt.hold(False)
            
        
