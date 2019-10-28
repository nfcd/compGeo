# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 13:19:58 2019

@author: A
"""
import math


def ZeroTwoPi(a):
    '''
    a = input azimuth
    b = output azimuth from 0 to 2*pi
    
    NOTE: Azimuths a and b are input/output in radians 
    '''
    b=a
    twopi = 2*math.pi
    if b < 0:
        b += twopi
    elif b >= twopi:
        b -= twopi
    return b
