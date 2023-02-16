#import math

import numpy as np
from sphtocart import sphtocart
from carttosph import carttosph
from zerotwopi import zerotwopi


def pole(trd, plg, k):

    '''
    function [trd1,plg1] = pole(trd,plg,k)
    pole returns the pole to a plane or the plane which correspond to a pole 
    
       k is an integer that tells the program what to calculate. 
    
       If k = 0, [trd1,plg1] = Pole(trd,plg,k) returns the strike 
       (trd1) and dip (plg1) of a plane, given the trend (trd) 
       and plunge (plg) of its pole.
    
       If k = 1, [trd1,plg1] = Pole(trd,plg,k) returns the trend
       (trd1) and plunge (plg1) of a pole, given the strike (trd)
       and dip (plg) of its plane.
    
       NOTE: Input/Output angles are in radians. Input/Output strike 
       and dip are in right hand rule
    
      Pole uses functions zerotwopi, sphtocart and carttosph
    
    MATLAB script written by Nestor Cardozo for the book Structural 
    Geology Algorithms by Allmendinger, Cardozo, & Fisher, 2011. If you use
    this script, please cite this as "Cardozo in Allmendinger et al. (2011)"
    
    Converted to python by Ben Melosh Sept 28th 2022
    '''

    #Some constants
    east = np.pi/2.0

    #Calculate plane given its pole
    if k == 0:
        if plg >= 0.0:
            plg1 = east - plg
            dipaz = trd - np.pi
        else:
            plg1 = east + plg
            dipaz = trd
        #Calculate trd1 and make sure it is between 0 and 2*pi
        trd1 = zerotwopi(dipaz - east)

        return [trd1,plg1]

    #Else calculate pole given its plane
    elif k == 1:
        [cn,ce,cd] = sphtocart(trd,plg,k)
        [trd1,plg1] = carttosph(cn,ce,cd)

        return [trd1,plg1]