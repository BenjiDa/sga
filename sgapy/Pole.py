import math
from SphToCart import SphToCart
from CartToSph import CartToSph
import numpy as np

# trd = 132
# plg = 43

def Pole(trd=None,plg=None,k=None,*args,**kwargs):
    # varargin = Pole.varargin
    # nargin = Pole.nargin

    #Pole returns the pole to a plane or the plane which correspond to a pole
    
    #   k is an integer that tells the program what to calculate.
    
    #   If k = 0, [trd1,plg1] = Pole(trd,plg,k) returns the strike 
#   (trd1) and dip (plg1) of a plane, given the trend (trd) 
#   and plunge (plg) of its pole.
    
    #   If k = 1, [trd1,plg1] = Pole(trd,plg,k) returns the trend
#   (trd1) and plunge (plg1) of a pole, given the strike (trd)
#   and dip (plg) of its plane.
    
    #   NOTE: Input/Output angles are in radians. Input/Output strike 
#   and dip are in right hand rule
    
    #   Pole uses functions ZeroTwoPi, SphToCart and CartToSph
    
    #MATLAB script written by Nestor Cardozo for the book Structural 
#Geology Algorithms by Allmendinger, Cardozo, & Fisher, 2011. If you use
#this script, please cite this as "Cardozo in Allmendinger et al. (2011)"
    #Some constants
    east= math.pi / 2.0
# Pole.m:24
    #Calculate plane given its pole
    if k == 0:
        if plg >= 0.0:
            plg1=east - plg
# Pole.m:29
            dipaz=trd - pi
# Pole.m:30
        else:
            plg1=east + plg
# Pole.m:32
            dipaz=copy(trd)
# Pole.m:33
        #Calculate trd1 and make sure it is between 0 and 2*pi
        trd1=ZeroTwoPi(dipaz - east)
# Pole.m:36
        #Else calculate pole given its plane
    else:
        if k == 1:
            cn,ce,cd=SphToCart(trd,plg,k,nargout=3)
# Pole.m:39
            trd1,plg1=CartToSph(cn,ce,cd,nargout=2)
# Pole.m:40
    
    return trd1,plg1
    
if __name__ == '__main__':
    pass
    