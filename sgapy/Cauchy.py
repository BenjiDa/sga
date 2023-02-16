import numpy as np

from dircosaxes import dircosaxes
from sphtocart import sphtocart


def cauchy(stress,tX1,pX1,tX3,strike,dip):

    '''
    Given the stress tensor in a X1,X2,X3 coordinate system of any 
    orientation, cauchy computes the X1,X2,X3 tractions on an arbitrarily
    oriented plane 
    
       USE: [T,pT] = cauchy(stress,tX1,pX1,tX3,strike,dip)
    
       stress = Symmetric 3 x 3 stress tensor
       tX1 = trend of X1
       pX1 = plunge of X1
       tX3 = trend of X3
       strike = strike of plane
       dip = dip of plane
       T = 1 x 3 vector with tractions in X1, X2 and X3
       pT = 1 x 3 vector with direction cosines of pole to plane transformed
            to X1,X2,X3 coordinates
    
       NOTE = Plane orientation follows the right hand rule 
              Input/Output angles are in radians
    
       cauchy uses functions dircosaxes, sphtocart, and carttosph
    
    MATLAB script written by Nestor Cardozo for the book Structural 
    Geology Algorithms by Allmendinger, Cardozo, & Fisher, 2011. If you use
    this script, please cite this as "Cardozo in Allmendinger et al. (2011)"
    
    Converted to python by bmelosh Sept 30th 2022
    '''

    #Compute direction cosines of X1,X2,X3
    dC = dircosaxes(tX1,pX1,tX3)

    #Calculate direction cosines of pole to plane
    p = np.zeros((1,3))
    [p[0][0],p[0][1],p[0][2]] = sphtocart(strike,dip,1)

    #Transform pole to plane to stress coordinates X1,X2,X3
    #The transformation matrix is just the direction cosines of X1,X2,X3
    pT = np.zeros((1,3))
    # for i in range(0,3): #np.arange(0,3).reshape(-1):
    #     for j in range(0,3):#np.arange(0,3).reshape(-1):
    #         pT[0][i] = np.dot(dC[i][j], p[0][j]) + pT[0][i]
    pT = p @ dC.transpose() + pT 

    #Convert transformed pole to unit vector
    r = np.sqrt(pT[0][0]*pT[0][0]+pT[0][1]*pT[0][1]+pT[0][2]*pT[0][2])
    #for i in np.arange(0,3).reshape(-1):
    for i in range(0,3):
        pT[0][i] = pT[0][i]/r

    #Calculate the tractions in stress coordinates X1,X2,X3
    T = np.zeros((1,3)) #Initialize T
    #Compute tractions using Cauchy's law (Eq. 6.7b)
    T = pT @ stress.transpose() + T

    return [T,pT]