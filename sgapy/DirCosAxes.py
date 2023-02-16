import numpy as np

from sphtocart import sphtocart


def dircosaxes(tX1,pX1,tX3):

    '''
    dircosaxes calculates the direction cosines of a right handed, orthogonal
    X1,X2,X3 cartesian coordinate system of any orientation with respect to 
    North-East-Down
    
       USE: dC = dircosaxes(tX1,pX1,tX3)
    
       tX1 = trend of X1
       pX1 = plunge of X1
       tX3 = trend of X3
       dC = 3 x 3 matrix containing the direction cosines of X1 (row 1),
            X2 (row 2), and X3 (row 3)
    
       Note: Input angles should be in radians
    
       dircosaxes uses function sphtocart
    
    MATLAB script written by Nestor Cardozo for the book Structural 
    Geology Algorithms by Allmendinger, Cardozo, & Fisher, 2011. If you use
    this script, please cite this as "Cardozo in Allmendinger et al. (2011)"
    
    Converted to python by bmelosh Sept 28th 2022
    '''

    #Some constants
    east = np.pi/2.0
    west = 1.5*np.pi

    #Initialize matrix of direction cosines
    dC = np.zeros((3,3))

    #Direction cosines of X1
    [dC[0,0],dC[0,1],dC[0,2]] = sphtocart(tX1,pX1,0)

    #Calculate plunge of axis 3
    #If axis 1 is horizontal
    if pX1 == 0.0:
        if np.abs(tX1-tX3) == east or np.abs(tX1-tX3) == west: ##Check this line
            pX3 = 0.0
        else:
            pX3 = east

    #Else
    else:
        #From Equation 2.14 and with theta equal to 90 degrees
        pX3 = np.arctan(-(dC[0,0]*np.cos(tX3)+dC[0,1]*np.sin(tX3))/dC[0,2])

    #%Direction cosines of X3
    [dC[2,0],dC[2,1],dC[2,2]] = sphtocart(tX3,pX3,0)

    #Compute direction cosines of X2 by the cross product of X3 and X1
    dC[1,0] = dC[2,1]*dC[0,2] - dC[2,2]*dC[0,1]
    dC[1,1] = dC[2,2]*dC[0,0] - dC[2,0]*dC[0,2]
    dC[1,2] = dC[2,0]*dC[0,1] - dC[2,1]*dC[0,0]

    # Convert X2 to a unit vector
    r = np.sqrt(dC[1,0]*dC[1,0]+dC[1,1]*dC[1,1]+dC[1,2]*dC[1,2])
    for i in range(1,3):
        dC[1,i] = dC[1,i]/r

    return dC

