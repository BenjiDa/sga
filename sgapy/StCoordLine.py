import numpy as np

def stcoordline(trd, plg, sttype):

    '''
    stcoordline computes the coordinates of a line 
    in an equal angle or equal area stereonet of unit radius
     
        USE: [xp,yp] = stcoordline(trd,plg,sttype)
     
        trd = trend of line
        plg = plunge of line
        sttype = An integer indicating the type of stereonet. 0 for equal angle
                 and 1 for equal area
        xp and yp = Coordinates of the line in the stereonet plot
     
        NOTE: trend and plunge should be entered in radians
     
        stcoordlLine uses function zerotwopi
     
     MATLAB script written by Nestor Cardozo for the book Structural 
     Geology Algorithms by Allmendinger, Cardozo, & Fisher, 2011. If you use
     this script, please cite this as "Cardozo in Allmendinger et al. (2011)"
    
     Converted to python by Ben Melosh Sept 28th, 2022
     '''


    # Some constants
    piS4 = np.pi/4.0
    s2 = np.sqrt(2.0)
    plgS2 = plg/2.0

    # Take care of negative plunges
    if plg < 0.0:
        trd = zerotwopi(trd+np.pi)
        plg = -plg

    # Equal angle stereonet: From Equation 1.5 above
    # Also see Pollard and Fletcher (2005), eq.2.72
    if sttype == 0:
        xp = np.tan(piS4 - plgS2)*np.sin(trd)
        yp = np.tan(piS4 - plgS2)*np.cos(trd)
    # Equal area stereonet: From Equation 1.6 above
    # Also see Pollard and Fletcher (2005), eq.2.90
    elif sttype == 1:
        xp = s2*np.sin(piS4 - plgS2)*np.sin(trd)
        yp = s2*np.sin(piS4 - plgS2)*np.cos(trd)


    return [xp, yp]