import numpy as np


'''
This file contains functions from the Cardozo and Allmendinger 2011 text book, 
"Structural Geology and Algorithms" converted into python.
'''


def zerotwopi(a):

    '''
     zerotwopi constrains azimuth to lie between 0 and 2*pi radians
    
       b = ZeroTwoPi(a) returns azimuth b (from 0 to 2*pi)
       for input azimuth a (which may not be between 0 to 2*pi)
    
       NOTE: Azimuths a and b are input/output in radians 
    
    MATLAB script written by Nestor Cardozo for the book Structural 
    Geology Algorithms by Allmendinger, Cardozo, & Fisher, 2011. If you use
    this script, please cite this as "Cardozo in Allmendinger et al. (2011)"
    
    Converted to python by bmelosh Sept 27th 2022
    '''

    b=a

    twopi = 2.0*np.pi

    if b < 0.0:
        b = b + twopi
    elif b >= twopi:
        b = b - twopi

    return b

def sphtocart(trd,plg,k):
    

    '''
    sphtocart converts from spherical to cartesian coordinates 
    
       [cn,ce,cd] = SphToCart(trd,plg,k) returns the north (cn), 
       east (ce), and down (cd) direction cosines of a line.
    
       k is an integer to tell whether the trend and plunge of a line 
       (k = 0) or strike and dip of a plane in right hand rule 
       (k = 1) are being sent in the trd and plg slots. In this 
       last case, the direction cosines of the pole to the plane 
       are returned
    
       NOTE: Angles should be entered in radians 
    
    MATLAB script written by Nestor Cardozo for the book Structural 
    Geology Algorithms by Allmendinger, Cardozo, & Fisher, 2011. If you use
    this script, please cite this as "Cardozo in Allmendinger et al. (2011)"*/
    
    Converted to python by bmelosh Sept 28th 2022
    '''

    #If line (see Table 2.1)
    if k == 0:
        cd = np.sin(plg)
        ce = np.cos(plg) * np.sin(trd)
        cn = np.cos(plg) * np.cos(trd) 
    #Else pole to plane (see Table 2.1)
    elif k == 1:
        cd = np.cos(plg)
        ce = -1*np.sin(plg) * np.cos(trd)
        cn = np.sin(plg) * np.sin(trd)

    return [cn,ce,cd]


def carttosph(cn,ce,cd):

    '''
    carttosph converts from cartesian to spherical coordinates 
    
       [trd,plg] = carttosph(cn,ce,cd) returns the trend (trd)
       and plunge (plg) of a line for input north (cn), east (ce), 
       and down (cd) direction cosines
    
       NOTE: Trend and plunge are returned in radians
    
       CartToSph uses function ZeroTwoPi
    
    MATLAB script written by Nestor Cardozo for the book Structural 
    Geology Algorithms by Allmendinger, Cardozo, & Fisher, 2011. If you use
    this script, please cite this as "Cardozo in Allmendinger et al. (2011)"
    
    Converted to python by Ben Melosh Sept 28th 2022
    '''



    # Plunge (see Table 2.1)
    plg = np.arcsin(cd)

    # Trend
    # If north direction cosine is zero, trend is east or west
    # Choose which one by the sign of the east direction cosine
    if cn == 0.0: 
        if ce < 0.0:
            trd = 3.0/2.0*np.pi #% trend is west
        else:
            trd = np.pi/2.0 # trend is east

    # Else use Table 2.1
    else:
        trd = np.arctan(ce/cn) 
        if cn < 0.0:
            #Add pi 
            trd = trd+np.pi

        # Make sure trd is between 0 and 2*pi
        trd = zerotwopi(trd)

    return [trd,plg]



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



def principalstress(stress,tX1,pX1,tX3):


    '''
    Given the stress tensor in a X1,X2,X3 coordinate system of any 
    orientation, principalstress calculates the principal stresses and their
    orientations (trend and plunge) 
    
       USE: [pstress,dCp] = principalstress(stress,tX1,pX1,tX3)
    
       stress = Symmetric 3 x 3 stress tensor
       tX1 = trend of X1
       pX1 = plunge of X1
       tX3 = trend of X3
       pstress = 3 x 3 matrix containing the magnitude (column 1), trend
                 (column 2), and plunge (column 3) of the maximum (row 1),
                 intermediate (row 2), and minimum (row 3) principal stresses
       dCp = 3 x 3 matrix with direction cosines of the principal stress
             directions: Max. (row 1), Int. (row 2), and Min. (row 3)
    
       NOTE: Input/Output angles are in radians
    
       principalstress uses functions dircosaxes and carttosph
    
    MATLAB script written by Nestor Cardozo for the book Structural 
    Geology Algorithms by Allmendinger, Cardozo, & Fisher, 2011. If you use
    this script, please cite this as "Cardozo in Allmendinger et al. (2011)"

    Converted to python by bmelosh Sept 30th 2022
    '''

    #Compute direction cosines of X1,X2,X3
    dC = dircosaxes(tX1,pX1,tX3)

    #Calculate the eigenvalues and eigenvectors of the stress tensor. Use
    #MATLAB function eig. D is a diagonal matrix of eigenvalues
    #(i.e. principal stress magnitudes), and V is a full matrix whose columns
    #are the corresponding eigenvectors (i.e. principal stress directions)
    #[V,D] = eig(stress);
    #[V,D] = np.linalg.eig(stress) ## my first attempt
    w,v = np.linalg.eig(stress) ## my second attempt


    #Initialize pstress
    if np.iscomplex(v).any() == True:
        pstress = np.zeros((3,3), dtype = np.complex_)
    else:
        pstress = np.zeros((3,3))

    #Fill principal stress magnitudes
    pstress[0][0] = w[2] #Maximum principal stress
    pstress[1][0] = w[1] #Intermediate principal stress
    pstress[2][0] = w[0] #Minimum principal stress

    #The direction cosines of the principal stress tensor are given with
    #respect to X1,X2,X3 stress coordinate system, so they need to be
    #transformed to the North-East-Down coordinate system (e.g. Eq. 3.9)
    if np.iscomplex(v).any() == True:
        tV = np.zeros((3,3), dtype = np.complex_)
    else:
        tV = np.zeros((3,3))

    # I can optimize this better!
    for i in range(0,3):
        for j in range(0,3):
            for k in range(0,3):
                tV[j][i] = np.dot(dC[k][j],v[k][i]) + tV[j][i]

    #Initialize dCp
    if np.iscomplex(v).any() == True:
        dCp = np.zeros((3,3), dtype = np.complex_)
    else:
        dCp = np.zeros((3,3))

    #Trend and plunge of maximum principal stress direction
    dCp[0] = [tV[0][2],tV[1][2],tV[2][2]]
    [pstress[0][1],pstress[0][2]] = carttosph(tV[0][2],tV[1][2],tV[2][2])

    #Trend and plunge of intermediate principal stress direction
    dCp[1] = [tV[0][1],tV[1][1],tV[2][1]]
    [pstress[1][1],pstress[1][2]] = carttosph(tV[0][1],tV[1][1],tV[2][1])

    #Trend and plunge of minimum principal stress direction
    dCp[2] = [tV[0][0],tV[1][0],tV[2][0]]
    [pstress[2][1],pstress[2][2]] = carttosph(tV[0][0],tV[1][0],tV[2][0])


    return [pstress,dCp]


def shearonplane(stress,tX1,pX1,tX3,strike,dip):

    '''
    shearonplane calculates the direction and magnitudes of the normal
    and shear tractions on an arbitrarily oriented plane
    
       USE: [TT,dCTT,R] = shearonplane(stress,tX1,pX1,tX3,strike,dip)
    
       stress = 3 x 3 stress tensor
       tX1 = trend of X1
       pX1 = plunge of X1
       tX3 = trend of X3
       strike = strike of plane
       dip = dip of plane
       TT = 3 x 3 matrix with the magnitude (column 1), trend (column 2) and 
           plunge (column 3) of: normal traction on the plane (row 1), 
           minimum shear traction (row 2), and maximum shear traction (row 3)
       dCTT = 3 x 3 matrix with the direction cosines of unit vectors parallel
             to: normal traction on the plane (row 1), minimum shear traction
             (row 2), and maximum shear traction (row 3)
       R = Stress ratio
    
       NOTE = Input stress tensor does not need to be along principal stress
              directions
              Plane orientation follows the right hand rule 
              Input/Output angles are in radians
    
       shearonplane uses functions principalstress, cauchy, carttosph, 
       and shptocart
    
    MATLAB script written by Nestor Cardozo for the book Structural 
    Geology Algorithms by Allmendinger, Cardozo, & Fisher, 2011. If you use
    this script, please cite this as "Cardozo in Allmendinger et al. (2011)"
    
    Converted to python by Ben Melosh Oct 3rd 2022
    '''

    # are there real or complex numbers
    w,v = np.linalg.eig(stress) 

    # Initialize TT and dCTT
    if np.iscomplex(v).any() == True:
        TT = np.zeros((3,3), dtype = np.complex_)
    else:
        TT = np.zeros((3,3))
    if np.iscomplex(v).any() == True:
        dCTT = np.zeros((3,3), dtype = np.complex_)
    else:
        dCTT = np.zeros((3,3))

    # Compute principal stresses and principal stress directions
    [pstress,dCp] = principalstress(stress,tX1,pX1,tX3)

    # Update stress vector so that it is along principal stress directions
    stress = np.zeros((3,3))
    stress = np.diag(pstress[:,0])

    #  New
    #  Calculate direction cosines of pole to plane
    if np.iscomplex(v).any() == True:
        p = np.zeros((1,3), dtype = np.complex_)
    else:
        p = np.zeros((1,3))
    [p[0][0],p[0][1],p[0][2]] = sphtocart(strike,dip,1)

    # Transform pole to plane to  principal stress coordinates
    if np.iscomplex(v).any() == True:
        pT = np.zeros((1,3), dtype = np.complex_)
    else:
        pT = np.zeros((1,3))
    # for i = 1:3
    #     for j = 1:3
    #         pT(i) = dCp(i,j)*p(j) + pT(i);
    pT = p @ dCp.transpose() + pT

    # Calculate the tractions in principal stress coordinates
    if np.iscomplex(v).any() == True:
        T = np.zeros((1,3), dtype = np.complex_)
    else:
        T = np.zeros((1,3))

    # Compute tractions using Cauchy's law
    #for i = 1:3
     #   for j = 1:3
     #       T(i) = stress(i,j)*pT(j) + T(i);

     # for i,j in range(1:3):
     #    T[i] = stress(i,j)*pT(j) + T(i)
    T = pT @ stress.transpose() + T


    # Find the B axis by the cross product of T cross pT and convert to
    # direction cosines (Eq 6.27)
    if np.iscomplex(v).any() == True:
        B = np.zeros((1,3), dtype = np.complex_)
    else:
        B = np.zeros((1,3))
    B[0][0] = T[0][1]*pT[0][2] - T[0][2]*pT[0][1]
    B[0][1] = T[0][2]*pT[0][0] - T[0][0]*pT[0][2]
    B[0][2] = T[0][0]*pT[0][1] - T[0][1]*pT[0][0]


    # Find the shear direction by the cross product of pT cross B. This will
    # give S in right handed coordinates (Eq. 6.27)
    if np.iscomplex(v).any() == True:
        S = np.zeros((1,3), dtype = np.complex_)
    else:
        S = np.zeros((1,3))
    S[0][0] = pT[0][1]*B[0][2] - pT[0][2]*B[0][1]
    S[0][1] = pT[0][2]*B[0][0] - pT[0][0]*B[0][2]
    S[0][2] = pT[0][0]*B[0][1] - pT[0][1]*B[0][0]

    # New: Convert B and S to unit vectors
    rB = np.sqrt(B[0][0]*B[0][0]+B[0][1]*B[0][1]+B[0][2]*B[0][2])
    rS = np.sqrt(S[0][0]*S[0][0]+S[0][1]*S[0][1]+S[0][2]*S[0][2])

    for i in range(0,3):
        B[0][i] = B[0][i]/rB
        S[0][i] = S[0][i]/rS

    # Now we can write the transformation matrix from principal stress
    # coordinates to plane coordinates (Eq. 6.28)
    if np.iscomplex(v).any() == True:
        aa = np.zeros((3,3), dtype = np.complex_)
    else:
        aa = np.zeros((3,3))
    aa[0] = [pT[0][0],pT[0][1],pT[0][2]]
    aa[1] = [B[0][0],B[0][1],B[0][2]]
    aa[2] = [S[0][0],S[0][1],S[0][2]]
    
    # Calculate stress ratio (Eq. 6.32)
    R = (stress[1][1] - stress[0][0])/(stress[2][2]-stress[0][0])

    # Calculate magnitude of normal and shear tractions (Eq. 6.31)
    for i in range(0,3):
        TT[i][0] = stress[0][0]*aa[0][0]*aa[i][0] + stress[1][1]*aa[0][1]*aa[i][1] + stress[2][2]*aa[0][2]*aa[i][2]

    # To get the orientation of the tractions in North-East-Down coordinates, we
    # need to do a vector transformation between principal stress and
    # North-East-Down coordinates. The transformation matrix are just the
    # direction cosines of the principal stresses in North-East-Down coordinates
    # (Eq. 6.29)

    # for i in range(0,3):
    #     for j in range(0,3):
    #         dCTT[0][i] = dCp[j][i]*pT[0][j] + dCTT[0][i]
    #         dCTT[1][i] = dCp[j][i]*B[0][j] + dCTT[1][i]
    #         dCTT[2][i] = dCp[j][i]*S[0][j] + dCTT[2][i]
    dCTT[0] = pT @ dCp + dCTT[0]
    dCTT[1] = B @ dCp + dCTT[1]
    dCTT[2] = S @ dCp + dCTT[2] 

    #Trend and plunge of traction on plane
    [TT[0][1],TT[0][2]] = carttosph(dCTT[0][0],dCTT[0][1],dCTT[0][2])
    #Trend and plunge of minimum shear direction
    [TT[1][1],TT[1][2]] = carttosph(dCTT[1][0],dCTT[1][1],dCTT[1][2])
    #Trend and plunge of maximum shear direction
    [TT[2][1],TT[2][2]] = carttosph(dCTT[2][0],dCTT[2][1],dCTT[2][2])


    return [TT,dCTT,R]