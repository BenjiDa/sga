import numpy as np

from dircosaxes import dircosaxes
from carttosph import carttosph


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