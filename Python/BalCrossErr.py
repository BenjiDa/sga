# Autogenerated with SMOP 
from smop.core import *
# BalCrossErr.m

    
@function
def BalCrossErr(strat=None,vert=None,kk=None,*args,**kwargs):
    varargin = BalCrossErr.varargin
    nargin = BalCrossErr.nargin

    #BalCrossErr computes the shortening error in area balanced cross sections.
#The algorithm was originally written by Phoebe A. Judge
    
    #   USE: [short,shortp,defa,inw] = BalCrossErr(strat,vert,kk)
    
    #   strat= 1 x 5 vector with east stratigraphic thickness (entry 1),
#          west strat. thickness (entry 2), error on east strat
#          thickness (entry 3), error on west strat thickness 
#          (entry 4), and error on final width (entry 5)
#   vert= number of vertices x 5 vector with x coordinates of vertices 
#         (column 1), y coords of vertices (column 2), errors in x
#         coords of vertices (column 3), errors in y coords of vertices 
#         (column 4), and vertices tags (column 5). The vertices tags are 
#         as follows: 1 = Vertex at decollement, 2 = Vertex at surface, 
#         3 = Vertex at subsurface, 4 = Vertex at eroded hanging-wall
#         cutoff
#   kk = A flag to indicate wether the program computes total errrors 
#        (kk = 0), errors due to stratigraphy only (kk = 1), errors due to
#        vertices at decollement only (kk = 2), errors due to vertices in
#        eroded hanging walls only (kk = 3), errors due to surface vertices
#        (kk = 4), or errors due to subsurface vertices (kk = 5)
#   short = Shortening magnitude and its gaussian and maximum errors
#   shortp = Shortening percentage and its gaussian and maximum errrors
#   defa = Deformed area and its gaussian and maximum errors
#   inw = Initial width and its gaussian and maximum errors
#   
#   NOTE: The user selects the length units of the problem. Typical length
#         units are kilometers
    
    #MATLAB script written by Nestor Cardozo for the book Structural 
#Geology Algorithms by Allmendinger, Cardozo, & Fisher, 2011. If you use
#this script, please cite this as "Cardozo in Allmendinger et al. (2011)"
    
    #Stratigraphic thicknesses
    E1=strat[1]
# BalCrossErr.m:36
    
    W1=strat[2]
# BalCrossErr.m:37
    
    dE1=strat[3]
# BalCrossErr.m:38
    
    dW1=strat[4]
# BalCrossErr.m:39
    
    dx2=strat[5]
# BalCrossErr.m:40
    
    #Vertices
    X=vert[:,1]
# BalCrossErr.m:43
    
    Y=vert[:,2]
# BalCrossErr.m:44
    
    dX=vert[:,3]
# BalCrossErr.m:45
    
    dY=vert[:,4]
# BalCrossErr.m:46
    
    Loc=vert[:,5]
# BalCrossErr.m:47
    
    n=size(vert,1)
# BalCrossErr.m:48
    
    #If only errors due to stratigraphy
    if kk == 1:
        dx2=0.0
# BalCrossErr.m:52
        dX=dot(dX,0.0)
# BalCrossErr.m:53
        dY=dot(dY,0.0)
# BalCrossErr.m:54
        #If only errors due to vertices
    else:
        if kk > 1:
            dE1=0.0
# BalCrossErr.m:57
            dW1=0.0
# BalCrossErr.m:58
            dx2=0.0
# BalCrossErr.m:59
            for i in arange(1,n).reshape(-1):
                #If only errors due to decollement vertices
                if kk == 2:
                    if Loc[i] != 1:
                        dX[i]=0.0
# BalCrossErr.m:64
                        dY[i]=0.0
# BalCrossErr.m:65
                    #If only errors due to eroded hanging walls
                else:
                    if kk == 3:
                        if Loc[i] != 4:
                            dX[i]=0.0
# BalCrossErr.m:70
                            dY[i]=0.0
# BalCrossErr.m:71
                        #If only errors due to surface vertices
                    else:
                        if kk == 4:
                            if Loc[i] != 2:
                                dX[i]=0.0
# BalCrossErr.m:76
                                dY[i]=0.0
# BalCrossErr.m:77
                            #If only errors due to subsurface vertices
                        else:
                            if kk == 5:
                                if Loc[i] != 3:
                                    dX[i]=0.0
# BalCrossErr.m:82
                                    dY[i]=0.0
# BalCrossErr.m:83
    
    #Initialize output variables
    short=zeros(1,3)
# BalCrossErr.m:90
    shortp=zeros(1,3)
# BalCrossErr.m:90
    defa=zeros(1,3)
# BalCrossErr.m:91
    inw=zeros(1,3)
# BalCrossErr.m:91
    #Deformed area
#Calculate area of deformed state
    aX=matlabarray(cat([X],[X[1]]))
# BalCrossErr.m:95
    aY=matlabarray(cat([Y],[Y[1]]))
# BalCrossErr.m:96
    XArea=dot(0.5,(multiply(aX[1:n],aY[2:n + 1]) - multiply(aX[2:n + 1],aY[1:n])))
# BalCrossErr.m:97
    defa[1]=(abs(sum(XArea)))
# BalCrossErr.m:98
    #Calculate gaussian uncertainty of deformed area
    aX=matlabarray(cat([X[n]],[aX]))
# BalCrossErr.m:100
    aY=matlabarray(cat([Y[n]],[aY]))
# BalCrossErr.m:101
    dAx=dot(0.5,(aY[3:n + 2] - aY[1:n]))
# BalCrossErr.m:102
    dAy=dot(0.5,(aX[3:n + 2] - aX[1:n]))
# BalCrossErr.m:103
    delAx=(multiply(dAx,dX)) ** 2
# BalCrossErr.m:104
    delAy=(multiply(dAy,dY)) ** 2
# BalCrossErr.m:105
    #Sum the X and Y components
    SdelAx=sum(delAx)
# BalCrossErr.m:107
    SdelAy=sum(delAy)
# BalCrossErr.m:108
    #take the square root of the sum of individual components
    defa[2]=sqrt(SdelAx + SdelAy)
# BalCrossErr.m:110
    #Calculate maximum uncertainty of deformed area
    dAxM=abs(dAx)
# BalCrossErr.m:112
    dAyM=abs(dAy)
# BalCrossErr.m:112
    delAxM=multiply(dAxM,dX)
# BalCrossErr.m:113
    delAyM=multiply(dAyM,dY)
# BalCrossErr.m:114
    #sum the X and Y components
    SdelAxM=sum(delAxM)
# BalCrossErr.m:116
    SdelAyM=sum(delAyM)
# BalCrossErr.m:117
    #Add everything together to get the maximum uncertainty in Area
    defa[3]=SdelAxM + SdelAyM
# BalCrossErr.m:119
    #Original width
#Calculate the original width assuming constant Area
    inw[1]=defa[1] / (((E1) / 2) + ((W1) / 2))
# BalCrossErr.m:123
    #Calculate final width from the imported polygon
    x21=max(X) - min(X)
# BalCrossErr.m:125
    #Calculate gaussian uncertainty of the original width
    ddA=1 / (((E1) / 2) + ((W1) / 2))
# BalCrossErr.m:127
    
    ddE1=- (dot(2,defa[1])) / ((E1) ** 2 + (dot(dot(2,E1),W1)) + (W1) ** 2)
# BalCrossErr.m:128
    
    ddW1=- (dot(2,defa[1])) / ((E1) ** 2 + (dot(dot(2,E1),W1)) + (W1) ** 2)
# BalCrossErr.m:129
    
    inw[2]=sqrt(((dot(ddA,defa[2])) ** 2) + ((dot(ddE1,dE1)) ** 2) + ((dot(ddW1,dW1)) ** 2))
# BalCrossErr.m:130
    #Calculate maximum uncertainty of the original width
    inw[3]=((abs(dot(ddA,defa[3]))) + (abs(dot(ddE1,dE1))) + (abs(dot(ddW1,dW1))))
# BalCrossErr.m:132
    #Shortening
#Calculate shortening
    short[1]=inw[1] - x21
# BalCrossErr.m:136
    #Calculate gaussian uncertainty in shortening
    short[2]=sqrt((inw[2]) ** 2 + (dx2) ** 2)
# BalCrossErr.m:138
    #Calculate maximum uncertainty in shortening
    short[3]=(inw[3] + dx2)
# BalCrossErr.m:140
    #Calculate percent shortening
    shortp[1]=dot((1 - (x21 / inw[1])),100)
# BalCrossErr.m:142
    #Calculate gaussian uncertainty of percent shortening
    ddx1=x21 / ((inw[1]) ** 2)
# BalCrossErr.m:144
    
    ddx2=- 1 / inw[1]
# BalCrossErr.m:145
    
    shortp[2]=dot(sqrt(((dot(ddx1,inw[2])) ** 2) + ((dot(ddx2,dx2)) ** 2)),100)
# BalCrossErr.m:146
    #Calculate maximum uncertainty in shortening
    shortp[3]=dot((abs(dot(ddx1,inw[3])) + (abs(dot(ddx2,dx2)))),100)
# BalCrossErr.m:148
    return short,shortp,defa,inw
    
if __name__ == '__main__':
    pass
    