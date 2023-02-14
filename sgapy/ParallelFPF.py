# Autogenerated with SMOP 
from smop.core import *
# ParallelFPF.m

    
@function
def ParallelFPF(yp=None,psect=None,pramp=None,pslip=None,*args,**kwargs):
    varargin = ParallelFPF.varargin
    nargin = ParallelFPF.nargin

    #ParallelFPF plots the evolution of a simple step, parallel
#fault propagation fold
    
    #   USE: frames = ParallelFPF(yp,psect,pramp,pslip)
    
    #   yp = Datums or vertical coordinates of undeformed, horizontal beds
#   psect = A 1 x 2 vector containing the extent of the section, and the 
#           number of points in each bed
#   pramp = A 1 x 2 vector containing the x coordinate of the lower bend in
#           the decollement, and the ramp angle
#   pslip = A 1 x 2 vector containing the total and incremental slip
#   frames = An array structure containing the frames of the fold evolution
#            You can play the movie again just by typing movie(frames)
#   
#   NOTE: Input ramp angle should be in radians
    
    #   ParallelFPF uses function SuppeEquationTwo
    
    #MATLAB script written by Nestor Cardozo for the book Structural 
#Geology Algorithms by Allmendinger, Cardozo, & Fisher, 2011. If you use
#this script, please cite this as "Cardozo in Allmendinger et al. (2011)"
    
    # Base of layers
    base=yp[1]
# ParallelFPF.m:25
    #Extent of section and number of points in each bed
    extent=psect[1]
# ParallelFPF.m:28
    npoint=psect[2]
# ParallelFPF.m:28
    #Make undeformed beds geometry: This is a grid of points along the beds
    xp=arange(0.0,extent,extent / npoint)
# ParallelFPF.m:30
    XP,YP=meshgrid(xp,yp,nargout=2)
# ParallelFPF.m:31
    #Fault geometry and slip
    xramp=pramp[1]
# ParallelFPF.m:34
    ramp=pramp[2]
# ParallelFPF.m:34
    slip=pslip[1]
# ParallelFPF.m:35
    sinc=pslip[2]
# ParallelFPF.m:35
    #Number of slip increments
    ninc=round(slip / sinc)
# ParallelFPF.m:37
    #Solve model parameters
#Solve first equation in Eq. 11.20 by minimizing SuppeEquationTwo
    options=optimset('display','off')
# ParallelFPF.m:41
    gamstar=fzero('SuppeEquationTwo',0.5,options,ramp)
# ParallelFPF.m:42
    #Solve second equation in Eq. 11.20
    gam1=pi / 2.0 - ramp / 2.0
# ParallelFPF.m:44
    #Solve third equation in Eq. 11.20
    gam=pi / 2.0 + gamstar - gam1
# ParallelFPF.m:46
    #Solve fourth equation in Eq. 11.20
    bet2=pi - dot(2.0,gamstar)
# ParallelFPF.m:48
    #Other angle for computation
    kap=pi - bet2 + ramp
# ParallelFPF.m:50
    #Eq. 11.21
    lbrat=1.0 / (1.0 - sin(ramp) / sin(dot(2.0,gam) - ramp))
# ParallelFPF.m:52
    #Eq. 11.23
    R1=sin(gam1 + ramp) / sin(gam1 + gam)
# ParallelFPF.m:54
    R2=sin(bet2) / sin(bet2 - ramp + gam)
# ParallelFPF.m:55
    #From the origin of each bed compute the number of points that are in the
#hanging wall. These points are the ones that will move. Notice that this
#has to bee done for each slip increment, since the fault propagates
    hwid=zeros(ninc,size(yp,2))
# ParallelFPF.m:60
    for i in arange(1,ninc).reshape(-1):
        uplift=dot(dot(dot(lbrat,i),sinc),sin(ramp))
# ParallelFPF.m:62
        for j in arange(1,size(yp,2)).reshape(-1):
            if yp[j] - base <= uplift:
                hwid[i,j]=0
# ParallelFPF.m:65
                for k in arange(1,size(xp,2)).reshape(-1):
                    if xp[k] <= xramp + (yp[j] - base) / tan(ramp):
                        hwid[i,j]=hwid[i,j] + 1
# ParallelFPF.m:68
            else:
                hwid[i,j]=size(xp,2)
# ParallelFPF.m:72
    
    #Deform beds: Apply velocity fields of Eq. 11.22
#Loop over slip increments
    for i in arange(1,ninc).reshape(-1):
        # Compute uplift
        lb=dot(dot(lbrat,i),sinc)
# ParallelFPF.m:81
        uplift=dot(lb,sin(ramp))
# ParallelFPF.m:82
        lbh=dot(lb,cos(ramp))
# ParallelFPF.m:83
        ef=uplift / sin(dot(2.0,gamstar))
# ParallelFPF.m:85
        xt=xramp + lbh
# ParallelFPF.m:87
        yt=base + uplift
# ParallelFPF.m:88
        xe=xt + dot(ef,cos(kap))
# ParallelFPF.m:90
        ye=yt + dot(ef,sin(kap))
# ParallelFPF.m:91
        for j in arange(1,size(XP,1)).reshape(-1):
            #Loop over number of hanging wall points in each bed
            for k in arange(1,hwid[i,j]).reshape(-1):
                #If point is in domain 1
                if XP[j,k] < xramp - (YP[j,k] - base) / tan(gam1):
                    XP[j,k]=XP[j,k] + sinc
# ParallelFPF.m:98
                else:
                    # if y lower than y at e
                    if YP[j,k] < ye:
                        #If point is in domain 2
                        if XP[j,k] < xt + (YP[j,k] - yt) / tan(kap):
                            XP[j,k]=XP[j,k] + dot(sinc,cos(ramp))
# ParallelFPF.m:104
                            YP[j,k]=YP[j,k] + dot(sinc,sin(ramp))
# ParallelFPF.m:105
                        else:
                            #If point is in domain 4
                            if XP[j,k] < xt + (YP[j,k] - yt) / tan(gam):
                                XP[j,k]=XP[j,k] + dot(dot(sinc,R2),cos(gam))
# ParallelFPF.m:109
                                YP[j,k]=YP[j,k] + dot(dot(sinc,R2),sin(gam))
# ParallelFPF.m:110
                        # if y higher than y at e
                    else:
                        #If point is in domain 2
                        if XP[j,k] < xe - (YP[j,k] - ye) / tan(gam1):
                            XP[j,k]=XP[j,k] + dot(sinc,cos(ramp))
# ParallelFPF.m:117
                            YP[j,k]=YP[j,k] + dot(sinc,sin(ramp))
# ParallelFPF.m:118
                        else:
                            #If point is in domain 3
                            if XP[j,k] < xe + (YP[j,k] - ye) / tan(gam):
                                XP[j,k]=XP[j,k] + dot(dot(sinc,R1),cos(gam))
# ParallelFPF.m:122
                                YP[j,k]=YP[j,k] + dot(dot(sinc,R1),sin(gam))
# ParallelFPF.m:123
                            else:
                                #If point is in domain 4
                                if XP[j,k] < xt + (YP[j,k] - yt) / tan(gam):
                                    XP[j,k]=XP[j,k] + dot(dot(sinc,R2),cos(gam))
# ParallelFPF.m:127
                                    YP[j,k]=YP[j,k] + dot(dot(sinc,R2),sin(gam))
# ParallelFPF.m:128
        #Plot increment
    #Fault
        xf=matlabarray(cat(0,xramp,xramp + lbh))
# ParallelFPF.m:139
        yf=matlabarray(cat(base,base,uplift + base))
# ParallelFPF.m:140
        plot(xf,yf,'r-','LineWidth',2)
        hold('on')
        for j in arange(1,size(yp,2)).reshape(-1):
            #If beds cut by the fault
            if yp[j] - base <= uplift:
                plot(XP[j,1:1:hwid[i,j]],YP[j,1:1:hwid[i,j]],'k-')
                plot(XP[j,hwid[i,j] + 1:1:size(xp,2)],YP[j,hwid[i,j] + 1:1:size(xp,2)],'k-')
            else:
                plot(XP[j,:],YP[j,:],'k-')
        #Plot settings
        text(dot(0.8,extent),dot(1.75,max(yp)),strcat('Slip = ',num2str(dot(i,sinc))))
        axis('equal')
        axis(cat(0,extent,0,dot(2.0,max(yp))))
        hold('off')
        frames[i]=getframe
# ParallelFPF.m:161
    
    return frames
    
if __name__ == '__main__':
    pass
    