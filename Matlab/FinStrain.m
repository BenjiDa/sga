function [eps,pstrains,dilat,maxsh] = FinStrain(e,frame)
%FinStrain computes finite strain from an input displacement
%gradient tensor
%
%   USE: [eps,pstrains,dilat,maxsh] = FinStrain(e,frame)
%
%   e = 3 x 3 Lagrangian or Eulerian displacement gradient tensor
%   frame = Reference frame. Enter 0 for undeformed (Lagrangian) state, or
%           1 for deformed (Eulerian) state
%   eps = 3 x 3 Lagrangian or Eulerian strain tensor
%   pstrains = 3 x 3 matrix with magnitude (column 1), trend (column 2) and
%              plunge (column 3) of maximum (row 1), intermediate (row 2),
%              and minimum (row 3) elongations
%   dilat = dilatation
%   maxsh = 1 x 2 vector with max. shear strain and orientation with 
%           respect to maximum principal strain direction. Only valid in 2D
%
%   NOTE: Output angles are in radians
%
%   FinStrain uses function CartToSph
%
%MATLAB script written by Nestor Cardozo for the book Structural 
%Geology Algorithms by Allmendinger, Cardozo, & Fisher, 2011. If you use
%this script, please cite this as "Cardozo in Allmendinger et al. (2011)"

%Initialize variables
eps = zeros(3,3);
pstrains = zeros(3,3);
maxsh = zeros(1,2);

%Compute strain tensor (Eqs. 9.4 and 9.5)
for i=1:3
    for j=1:3
        eps(i,j)=0.5*(e(i,j)+e(j,i));
        for k=1:3
            %If undeformed reference frame: Lagrangian strain tensor
            if frame == 0
                eps(i,j) = eps(i,j) + 0.5*(e(k,i)*e(k,j));
            %If deformed reference frame: Eulerian strain tensor
            elseif frame == 1
                eps(i,j) = eps(i,j) - 0.5*(e(k,i)*e(k,j));
            end
        end
    end
end

%Compute principal elongations and orientations. Here we use the MATLAB
%function eig
[V,D] = eig(eps);

%Principal elongations
for i=1:3
    ind = 4-i;
    %Magnitude
    %If undeformed reference frame: Lagrangian strain tensor (Eq. 9.14)
    if frame == 0
        pstrains(i,1) = sqrt(1.0+2.0*D(ind,ind))-1.0;
    %If deformed reference frame: Eulerian strain tensor (Eq. 9.16)
    elseif frame == 1
        pstrains(i,1) = sqrt(1.0/(1.0-2.0*D(ind,ind)))-1.0;       
    end
    %Orientations
    [pstrains(i,2),pstrains(i,3)] = CartToSph(V(1,ind),V(2,ind),V(3,ind));
end

%dilatation (Eq. 9.18)
dilat = (1.0+pstrains(1,1))*(1.0+pstrains(2,1))*(1.0+pstrains(3,1)) - 1.0;

%Maximum shear strain: This only works if plane strain
lmax = (1.0+pstrains(1,1))^2; %Maximum quadratic elongation
lmin = (1.0+pstrains(3,1))^2; %Minimum quadratic elongation
%Maximum shear strain: Ramsay (1967) Eq. 3.46
maxsh(1,1) = (lmax-lmin)/(2.0*sqrt(lmax*lmin));
%Angle of maximum shear strain with respect to maximum principal strain
%Ragan (1967) Eq. 3.45
%If undeformed reference frame
if frame == 0
    maxsh(1,2) = pi/4.0;
%If deformed reference frame
elseif frame == 1
    maxsh(1,2) = atan(sqrt(lmin/lmax));
end

end