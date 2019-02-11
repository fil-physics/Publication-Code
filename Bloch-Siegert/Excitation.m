%%***************************************************************************************%%
%
%	06.02.18
%   Nadège Corbin, Wellcome Centre for Human Neuroimaging.
%
%%***************************************************************************************%%

function ExcMat=Excitation(alpha,Phi)
% Transformation matrix of the rotation with an angle alpha [rad] about an axis in the
% transverse plane(x,y) with an azymuthal angle of Phi.
Rz = Rotz(-Phi);
Rx = Rotx(alpha);
ExcMat = inv(Rz)*Rx*Rz;