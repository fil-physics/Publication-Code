%%***************************************************************************************%%
%
%	06.02.18
%   Nadège Corbin, Wellcome Centre for Human Neuroimaging.
%
%%***************************************************************************************%%

function Rz=Rotz(alpha)
% Transformation matrix of a rotation about the z axis with an angle alpha
% in radian 
Rz=[cos(alpha) -sin(alpha) 0; sin(alpha) cos(alpha) 0; 0 0 1];
end
