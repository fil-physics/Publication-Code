%%***************************************************************************************%%
%
%	06.02.18
%   Nadège Corbin, Wellcome Centre for Human Neuroimaging.
%
%%***************************************************************************************%%

function Ry=Roty(alpha)
% Transformation matrix of a rotation about the y axis with an angle alpha
% in radian 
Ry=[cos(alpha) 0 sin(alpha);0 1 0; -sin(alpha) 0 cos(alpha)];
end
