%%***************************************************************************************%%
%
%	06.02.18
%   Nadège Corbin, Wellcome Centre for Human Neuroimaging.
%
%%***************************************************************************************%%

function Rx=Rotx(alpha)
% Transformation matrix of a rotation about the x axis with an angle alpha
% in radian 
Rx=[1 0 0; 0 cos(alpha) -sin(alpha); 0 sin(alpha) cos(alpha)];
end
