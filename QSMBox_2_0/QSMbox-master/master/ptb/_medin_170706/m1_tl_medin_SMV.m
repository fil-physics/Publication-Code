function y = m1_tl_medin_SMV( iFreq, matrix_size, voxel_size, radius )
% Spherical Mean Value operator
%   y = m1_tl_medin_SMV( iFreq, matrix_size, voxel_size, radius )
%   
%   OUTPUT
%   y               SMV operator
% 
%   INPUT
%   iFreq           Input B_local image
%   matrix_size     Matrix size
%   voxel_size      Voxel resolution in mm
%   radius          Sphere radius in mm
%
%   Created by Tian Liu in 2010
%   Last modified by Tian Liu on 2013.07.24
%   Last modified by Julio Acosta-Cabronero on 2014.10.31

if (nargin<4)
    radius = round(6/max(voxel_size))*max(voxel_size); % default radius
end

y = ifftn( fftn(iFreq).*...
               m1_tl_medin_sphere_kernel(matrix_size,voxel_size,radius) );
