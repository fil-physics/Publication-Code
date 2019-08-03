function y = m1_tl_medin_SMV_kernel(matrix_size,voxel_size,radius)
% Generate a kernel that performs the removal of the spherical mean value
%   y = m1_tl_medin_SMV_kernel(matrix_size,voxel_size, radius)
%   
%   OUTPUT
%   y               kernel
% 
%   INPUT
%   matrix_size     FOV dimensions
%   voxel_size      voxel size in mm
%   radius          sphere radius in mm
%
%   Created by Tian Liu in 2010
%   Modified by Tian on 2011.02.01
%   Modified by Tian on 2011.03.14 The sphere is now rendered.
%   Last modified by Tian Liu on 2013.07.23

y = 1-m1_tl_medin_sphere_kernel(matrix_size,voxel_size,radius);
