function [ mat_disp ] = m1_rot_mat_fig( mat )
% DESCRIPTION
%  Rotate 3D matrix for display (right is right, left is left)
% 
% SINTAX
%  [ matrix4display ] = m1_rot_mat_fig( 3D_image_matrix )
% 
% Created by Julio Acosta-Cabronero

dim = size( mat );

for x=1:dim(3)
    for y=1:dim(1)
        ttt(y,:,x)=mat(dim(1)-y+1,:,x);
    end
    mat_disp(:,:,x) = rot90( ttt(:,:,x) );
end
