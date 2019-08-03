function [ima_SMV,ima_edge] = m1_SMVfilt(ima,res,radius)
% Spherical mean value filtering and edge detection (if radius=1)
%
% SINTAX
%  [ima_SMV,ima_edge] = m1_SMVfilt(ima,res,radius)
%
% Based on code from Tian Liu (MEDI toolbox, Cornell)
% Adapted by Julio Acosta-Cabronero, 23 Jun 2017

ima_SMV  = SMV(ima,res,radius);
ima_edge = ima-ima_SMV;

function ima_SMV = SMV(ima,res,radius)

ima_SMV = ifftn(fftn(ima).*sphere_kernel(size(ima),res,radius));

function S = sphere_kernel(matrix_size,voxel_size,radius)

[Y,X,Z]=meshgrid(-matrix_size(2)/2:matrix_size(2)/2-1,...
                 -matrix_size(1)/2:matrix_size(1)/2-1,...
                 -matrix_size(3)/2:matrix_size(3)/2-1);

X = X*voxel_size(1);
Y = Y*voxel_size(2);
Z = Z*voxel_size(3);

Sphere_out = (   max(abs(X)-0.5*voxel_size(1),0).^2 ... 
                +max(abs(Y)-0.5*voxel_size(2),0).^2 ...
                +max(abs(Z)-0.5*voxel_size(3),0).^2 )>radius^2;

Sphere_in = (  (abs(X)+0.5*voxel_size(1)).^2 ... 
              +(abs(Y)+0.5*voxel_size(2)).^2 ...
              +(abs(Z)+0.5*voxel_size(3)).^2 )<=radius^2; 

            
Sphere_mid = zeros(matrix_size);

split = 10; %such that error is controlled at <1/(2*10)
[X_v,Y_v,Z_v] = meshgrid(-split+0.5:split-0.5, -split+0.5:split-0.5, -split+0.5:split-0.5);
X_v = X_v/(2*split);
Y_v = Y_v/(2*split);
Z_v = Z_v/(2*split);

shell = 1-Sphere_in-Sphere_out;
X = X(shell==1);
Y = Y(shell==1);
Z = Z(shell==1);
shell_val = zeros(size(X));

for i = 1:length(X)
    xx = X(i);
    yy = Y(i);
    zz = Z(i);
   
    occupied = ( (xx+X_v*voxel_size(1)).^2+...
        (yy+Y_v*voxel_size(2)).^2+...
        (zz+Z_v*voxel_size(3)).^2)<=radius^2;
    shell_val(i) = sum(occupied(:))/numel(X_v);
end

Sphere_mid(shell==1) = shell_val;

Sphere = Sphere_in+Sphere_mid;    
Sphere = Sphere/sum(Sphere(:));
S = fftn(fftshift(Sphere));
