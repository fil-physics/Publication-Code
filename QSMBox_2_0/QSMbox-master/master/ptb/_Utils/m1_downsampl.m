function ima_lowres = m1_downsampl(ima,downsampl_factor)
% Downsample image resolution
% 
% ima_lowres = m1_downsampl(ima,downsampl_factor)
%
% Created by Julio Acosta-Cabronero

% clear all
% close all
% 
% [ima,st]=jniiload('SMVII4_nfm');
% 
% load ptb_res.txt
% res = ptb_res; 
% clear ptb_res 
% 
% SMV_radius = 4;
% downsampl_factor        = floor(SMV_radius/mean(res(:)));

matrix_size             = size(ima);
kcentre                 = floor(matrix_size/2)+1;
lb                      = kcentre - floor(matrix_size/downsampl_factor/2) ;
ub                      = lb      + ceil(matrix_size/downsampl_factor) - 1;
new_matrix_size         = ub-lb+1;
dim1                    = new_matrix_size(1);
dim2                    = new_matrix_size(2);
dim3                    = new_matrix_size(3);
kfdoff                  = 30; % # of voxels where the filter will drop off from 1 to 0.
[kf1,kf2,kf3]           = ndgrid(m1_cosfilt(dim1,kfdoff),m1_cosfilt(dim2,kfdoff),m1_cosfilt(dim3,kfdoff));
kf                      = min(min(kf1,kf2),kf3); clear kf1 kf2 kf3

A = ifftshift(ima);                         clear ima
B = fftn(A);                                clear A
C = fftshift(B);                            clear B
K = 1/sqrt(length(C(:)))*C;                 clear C
K = K(lb(1):ub(1),lb(2):ub(2),lb(3):ub(3));
K_real = kf.*real(K);
K_imag = kf.*imag(K);                       clear kf
K = complex(K_real,K_imag);                 clear K_real K_imag
A = ifftshift(K);                           clear K
B = ifftn(A);                               clear A
C = fftshift(B);                            clear B
ima_lowres = sqrt(length(C(:)))*real(C);    clear C

% jniisave(ima_lowres,st,'SMVII4_nfm_lowres')