function ima_hires = m1_upsampl(ima,upsampl_vec)
% Upsample image resolution
% 
% ima_hires = m1_upsampl(ima,upsampl_vec)
%
% Created by Julio Acosta-Cabronero

% clear all
% close all
% 
% [ima,st]=jniiload('SMVII4_nfm_lowres');
% 
% upsampl_vec = [4 4 4];

matrix_size             = size(ima);
dim1                    = matrix_size(1);
dim2                    = matrix_size(2);
dim3                    = matrix_size(3);
kfdoff                  = 30; % # of voxels where the filter will drop off from 1 to 0.
[kf1,kf2,kf3]           = ndgrid(m1_cosfilt(dim1,kfdoff),m1_cosfilt(dim2,kfdoff),m1_cosfilt(dim3,kfdoff));
kf                      = min(min(kf1,kf2),kf3); clear kf1 kf2 kf3

new_matrix_size         = matrix_size.*upsampl_vec;
pad                     = (new_matrix_size-matrix_size)/2;
rem_pad                 = rem(new_matrix_size-matrix_size,2);
rem_pad(rem_pad==.5)    = 1;
pad                     = floor(pad);

A = ifftshift(ima);                         clear S
B = fftn(A);                                clear A
C = fftshift(B);                            clear B
K = 1/sqrt(length(C(:)))*C;                 clear C
K_real = kf.*real(K);
K_real_pad = padarray(K_real,pad,'both');   clear K_real
if sum(rem_pad(:))>0
    K_real_pad = padarray(K_real_pad,rem_pad,'pre');
end
K_imag = kf.*imag(K);                       clear K
K_imag_pad = padarray(K_imag,pad,'both');   clear K_imag kf
if sum(rem_pad(:))>0
    K_imag_pad = padarray(K_imag_pad,rem_pad,'pre');
end
K_pad = complex(K_real_pad,K_imag_pad);     clear K_*_pad
A = ifftshift(K_pad);                       clear K_pad
B = ifftn(A);                               clear A
C = fftshift(B);                            clear B
ima_hires = sqrt(length(C(:)))*real(C);     clear C

% jniisave(ima_hires,st,'SMVII4_nfm_hires')