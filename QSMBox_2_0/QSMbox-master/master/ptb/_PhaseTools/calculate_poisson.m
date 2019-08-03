clear all
close all

addpath ~/QSMbox/master/ptb/_medin_170706
addpath ~/QSMbox/master/ptb/_Utils
addpath ~/QSMbox/master/ptb/_spm12b_nifti
addpath ~/QSMbox/master/ptb/_PhaseTools

load ptb_res.txt
res = ptb_res; 
clear ptb_res 

load ptb_rad2ppm_ref.txt
rad2ppm = ptb_rad2ppm_ref;
clear ptb_rad2ppm_ref

[phase,st] = jniiload('nfm_0');
% [phase,st]  = jniiload(['SMVII' num2str(mean(res(:))) '_nfm']);
[msk]       = jniiload('Mask');

w=0.1/10;
for x=1:20
    chiname              = dir(['qsm_INTEGRAL_*sweep' num2str(x) '.nii']);
    chiroot              = m1_fname_split(chiname.name);
    chi                  = jniiload(chiroot);
    [diff_poisson,left,right]         = m1_phase_poisson(msk.*phase*rad2ppm,msk.*chi,res);
%     diff_poisson_norm(x) = norm(diff_poisson(msk>0),2) / numel(diff_poisson(msk>0))
    diff_poisson_norm(x) = norm(diff_poisson(msk>0),2) / norm(left(msk>0),2)
    
%     figure(1); jfig(diff_poisson,{-w,w})
%     figure(2); subplot 211; jfig(left,{-w,w}); subplot 212; jfig(right,{-w,w})
%     jniisave(diff_poisson,st,['diff_poisson_' num2str(9+x)])
end

save 'diff_poisson_norm_nfm_0.txt' 'diff_poisson_norm' '-ascii'

% !fslmerge -t all_diff_poisson diff_poisson_*