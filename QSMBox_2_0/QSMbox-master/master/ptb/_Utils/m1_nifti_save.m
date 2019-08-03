function m1_nifti_save(new_mat,ima_struct,fileroot)
% Save MATLAB array as NIFTI(.gz) file using SPM
%
% SINTAX
%  m1_nifti_save(new_mat,ima_struct,fileroot)
%
% OUTPUT (file)
%  [fileroot '.nii.gz']
% 
% DEPENDENCIES
%  SPM
%
% Created by Julio Acosta-Cabronero

% If image structure was derived from 4D image, use structure of first
%  image only to save 3D output data 

disp(['Save ' fileroot '.nii'])

if length(size(new_mat))==3
    ima_struct = ima_struct(1);
end

fname = [fileroot '.nii'];   
ima_struct.fname = fname;
ima_struct.private.dat.fname = fname;

ima_struct.dim = size(new_mat);
ima_struct.private.dat.dim = size(new_mat);

ima_struct.descrip = 'QSMbox';
ima_struct.private.descrip = 'QSMbox';

spm_write_vol(ima_struct,new_mat);

% gzip(fname)
% if exist([fname '.gz'])
%     delete(fname)
% end
