function [mat,ima_struct] = m1_nifti_load(froot)
% Load NIFTI(.gz) file using SPM
%
% SINTAX
%  [matrix, image_structure] = m1_nifti_load(fileroot_no_extension)
%
% DEPENDENCIES
%  SPM
%
% Created by Julio Acosta-Cabronero

fname0 = dir([froot '.nii*']);
fname = fname0(1).name;
if length(fname0) > 1
    disp(['More than one image with that fileroot - using: ' fname0(1).name])
end

disp(['Load ' fname])

[froot,ext] = m1_fname_split(fname);

if strcmpi(ext,'gz')
    gunzip(fname)
    delete(fname)
    [froot, ext] = m1_fname_split(froot);
    fname = [froot '.' ext];
%     gz_flag = 1;
end

if strcmpi(ext,'nii')
    ima_struct = spm_vol(fname);
    mat = spm_read_vols(ima_struct);
    mat(isnan(mat)==1) = 0;
%     if exist('gz_flag')
%         gzip(fname)
%         if exist([fname '.gz'])
%             delete(fname)
%         end
%     end
end
