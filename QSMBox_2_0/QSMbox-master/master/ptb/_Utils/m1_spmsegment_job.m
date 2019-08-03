load fname.mat;
load settings.mat;
spmdir = fileread('~/.QSMbox/ptb_spmpath.m');
spmdir = spmdir(9:end);
%-----------------------------------------------------------------------
% Job below saved on 11-Apr-2014 16:38:31 by cfg_util (rev $Rev: 5797 $)
% spm SPM - SPM12b (5918)
%-----------------------------------------------------------------------
matlabbatch{1}.spm.spatial.preproc.channel.vols 	= {[fname ',1']};
if strcmpi(settings,'3T_default')
    matlabbatch{1}.spm.spatial.preproc.channel.biasreg 	= 0.001;
    matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 60;
elseif strcmpi(settings,'7T_default')
    matlabbatch{1}.spm.spatial.preproc.channel.biasreg 	= 1e-05;
    matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 30;
end
matlabbatch{1}.spm.spatial.preproc.channel.write 	= [0 1];
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm 	= {[spmdir '/tpm/TPM.nii,1']};
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus 	= 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native 	= [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped 	= [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm 	= {[spmdir '/tpm/TPM.nii,2']};
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus 	= 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native 	= [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped 	= [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm 	= {[spmdir '/tpm/TPM.nii,3']};
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus 	= 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native 	= [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped 	= [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm 	= {[spmdir '/tpm/TPM.nii,4']};
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus 	= 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native 	= [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped 	= [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm 	= {[spmdir '/tpm/TPM.nii,5']};
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus 	= 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native 	= [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped 	= [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm 	= {[spmdir '/tpm/TPM.nii,6']};
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus 	= 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native 	= [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped 	= [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf 		= 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup 	= 1;
matlabbatch{1}.spm.spatial.preproc.warp.reg 		= [0 0.001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg 		= 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm 		= 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp 		= 3;
matlabbatch{1}.spm.spatial.preproc.warp.write 		= [0 0];
%
delete('fname.mat');
delete('settings.mat');
