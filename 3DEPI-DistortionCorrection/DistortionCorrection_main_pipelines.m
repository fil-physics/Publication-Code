%% 3DEPI distortion correction piplines at 7T for in-vivo data
% No-correction, B0 field mapping and reversed-PE (fsl_topup) techniques
% V Malekian, FIL Physics

clc;clear;close all;
addpath('/export/home/vmalekian/Matlab_lib/spm12/')

%% Topup LSR and Jacobian

%%3DEPI fsl topup
for sub=1:22

    dir1=  ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/'];
    unix(['fslmerge -t ' dir1 'AP_PA ' dir1 'AP ' dir1 'PA'])
    unix(['topup --imain=' dir1 'AP_PA --datain=param.txt --config=b02b0.cnf --out=' dir1 'topup_AP_PA --fout=' dir1 'topup_field --iout=' dir1 'topup_unwarp --jacout=' dir1 'jac --rbmout=' dir1 'xfm --dfout=' dir1 'topup_warpfield'])
    unix(['applytopup --imain=' (dir1) 'AP --topup=' (dir1) 'topup_AP_PA --datain=param.txt --inindex=1 --out=' (dir1) 'AP_dc_def_jac --method=jac']) %% applytopup with Jacobian correction
    unix(['applytopup --imain=' (dir1) 'AP,' (dir1) 'PA --topup=' (dir1) 'topup_AP_PA --datain=param.txt --inindex=1,2 --out=' dir1 'AP_dc_def']) %% applytopup with LSR correction
    unix(['fslmaths ' dir1 'AP_dc_def -thr 0 ' dir1 'AP_dc_def'])
    unix(['fslmaths ' dir1 'AP_dc_def_jac -thr 0 ' dir1 'AP_dc_def_jac'])
   
end

%%MT-3DEPI fsl topup
for sub=1:22
    dir1=  ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/nc_epi3d_v2d_MTw_WholeBrain_0005/'];

    unix(['fslmerge -t ' dir1 'AP_PA ' dir1 'AP ' dir1 'PA'])
    unix(['topup --imain=' dir1 'AP_PA --datain=param.txt --config=b02b0.cnf --out=' dir1 'topup_AP_PA --fout=' dir1 'topup_field --iout=' dir1 'topup_unwarp --jacout=' dir1 'jac --rbmout=' dir1 'xfm --dfout=' dir1 'topup_warpfield'])
    unix(['applytopup --imain=' (dir1) 'AP --topup=' (dir1) 'topup_AP_PA --datain=param.txt --inindex=1 --out=' (dir1) 'AP_dc_def_jac --method=jac']) %%applytopup with Jacobian correction
    unix(['applytopup --imain=' (dir1) 'AP,' (dir1) 'PA --topup=' (dir1) 'topup_AP_PA --datain=param.txt --inindex=1,2 --out=' dir1 'AP_dc_def']) %%applytopup with LSR correction
    unix(['fslmaths ' dir1 'AP_dc_def -thr 0 ' dir1 'AP_dc_def'])
    unix(['fslmaths ' dir1 'AP_dc_def_jac -thr 0 ' dir1 'AP_dc_def_jac'])
   
end



%% MP2RAGE GM Segmentation 

%%Reorient to standard coordinate
for sub=1:22
    dir_struct=  ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/MP2RAGE_WBIC_0pt65_PAT3_PF68_240Hz_UNI/'];
    dir_struct1=  ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/MP2RAGE_WBIC_0pt65_PAT3_PF68_240Hz_INV1/'];
    dir_struct2=  ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/MP2RAGE_WBIC_0pt65_PAT3_PF68_240Hz_INV2/'];

    unix(['fslreorient2std ' (dir_struct) 'uni.nii ' (dir_struct) 'uni_reo.nii.gz'])
    unix(['gunzip ' (dir_struct) 'uni_reo.nii.gz'])
    unix(['fslreorient2std ' (dir_struct1) 'inv1.nii ' (dir_struct1) 'inv1_reo.nii.gz'])
    unix(['gunzip ' (dir_struct1) 'inv1_reo.nii.gz'])
    unix(['fslreorient2std ' (dir_struct2) 'inv2.nii ' (dir_struct2) 'inv2_reo.nii.gz'])
    unix(['gunzip ' (dir_struct2) 'inv2_reo.nii.gz'])
    
end

%%SPM MP2RAGE background noise removal and tissue segmentation
for sub=1:22
    
matlabbatch{1}.spm.tools.mp2rage.rmbg.INV1 = { ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/MP2RAGE_WBIC_0pt65_PAT3_PF68_240Hz_INV1/inv1_reo.nii,1']};
matlabbatch{1}.spm.tools.mp2rage.rmbg.INV2 = { ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/MP2RAGE_WBIC_0pt65_PAT3_PF68_240Hz_INV2/inv2_reo.nii,1']};
matlabbatch{1}.spm.tools.mp2rage.rmbg.UNI =  { ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/MP2RAGE_WBIC_0pt65_PAT3_PF68_240Hz_UNI/uni_reo.nii,1']};

matlabbatch{1}.spm.tools.mp2rage.rmbg.regularization = 1;
matlabbatch{1}.spm.tools.mp2rage.rmbg.output.prefix = 'clean_';
matlabbatch{1}.spm.tools.mp2rage.rmbg.show = 'yes';
matlabbatch{2}.spm.tools.MRI.MRTool_preproc.MRTool_brain.res_dir = {['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/MP2RAGE_WBIC_0pt65_PAT3_PF68_240Hz_UNI']};

matlabbatch{2}.spm.tools.MRI.MRTool_preproc.MRTool_brain.t1w(1) = cfg_dep('Remove background: Background free image', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{3}.spm.spatial.preproc.channel.vols(1) = cfg_dep('Remove background: Background free image', substruct('.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}, '.','val', '{}',{1}), substruct('.','files'));
matlabbatch{3}.spm.spatial.preproc.channel.biasreg = 0.01;
matlabbatch{3}.spm.spatial.preproc.channel.biasfwhm = 30;
matlabbatch{3}.spm.spatial.preproc.channel.write = [1 1];
matlabbatch{3}.spm.spatial.preproc.tissue(1).tpm = {'/export/data/vmalekian/spm12/tpm/TPM.nii,1'};
matlabbatch{3}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{3}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{3}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{3}.spm.spatial.preproc.tissue(2).tpm = {'/export/data/vmalekian/spm12/tpm/TPM.nii,2'};
matlabbatch{3}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{3}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{3}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{3}.spm.spatial.preproc.tissue(3).tpm = {'/export/data/vmalekian/spm12/tpm/TPM.nii,3'};
matlabbatch{3}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{3}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{3}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{3}.spm.spatial.preproc.tissue(4).tpm = {'/export/data/vmalekian/spm12/tpm/TPM.nii,4'};
matlabbatch{3}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{3}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{3}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{3}.spm.spatial.preproc.tissue(5).tpm = {'/export/data/vmalekian/spm12/tpm/TPM.nii,5'};
matlabbatch{3}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{3}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{3}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{3}.spm.spatial.preproc.tissue(6).tpm = {'/export/data/vmalekian/spm12/tpm/TPM.nii,6'};
matlabbatch{3}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{3}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{3}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{3}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{3}.spm.spatial.preproc.warp.cleanup = 2;
matlabbatch{3}.spm.spatial.preproc.warp.reg = [0 0.0001 0.5 0.05 0.2];
matlabbatch{3}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{3}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{3}.spm.spatial.preproc.warp.samp = 1;
matlabbatch{3}.spm.spatial.preproc.warp.write = [0 0];
matlabbatch{3}.spm.spatial.preproc.warp.vox = NaN;
matlabbatch{3}.spm.spatial.preproc.warp.bb = [NaN NaN NaN
                                              NaN NaN NaN];                                                                           
 spm_jobman('run', matlabbatch);                                                                     
end            

%%Brain mask creation
for k=1:sub
    dir_struct=  ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',k) '/MP2RAGE_WBIC_0pt65_PAT3_PF68_240Hz_UNI/'];
    unix(['fslmaths ' (dir_struct) 'c1clean_uni_reo.nii -add ' (dir_struct) 'c2clean_uni_reo.nii -add ' (dir_struct) 'c3clean_uni_reo.nii -bin ' (dir_struct) 'clean_uni_reo_masc'])
    unix(['fslmaths ' (dir_struct) 'clean_uni_reo_masc.nii -dilM ' (dir_struct) 'clean_uni_reo_masc_dil.nii'])
    unix(['fslmaths ' (dir_struct) 'clean_uni_reo_masc_dil.nii -fillh ' (dir_struct) 'clean_uni_reo_masc_dil_fill.nii'])
    unix(['fslmaths ' (dir_struct) 'clean_uni_reo_masc_dil_fill.nii -ero ' (dir_struct) 'clean_uni_reo_masck.nii'])
    unix(['fslmaths ' (dir_struct) 'clean_uni_reo -mul ' (dir_struct) 'clean_uni_reo_masck '  (dir_struct) 'clean_uni_reo_brain'])
end

%% MT-3DEPI GM segmentation
%%SPM segmentation

 for sub=1:22

    matlabbatch{1}.spm.spatial.preproc.channel.vols = {

          ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/nc_epi3d_v2d_MTw_WholeBrain_0005/AP.nii,1'] 
}

    matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.01;
    matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 30;
    matlabbatch{1}.spm.spatial.preproc.channel.write = [1 1];
    matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {'/export/home/vmalekian/Malab_lib/spm12/tpm/TPM.nii,1'};
    matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
    matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {'/export/home/vmalekian/Malab_lib/spm12/tpm/TPM.nii,2'};
    matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
    matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {'/export/home/vmalekian/Malab_lib/spm12/tpm/TPM.nii,3'};
    matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {'/export/home/vmalekian/Malab_lib/spm12/tpm/TPM.nii,4'};
    matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {'/export/home/vmalekian/Malab_lib/spm12/tpm/TPM.nii,5'};
    matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {'/export/home/vmalekian/Malab_lib/spm12/tpm/TPM.nii,6'};
    matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 2;
    matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.0005 0.25 0.025 0.1];
    matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
    matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
    matlabbatch{1}.spm.spatial.preproc.warp.samp = 1;
    matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];
    matlabbatch{1}.spm.spatial.preproc.warp.vox = NaN;
    matlabbatch{1}.spm.spatial.preproc.warp.bb = [NaN NaN NaN
        NaN NaN NaN];
    
     spm_jobman('run', matlabbatch);

 end

%%Brain extraction
for sub=1:22
    dir1= ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/'];
    dir2= ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/nc_epi3d_v2d_MTw_WholeBrain_0005/'];
    unix(['bet ' dir1 'AP ' dir1 'AP_brain -R -f 0.2 -g 0 -m'])
    unix(['bet ' dir2 'AP ' dir2 'AP_brain -R -f 0.2 -g 0 -m'])
end

  
%% Three correction pipelines in MP2RAGE space %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% 1-No correction, 
% 2-Reversed-PE, 
% 3-B0 field mapping 

%%Creating folder for each analysis pipeline
for sub=1:22

    unix(['mkdir -p /export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_FM_3'])
    unix(['mkdir -p /export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_noc_3'])
    unix(['mkdir -p /export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_top_3'])

end

%% 1-No correction & transform MP2RAGE data to MNI space using non-linear registeration
for sub=1:22
    
    dirs=  ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/MP2RAGE_WBIC_0pt65_PAT3_PF68_240Hz_UNI/'];
    dir1= ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/'];
    dir2= ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/nc_epi3d_v2d_MTw_WholeBrain_0005/'];
    diro = ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_noc_3/'];
    dir_reg_no = ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_no.feat/reg/'];

    unix(['fslmaths ' dirs 'c2clean_uni_reo.nii.gz -thr 0.5 -bin ' dirs 'c2clean_uni_reo_wmseg.nii.gz'])
    unix(['flirt -in ' dir1 'AP -ref ' dir2 'AP -out ' diro 'example_func2initial_highres -omat ' diro 'example_func2initial_highres.mat -cost corratio -dof 6 -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -interp spline']) 
    unix(['epi_reg --epi=' dir2 'AP --t1=' dirs 'clean_uni_reo --t1brain=' dirs 'clean_uni_reo_brain --out=' diro 'initial_highres2highres --wmseg=' dirs 'c2clean_uni_reo_wmseg.nii.gz']) %%initial_highres2highres (no correction)
    unix(['convert_xfm -omat ' diro 'example_func2highres.mat -concat ' diro 'initial_highres2highres.mat ' diro 'example_func2initial_highres.mat'])
    unix(['flirt -ref ' dirs 'clean_uni_reo_brain -in ' dir1 'AP -out ' diro 'example_func2highres -applyxfm -init ' diro 'example_func2highres.mat -interp spline']) %%fMRI2highres (no correction)
    unix(['flirt -in ' dirs 'clean_uni_reo_brain -ref /usr/local/fsl/data/standard/MNI152_T1_1mm_brain -out ' diro 'highres2standard -omat ' diro 'highres2standard.mat -cost corratio -dof 12 -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -interp spline'])

    unix(['fnirt --in=' dirs 'clean_uni_reo.nii --ref=' dir_reg_no 'standard.nii.gz --aff=' diro 'highres2standard.mat --cout=' diro 'highres2standard_warp']) 
    unix(['applywarp -i ' dirs 'clean_uni_reo_brain -r /usr/local/fsl/data/standard/MNI152_T1_1mm_brain -o ' diro 'highres2standard -w ' diro 'highres2standard_warp --interp=spline'])
    unix(['convert_xfm -inverse -omat ' diro 'standard2highres.mat ' diro 'highres2standard.mat'])
    unix(['convert_xfm -omat ' diro 'example_func2highres.mat -concat ' diro 'initial_highres2highres.mat ' diro 'example_func2initial_highres.mat'])
    unix(['flirt -ref ' dirs 'clean_uni_reo_brain -in ' dir1 'AP -out ' diro 'example_func2highres -applyxfm -init ' diro 'example_func2highres.mat -interp spline'])
    unix(['flirt -ref ' dirs 'clean_uni_reo_brain -in ' dir1 'AP_brain_mask -out ' dir1 'AP_brain_mask2highres -applyxfm -init ' dir_reg_no 'example_func2highres.mat -interp nearestneighbour'])

    unix(['convert_xfm -inverse -omat ' diro 'highres2example_func.mat ' diro 'example_func2highres.mat'])
    unix(['convert_xfm -omat ' diro 'example_func2standard.mat -concat ' diro 'highres2standard.mat ' diro 'example_func2highres.mat'])
    unix(['convertwarp --ref=/usr/local/fsl/data/standard/MNI152_T1_1mm_brain --premat=' diro 'example_func2highres.mat --warp1=' diro 'highres2standard_warp --out=' diro 'example_func2standard_warp'])
    unix(['applywarp --ref=/usr/local/fsl/data/standard/MNI152_T1_1mm_brain --in=' dir1 'AP --out=' diro 'example_func2standard --warp=' diro 'example_func2standard_warp --interp=spline'])
    unix(['convert_xfm -inverse -omat ' diro 'standard2example_func.mat ' diro 'example_func2standard.mat'])

end


%% 2-Reversed-PE 
for sub=1:22

    dirs=  ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/MP2RAGE_WBIC_0pt65_PAT3_PF68_240Hz_UNI/'];
    dir1= ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/'];
    dir2= ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/nc_epi3d_v2d_MTw_WholeBrain_0005/'];
    diro = ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_top_3/'];
    dir_nc = ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_noc_3/'];

    unix(['flirt -in ' dir1 'AP_dc_def -ref ' dir2 'AP_dc_def -out ' diro 'example_func2initial_highres -omat ' diro 'example_func2initial_highres.mat -cost corratio -dof 6 -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -interp spline'])
    unix(['epi_reg --epi=' dir2 'AP_dc_def --t1=' dirs 'clean_uni_reo --t1brain=' dirs 'clean_uni_reo_brain --out=' diro 'initial_highres2highres --wmseg=' dirs 'c2clean_uni_reo_wmseg.nii.gz']) %% initial_highres2highres (LSR)
    unix(['convert_xfm -omat ' diro 'example_func2highres.mat -concat ' diro 'initial_highres2highres.mat ' diro 'example_func2initial_highres.mat'])   
    unix(['flirt -ref ' dirs 'clean_uni_reo_brain -in ' dir1 'AP_dc_def_jac -out ' diro 'example_func2highres -applyxfm -init ' diro 'example_func2highres.mat -interp spline']) %% fMRI2highres (Jac)
    unix(['convert_xfm -omat ' diro 'example_func2highres.mat -concat ' diro 'initial_highres2highres.mat ' diro 'example_func2initial_highres.mat'])
    unix(['flirt -ref ' dirs 'clean_uni_reo_brain -in ' dir2 'AP_dc_def_jac -out ' diro 'initial_highres2highres_jac -applyxfm -init ' diro 'initial_highres2highres.mat -interp spline'])  %% initial_highres2highres (Jac)
    unix(['convertwarp --ref=' dirs 'clean_uni_reo_brain --postmat=' diro 'initial_highres2highres.mat --warp1=' dir2 'warpfield_01 --out=' diro 'initial_highres2highres_warp --relout']) %% inverse of this used for ..


end

%% 3-B0 field mapping
for sub=1:22
    
    dirs=  ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/MP2RAGE_WBIC_0pt65_PAT3_PF68_240Hz_UNI/'];
    dir1= ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/'];
    dir2= ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/nc_epi3d_v2d_MTw_WholeBrain_0005/'];
    diro = ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_FM_3/'];
    
    %%B0 field preparation & estimation
    unix(['bet ' diro 'magn1_reo.nii.gz ' diro 'magn1_reo_brain1.nii.gz -B -f 0.5 -g 0 -m'])
    unix(['fslmaths ' diro 'magn1_reo_brain1.nii.gz -ero ' diro 'magn1_reo_brain.nii.gz'])
    unix(['fsl_prepare_fieldmap SIEMENS ' diro 'phase_reo.nii.gz ' diro 'magn1_reo_brain.nii.gz ' diro 'fmap_rads.nii.gz 3.08'])
    unix(['fslmaths ' diro 'fmap_rads.nii.gz -div 6.28 ' diro 'fmap_hz.nii.gz'])
    unix(['fslmaths ' dirs 'c2clean_uni_reo.nii.gz -thr 0.5 -bin ' dirs 'c2clean_uni_reo_wmseg.nii.gz'])
    
    unix(['flirt -in ' dir1 'AP -ref ' dir2 'AP -out ' diro 'example_func2initial_highres -omat ' diro 'example_func2initial_highres.mat -cost corratio -dof 6 -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -interp spline'])
    unix(['epi_reg --epi=' dir2 'AP --t1=' dirs 'clean_uni_reo --t1brain=' dirs 'clean_uni_reo_brain --out=' diro 'initial_highres2highres --fmap=' diro 'fmap_rads --fmapmag=' diro 'magn1_reo --fmapmagbrain=' diro 'magn1_reo_brain --echospacing=0.0001525 --pedir=y- --wmseg=' dirs 'c2clean_uni_reo_wmseg.nii.gz']) %% initial_highres2highres(B0)
    unix(['flirt -ref ' dirs 'clean_uni_reo_brain -in ' dir2 'AP -out ' diro 'initial_highres_distorted2highres -applyxfm -init ' diro 'initial_highres2highres.mat -interp spline'])
    
    unix(['convert_xfm -inverse -omat ' diro 'highres2initial_highres.mat ' diro 'initial_highres2highres.mat'])
    unix(['applywarp -i ' diro 'initial_highres2highres -r ' dir2 'AP --premat=' diro 'highres2initial_highres.mat -o ' diro 'initial_highres_undistorted --interp=spline'])
    unix(['convert_xfm -inverse -omat ' diro 'highres2initial_highres.mat ' diro 'initial_highres2highres.mat'])
    unix(['convert_xfm -omat ' diro 'example_func2highres.mat -concat ' diro 'initial_highres2highres.mat ' diro 'example_func2initial_highres.mat'])
    unix(['convert_xfm -inverse -omat ' diro 'highres2example_func.mat ' diro 'example_func2highres.mat'])
    unix(['convertwarp --ref=' dirs 'clean_uni_reo_brain --premat=' diro 'example_func2initial_highres.mat --warp1=' diro 'initial_highres2highres_warp --out=' diro 'example_func2highres_warp --relout'])
    unix(['applywarp --ref=' dirs 'clean_uni_reo_brain --in=' dir1 'AP --out=' diro 'example_func2highres --warp=' diro 'example_func2highres_warp --interp=spline']) %% fMRI2highres (B0)
%     unix(['convert_xfm -inverse -omat ' diro 'highres2example_func.mat ' diro 'example_func2highres.mat'])
%     unix(['applywarp -i ' dir1 'AP -r ' dir1 'AP -o ' diro 'AP_FM -w ' diro 'example_func2highres_warp --postmat='  diro 'highres2example_func.mat --interp=spline'])
   
end        


%% Distort MP2RAGE data with three approaches in MT-3DEPI and 3DEPI spaces %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% MT-3EDPI space 

for sub=1:22
    dirs=  ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/MP2RAGE_WBIC_0pt65_PAT3_PF68_240Hz_UNI/'];
    dir2=  ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/nc_epi3d_v2d_MTw_WholeBrain_0005/'];
    dir_nc = ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_noc_3/'];
    dir_b0 = ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_FM_3/'];
    dir_bp = ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_top_3/'];
    

    %%No correction (MP2RAGE & GM)
     unix(['convert_xfm -inverse -omat ' dir_nc 'highres2initial_highres.mat ' dir_nc 'initial_highres2highres.mat'])
     unix(['applywarp --ref=' dir2 'AP --in=' dirs 'clean_uni_reo_brain --out=' dir_nc 'clean_uni_reo_brain2EPI --premat=' dir_nc 'highres2initial_highres.mat --interp=spline'])
     unix(['fslmaths ' dir_nc 'clean_uni_reo_brain2EPI -thr 0 ' dir_nc 'clean_uni_reo_brain2EPI'])
     unix(['applywarp --ref=' dir2 'AP --in=' dirs 'c1clean_uni_reo --out=' dir_nc 'c1brain2EPI --premat=' dir_nc 'highres2initial_highres.mat --interp=spline'])
     
    %%B0 filed mapping (MP2RAGE & GM)
    unix(['invwarp -w ' dir_b0 'initial_highres2highres_warp -o ' dir_b0 'highres2intial_highres_warp -r ' dir2 'AP'])
    unix(['applywarp --ref=' dir2 'AP --in=' dirs 'clean_uni_reo_brain --out=' dir_b0 'clean_uni_reo_brain2EPI --warp=' dir_b0 'highres2intial_highres_warp --interp=spline'])
    unix(['fslmaths ' dir_b0 'clean_uni_reo_brain2EPI -thr 0 ' dir_b0 'clean_uni_reo_brain2EPI'])
    unix(['applywarp --ref=' dir2 'AP --in=' dirs 'c1clean_uni_reo --out=' dir_b0 'c1brain2EPI --warp=' dir_b0 'highres2intial_highres_warp --interp=spline'])
     
    %%reversed-PE (MP2RAGE & GM)
    %%JAC intensity correction applied to MP2RAGE
    unix(['convert_xfm -inverse -omat ' dir_bp 'highres2initial_highres.mat ' dir_bp 'initial_highres2highres.mat'])
    unix(['applywarp --ref=' dir2 'AP --in=' dirs 'clean_uni_reo_brain --out=' dir_bp 'clean_uni_reo_brain2EPI_temp --premat=' dir_bp 'highres2initial_highres.mat --interp=spline'])
    unix(['fslmaths ' dir_bp 'clean_uni_reo_brain2EPI_temp -div ' dir2 'jac_01.nii.gz ' dir_bp 'clean_uni_reo_brain2EPI_temp_div'])
    unix(['invwarp -w ' dir2 'topup_warpfield_01 -o ' dir2 'warpfield_01_inv -r ' dir_bp 'clean_uni_reo_brain2EPI_temp_div'])
    unix(['applywarp -r ' dir_bp 'clean_uni_reo_brain2EPI_temp_div -i ' dir_bp 'clean_uni_reo_brain2EPI_temp_div --warp=' dir2 'warpfield_01_inv.nii -o ' dir_bp 'clean_uni_reo_brain2EPI_rev --interp=spline'])
    nix(['fslmaths ' dir_reg_bp 'clean_uni_reo_brain2EPI_rev -thr 0 ' dir_reg_bp 'clean_uni_reo_brain2EPI_rev'])
    %%No intensity correction applied to GM
    unix(['convertwarp --ref=' dirs 'clean_uni_reo_brain --postmat=' dir_bp 'initial_highres2highres.mat --warp1=' dir2 'topup_warpfield_01 --out=' dir_bp 'initial_highres2highres_warp --relout'])
    unix(['invwarp -w ' dir_bp 'initial_highres2highres_warp -o ' dir_bp 'highres2intial_highres_warp -r ' dir2 'AP'])
    unix(['applywarp --ref=' dir2 'AP --in=' dirs 'c1clean_uni_reo --out=' dir_bp 'c1brain2EPI --warp=' dir_bp 'highres2intial_highres_warp --interp=spline'])

end

%% 3EDPI space 

for sub=1:22
    
    dirs=  ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/MP2RAGE_WBIC_0pt65_PAT3_PF68_240Hz_UNI/'];
    dir1= ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/'];
    dir2= ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/nc_epi3d_v2d_MTw_WholeBrain_0005/'];
    dir_nc = ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_noc_3/'];
    dir_b0 = ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_FM_3/'];
    dir_bp = ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_top_3/'];

    %%no correction (MP2RAGE) 
     unix(['applywarp --ref=' dir1 'AP --in=' dirs 'clean_uni_reo_brain --out=' dir_nc 'clean_uni_reo_brain2fMREPI --premat=' dir_nc 'highres2example_func.mat --interp=spline'])
    %%B0 field mapping (MP2RAGE) 
    unix(['invwarp -w ' dir_b0 'example_func2highres_warp -o ' dir_b0 'example_func2highres_warp_inv -r ' dir1 'AP'])
    unix(['applywarp --ref=' dir1 'AP --in=' dirs 'clean_uni_reo_brain --out=' dir_b0 'clean_uni_reo_brain2fMREPI --warp=' dir_b0 'example_func2highres_warp_inv --interp=spline'])
    %%reversed-PE (MP2RAGE)    
    unix(['convert_xfm -inverse -omat ' dir_bp 'example_highres2func.mat ' dir_bp 'example_func2highres.mat'])
    unix(['applywarp --ref=' dir1 'AP --in=' dirs 'clean_uni_reo_brain --out=' dir_bp 'clean_uni_reo_brain2fMREPI_temp --premat=' dir_bp 'example_highres2func.mat --interp=spline'])
    unix(['fslmaths ' dir_bp 'clean_uni_reo_brain2fMREPI_temp -div ' dir1 'jac_01.nii.gz ' dir_bp 'clean_uni_reo_brain2fMREPI_temp_div'])
    unix(['invwarp -w ' dir1 'topup_warpfield_01 -o ' dir1 'warpfield_01_inv -r ' dir_bp 'clean_uni_reo_brain2fMREPI_temp_div'])
    unix(['applywarp -r ' dir_bp 'clean_uni_reo_brain2fMREPI_temp_div -i ' dir_bp 'clean_uni_reo_brain2fMREPI_temp_div --warp=' dir1 'warpfield_01_inv.nii -o ' dir_bp 'clean_uni_reo_brain2fMREPI_rev --interp=spline'])

end

%% Distorted MP2RAGE segmentation in MT-3DEPI space
%% SPM tissue segmentation

%%unzip
for sub=1:22

    dir_reg_nc= ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_noc_3/'];
    dir_reg_bp= ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_top_3/'];
    dir_reg_b0= ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_FM_3/'];
    
    string_nc = 'clean_uni_reo_brain2EPI.nii.gz'; %% no correction data in MT-3DEPI
    unix(['gunzip ' (dir_reg_nc) string_nc])
     
    string_bp = 'clean_uni_reo_brain2EPI_rev.nii.gz'; %% no correction data in MT-3DEPI
    unix(['gunzip ' (dir_reg_bp) string_bp])
  
    string_b0 = 'clean_uni_reo_brain2EPI.nii.gz'; %% no correction data in MT-3DEPI
    unix(['gunzip ' (dir_reg_b0) string_nc])

end

%%SPM segmentation 
LL=1:22;
for sub= LL(~ismember(LL,[4,7])) 

    matlabbatch{1}.spm.spatial.preproc.channel.vols = {

          ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_noc_3/clean_uni_reo_brain2EPI.nii,1'] 
          ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_top_3/clean_uni_reo_brain2EPI_rev.nii,1'] 
          ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_FM_3/clean_uni_reo_brain2EPI.nii,1'] 

}

matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.01;
matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 30;
matlabbatch{1}.spm.spatial.preproc.channel.write = [1 1];
matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {'/export/home/vmalekian/Malab_lib/spm12/tpm/TPM.nii,1'};
matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {'/export/home/vmalekian/Malab_lib/spm12/tpm/TPM.nii,2'};
matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {'/export/home/vmalekian/Malab_lib/spm12/tpm/TPM.nii,3'};
matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {'/export/home/vmalekian/Malab_lib/spm12/tpm/TPM.nii,4'};
matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {'/export/home/vmalekian/Malab_lib/spm12/tpm/TPM.nii,5'};
matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {'/export/home/vmalekian/Malab_lib/spm12/tpm/TPM.nii,6'};
matlabbatch{1}.spm.spatial.preproc.tissue(6).ngaus = 2;
matlabbatch{1}.spm.spatial.preproc.tissue(6).native = [0 0];
matlabbatch{1}.spm.spatial.preproc.tissue(6).warped = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.mrf = 1;
matlabbatch{1}.spm.spatial.preproc.warp.cleanup = 2;
matlabbatch{1}.spm.spatial.preproc.warp.reg = [0 0.0001 0.5 0.05 0.2];
matlabbatch{1}.spm.spatial.preproc.warp.affreg = 'mni';
matlabbatch{1}.spm.spatial.preproc.warp.fwhm = 0;
matlabbatch{1}.spm.spatial.preproc.warp.samp = 1;
matlabbatch{1}.spm.spatial.preproc.warp.write = [0 0];
matlabbatch{1}.spm.spatial.preproc.warp.vox = NaN;
matlabbatch{1}.spm.spatial.preproc.warp.bb = [NaN NaN NaN
                                              NaN NaN NaN];                                                                           
 spm_jobman('run', matlabbatch);                                                                     
end   

%% GM BOUNDRY MASKS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% GM boundy mask creation in MP2RAGE space

for sub=1:22

    dir_reg_nc= ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_noc_3/'];
    dir_seg_str= ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/MP2RAGE_WBIC_0pt65_PAT3_PF68_240Hz_UNI/'];

    unix(['fslmaths ' dir_seg_str 'c1clean_uni_reo.nii -thr 0.9 -bin ' dir_seg_str 'c1clean_uni_reo_bin']) %%  mask threshold is 0.9
    unix(['fslmaths ' dir_seg_str 'c1clean_uni_reo_bin -edge ' dir_seg_str 'c1clean_uni_reo_bin_edge'])
    unix(['fslmaths ' dir_seg_str 'c1clean_uni_reo_bin_edge -dilF ' dir_seg_str 'c1clean_uni_reo_bin_edge_dil'])
    unix(['fslmaths ' dir_seg_str 'c1clean_uni_reo_bin_edge_dil -thr 0.4 -bin ' dir_seg_str 'c1clean_uni_reo_bin_edge_dil_thr'])
    unix(['fslmaths ' dir_seg_str 'c1clean_uni_reo_bin_edge_dil_thr -mas ' dir_seg_str 'clean_uni_reo_mask ' dir_seg_str 'c1clean_uni_reo_bin_edge_dil_thr_mask'])

end
%% GM boundy mask creation in the MT-3EPI and 3DEPI spaces

%%Transform boundry mask into the MT-3DEPI space
for sub=1:22
    dir_reg_nc= ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_noc_3/'];
    dir_seg_str= ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/nc_epi3d_v2d_MTw_WholeBrain_0005/'];

    unix(['fslmaths ' dir_seg_str 'c1AP.nii -thr 0.9 -bin ' dir_seg_str 'c1AP_bin']) %%  mask threshold is 0.9
    unix(['fslmaths ' dir_seg_str 'c1AP_bin -edge ' dir_seg_str 'c1AP_bin_edge'])
    unix(['fslmaths ' dir_seg_str 'c1AP_bin_edge -dilF ' dir_seg_str 'c1AP_bin_edge_dil'])
    unix(['fslmaths ' dir_seg_str 'c1AP_bin_edge_dil -thr 0.4 -bin ' dir_seg_str 'c1AP_bin_edge_dil_thr'])
    unix(['fslmaths ' dir_seg_str 'c1AP_bin_edge_dil_thr -mas ' dir_seg_str 'AP_brain_mask ' dir_seg_str 'c1AP_bin_edge_dil_thr_mask'])
end

%%Transform boundry mask into the 3DEPI space
for sub=1:22
    dir1= ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/'];
    dir2= ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/nc_epi3d_v2d_MTw_WholeBrain_0005/'];
    dir_nc = ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_noc_3/'];

     unix(['applywarp --ref=' dir1 'AP --in=' dir2 'c1AP_bin_edge_dil_thr_mask --out=' dir_nc 'c1AP_bin_edge_dil_thr_mask_2fMREPI --premat=' dir_nc 'initialhighres2example_func.mat --interp=nn'])
end


%% ATLAS MASKS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Transform Harvard-Oxford cortical atlas into the MP2RAGE space
for sub=1:22
    dir_atlas = '/export/home/vmalekian/MT_Alice_dataset_20_revision/atlas/';
    dirs=  ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/MP2RAGE_WBIC_0pt65_PAT3_PF68_240Hz_UNI/'];
    diro_nc = ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_noc_3/'];
    unix(['invwarp --ref=' dirs 'clean_uni_reo_brain --warp=' diro_nc 'highres2standard_warp --out=' diro_nc 'highres2standard_warp_inv'])
    unix(['applywarp --ref=' dirs 'clean_uni_reo_brain --in=' dir_atlas 'HarvardOxford-cort-maxprob-thr0-1mm.nii.gz --warp=' diro_nc 'highres2standard_warp_inv --out=' diro_nc 'oxf_thr0_2highres_nonlinear --interp=nn'])

end

%% Transform Harvard-Oxford cortical atlas into the MT-3DEPI and 3DEPI space
for sub=1:22
    dir_atlas = '/export/home/vmalekian/MT_Alice_dataset_20_revision/atlas/';
    dirs=  ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/MP2RAGE_WBIC_0pt65_PAT3_PF68_240Hz_UNI/'];
    diro_nc = ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_noc_3/'];
    diro_b0 = ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_FM_3/'];
    diro_bp = ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_top_3/'];

    dir1= ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/'];
    dir2= ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/nc_epi3d_v2d_MTw_WholeBrain_0005/'];
    % Transform into the MT-3DEPI space
    unix(['convertwarp --ref=/usr/local/fsl/data/standard/MNI152_T1_1mm_brain --premat=' diro_nc 'initial_highres2highres.mat --warp1=' diro_nc 'highres2standard_warp --out=' diro_nc 'example_initialhighres2standard_warp'])
    unix(['invwarp --ref=' dir2 'AP --warp=' diro_nc 'example_initialhighres2standard_warp --out=' diro_nc 'example_initialhighres2standard_warp_inv'])
    unix(['applywarp --ref=' dir2 'AP --in=' dir_atlas 'HarvardOxford-cort-maxprob-thr0-1mm.nii.gz --warp=' diro_nc 'example_initialhighres2standard_warp_inv --out=' diro_nc 'oxf_thr0_2initialhighres --interp=nn'])

    % Transform into the 3DEPI space
    unix(['convertwarp --ref=/usr/local/fsl/data/standard/MNI152_T1_1mm_brain --premat=' diro_nc 'example_func2highres.mat --warp1=' diro_nc 'highres2standard_warp --out=' diro_nc 'example_func2standard_warp'])
    unix(['invwarp --ref=' dir1 'AP --warp=' diro_nc 'example_func2standard_warp --out=' diro_nc 'example_func2standard_warp_inv'])
    unix(['applywarp --ref=' dir1 'AP --in=' dir_atlas 'HarvardOxford-cort-maxprob-thr0-1mm.nii.gz --warp=' diro_nc 'example_func2standard_warp_inv --out=' diro_nc 'oxf_thr0_2func --interp=nn'])
 
end
%% ROI MASKS  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Two spherical ROIs (vmPFC & dmPFC) in MNI space

dirt ='/export/home/vmalekian/MT_Alice_dataset_20_revision/atlas/';
unix(['fslmaths /usr/local/fsl/data/standard/MNI152_T1_1mm  -mul 0 -add 1 -roi 91 1 174 1 64 1 0 1 ' dirt 'vmPFCpoint -odt float']) %%vmPFC ROI1
unix(['fslmaths ' dirt 'vmPFCpoint -kernel sphere 26 -fmean ' dirt 'vmPFCsphere -odt float'])
unix(['fslmaths ' dirt 'vmPFCsphere.nii.gz -thr 0.0000000001 -bin ' dirt 'vmPFCsphere_bin.nii.gz'])


unix(['fslmaths /usr/local/fsl/data/standard/MNI152_T1_1mm  -mul 0 -add 1 -roi 91 1 179 1 102 1 0 1 ' dirt 'dmPFCpoint -odt float'])%%dmPFC ROI2
unix(['fslmaths ' dirt 'dmPFCpoint -kernel sphere 21 -fmean ' dirt 'dmPFCsphere -odt float'])
unix(['fslmaths ' dirt 'dmPFCsphere.nii.gz -thr 0.0000000001 -bin ' dirt 'dmPFCsphere_bin.nii.gz'])

%% MNI ROIS into MP2RAGE space
for sub=1:22
    dir_atlas = '/export/home/vmalekian/MT_Alice_dataset_20_revision/atlas/';
    dirs=  ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/MP2RAGE_WBIC_0pt65_PAT3_PF68_240Hz_UNI/'];
    diro_nc = ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_noc_3/'];
    unix(['applywarp --ref=' dirs 'clean_uni_reo_brain --in=' dir_atlas 'dmPFCsphere_bin --warp=' diro_nc 'highres2standard_warp_inv --out=' diro_nc 'dmPFC_highres --interp=nn'])
    unix(['applywarp --ref=' dirs 'clean_uni_reo_brain --in=' dir_atlas 'vmPFCsphere_bin --warp=' diro_nc 'highres2standard_warp_inv --out=' diro_nc 'vmPFC_highres --interp=nn'])
end

%% MNI ROIS into MT-3DEPI and 3DEPI spaces
for sub=1:22
    dir_atlas = '/export/home/vmalekian/MT_Alice_dataset_20_revision/atlas/';
    dirs=  ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/MP2RAGE_WBIC_0pt65_PAT3_PF68_240Hz_UNI/'];
    diro_nc = ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_noc_3/'];
    
    dir1= ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/'];
    dir2= ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/nc_epi3d_v2d_MTw_WholeBrain_0005/'];
    
    %%Two spherical ROIs into MT-3DEPI   
    unix(['applywarp --ref=' dir2 'AP --in=' dir_atlas 'dmPFCsphere_bin --warp=' diro_nc 'example_initialhighres2standard_warp_inv --out=' diro_nc 'dmPFC_2initialhighres --interp=nn'])
    unix(['applywarp --ref=' dir2 'AP --in=' dir_atlas 'vmPFCsphere_bin --warp=' diro_nc 'example_initialhighres2standard_warp_inv --out=' diro_nc 'vmPFC_2initialhighres --interp=nn'])
    %%Two spherical ROIs into 3DEPI   
    unix(['applywarp --ref=' dir1 'AP --in=' dir_atlas 'dmPFCsphere_bin --warp=' diro_nc 'example_func2standard_warp_inv --out=' diro_nc 'dmPFC_2func --interp=nn'])
    unix(['applywarp --ref=' dir1 'AP --in=' dir_atlas 'vmPFCsphere_bin --warp=' diro_nc 'example_func2standard_warp_inv --out=' diro_nc 'vmPFC_2func --interp=nn'])
end


%% Group-level MP2RAGE GM in MNI space (to use for display in relative CR maps) 
for sub=1:22
    
        dir_struct=  ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/MP2RAGE_WBIC_0pt65_PAT3_PF68_240Hz_UNI/'];
        diro_nc= ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_noc_3/'];
        dir_std= '/export/data/vmalekian/MT_Alice_dataset_20_revision/standard/';
        unix(['applywarp --ref=' (diro_nc) 'standard --in=' (dir_struct) 'c1clean_uni_reo.nii --warp=' (dir_struct) 'highres2standard_warp --out=' (dir_struct) 'c1clean_uni_reo2standard'])  
        unix(['cp ' (dir_struct) 'c1clean_uni_reo2standard.nii.gz ' (dir_std) 'c1clean_uni_reo2standard_no_' sprintf('%02d',sub) '.nii.gz'])
end

unix(['fslmerge -t ' (dir_std) 'c1clean_uni_reo2standard_no_total ' (dir_std) 'c1clean_uni_reo2standard_no_*.nii.gz'])
unix(['fslmaths ' (dir_std) 'c1clean_uni_reo2standard_no_total -Tmean ' (dir_std) 'c1clean_uni_reo2standard_no_total_mean'])



