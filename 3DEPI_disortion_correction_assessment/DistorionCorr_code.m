clc
clear 
close all

sub =22;


%% Topup LSR and Jacobian
for k=1:sub
    dir1=  ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',k) '/fMRI0003/'];

    unix(['fslmerge -t ' dir1 'AP_PA ' dir1 'AP ' dir1 'PA'])

    unix(['topup --imain=' dir1 'AP_PA --datain=param.txt --config=b02b0.cnf --out=' dir1 'topup_AP_PA --fout=' dir1 'topup_field --iout=' dir1 'topup_unwarp --jacout=' dir1 'jac --rbmout=' dir1 'xfm --dfout=' dir1 'topup_warpfield'])
    unix(['applytopup --imain=' (dir1) 'AP --topup=' (dir1) 'topup_AP_PA --datain=param.txt --inindex=1 --out=' (dir1) 'AP_dc_def_jac --method=jac']) %%Jacobian
    unix(['applytopup --imain=' (dir1) 'AP,' (dir1) 'PA --topup=' (dir1) 'topup_AP_PA --datain=param.txt --inindex=1,2 --out=' dir1 'AP_dc_def']) %%LSR
    unix(['fslmaths ' dir1 'AP_dc_def -thr 0 ' dir1 'AP_dc_def'])
    unix(['fslmaths ' dir1 'AP_dc_def_jac -thr 0 ' dir1 'AP_dc_def_jac'])
   
end

%% Field map
for k=1:sub
    dir_magn=  ['/export/data/vmalekian/MT_Alice_dataset_22sub/' sprintf('%02d',k) '/aa_b0_Sag_WE_20mm_0006/'];
    dir_phase= ['/export/data/vmalekian/MT_Alice_dataset_22sub/' sprintf('%02d',k) '/aa_b0_Sag_WE_20mm_0007/'];

    unix(['fslreorient2std ' (dir_magn) 'magn1.nii ' (dir_magn) 'magn1_reo.nii.gz'])
    unix(['fslreorient2std ' (dir_phase) 'phase.nii ' (dir_phase) 'phase_reo.nii.gz'])
    unix(['bet ' (dir_magn) 'magn1_reo.nii.gz ' (dir_magn) 'magn1_reo_brain1.nii.gz -B -f 0.5 -g 0 -m'])
    unix(['fslmaths ' (dir_magn) 'magn1_reo_brain1.nii.gz -ero ' (dir_magn) 'magn1_reo_brain.nii.gz'])
    unix(['fsl_prepare_fieldmap SIEMENS ' (dir_phase) 'phase_reo.nii.gz ' (dir_magn) 'magn1_reo_brain.nii.gz ' (dir_phase) 'fmap_rads.nii.gz 3.08'])
    unix(['fslmaths ' (dir_phase) 'fmap_rads.nii.gz -div 6.28 ' (dir_phase) 'fmap_hz.nii.gz'])
end

%% MP2RAGE GM Segmentation unisg SPM

%%Reorient
for k=1:sub
    dir_struct=  ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',k) '/MP2RAGE_WBIC_0pt65_PAT3_PF68_240Hz_UNI/'];
    dir_struct1=  ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',k) '/MP2RAGE_WBIC_0pt65_PAT3_PF68_240Hz_INV1/'];
    dir_struct2=  ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',k) '/MP2RAGE_WBIC_0pt65_PAT3_PF68_240Hz_INV2/'];

    unix(['fslreorient2std ' (dir_struct) 'uni.nii ' (dir_struct) 'uni_reo.nii.gz'])
    unix(['gunzip ' (dir_struct) 'uni_reo.nii.gz'])
    unix(['fslreorient2std ' (dir_struct1) 'inv1.nii ' (dir_struct1) 'inv1_reo.nii.gz'])
    unix(['gunzip ' (dir_struct1) 'inv1_reo.nii.gz'])
    unix(['fslreorient2std ' (dir_struct2) 'inv2.nii ' (dir_struct2) 'inv2_reo.nii.gz'])
    unix(['gunzip ' (dir_struct2) 'inv2_reo.nii.gz'])
    
end
%%SPM MP2RAGE tissue segmentation
for k=1:sub
    
matlabbatch{1}.spm.tools.mp2rage.rmbg.INV1 = { ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',k) '/MP2RAGE_WBIC_0pt65_PAT3_PF68_240Hz_INV1/inv1_reo.nii,1']};
matlabbatch{1}.spm.tools.mp2rage.rmbg.INV2 = { ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',k) '/MP2RAGE_WBIC_0pt65_PAT3_PF68_240Hz_INV2/inv2_reo.nii,1']};
matlabbatch{1}.spm.tools.mp2rage.rmbg.UNI =  { ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',k) '/MP2RAGE_WBIC_0pt65_PAT3_PF68_240Hz_UNI/uni_reo.nii,1']};

matlabbatch{1}.spm.tools.mp2rage.rmbg.regularization = 1;
matlabbatch{1}.spm.tools.mp2rage.rmbg.output.prefix = 'clean_';
matlabbatch{1}.spm.tools.mp2rage.rmbg.show = 'yes';
matlabbatch{2}.spm.tools.MRI.MRTool_preproc.MRTool_brain.res_dir = {['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',k) '/MP2RAGE_WBIC_0pt65_PAT3_PF68_240Hz_UNI']};

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

%%Brain extraction image calculator(i1+i2+i3)>0
for k=1:sub
    dir_struct=  ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',k) '/MP2RAGE_WBIC_0pt65_PAT3_PF68_240Hz_UNI/'];
    unix(['fslmaths ' (dir_struct) 'c1clean_uni_reo.nii -add ' (dir_struct) 'c2clean_uni_reo.nii -add ' (dir_struct) 'c3clean_uni_reo.nii -bin ' (dir_struct) 'clean_uni_reo_masc'])
    unix(['fslmaths ' (dir_struct) 'clean_uni_reo_masc.nii -dilM ' (dir_struct) 'clean_uni_reo_masc_dil.nii'])
    unix(['fslmaths ' (dir_struct) 'clean_uni_reo_masc_dil.nii -fillh ' (dir_struct) 'clean_uni_reo_masc_dil_fill.nii'])
    unix(['fslmaths ' (dir_struct) 'clean_uni_reo_masc_dil_fill.nii -ero ' (dir_struct) 'clean_uni_reo_masck.nii'])
    unix(['fslmaths ' (dir_struct) 'clean_uni_reo -mul ' (dir_struct) 'clean_uni_reo_masck '  (dir_struct) 'clean_uni_reo_brain'])
end


%% No correction  ; Reversed-PE and B0 field mapping 


%%Creating folders for analysis results

for sub=1:22

    unix(['mkdir -p /export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_FM_3'])
    unix(['mkdir -p /export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_noc_3'])
%     unix(['mkdir -p /export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_topjac_3'])
    unix(['mkdir -p /export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_top_3'])

end

%%No correction + MP2RAGE to NMI(1mm) transformation
for sub=1:22
    
    dirs=  ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/MP2RAGE_WBIC_0pt65_PAT3_PF68_240Hz_UNI/'];
    dir1= ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/'];
    dir2= ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/nc_epi3d_v2d_MTw_WholeBrain_0005/'];
    diro = ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_noc_3/'];
    dir_reg_no = ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_no.feat/reg/'];

    unix(['fslmaths ' dirs 'c2clean_uni_reo.nii.gz -thr 0.5 -bin ' dirs 'c2clean_uni_reo_wmseg.nii.gz'])
    unix(['flirt -in ' dir1 'AP -ref ' dir2 'AP -out ' diro 'example_func2initial_highres -omat ' diro 'example_func2initial_highres.mat -cost corratio -dof 6 -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -interp spline']) 
    unix(['epi_reg --epi=' dir2 'AP --t1=' dirs 'clean_uni_reo --t1brain=' dirs 'clean_uni_reo_brain --out=' diro 'initial_highres2highres --wmseg=' dirs 'c2clean_uni_reo_wmseg.nii.gz']) %%initial_highres2highres (no correction)
    unix(['convert_xfm -omat ' diro 'example_func2highres.mat -concat ' diro 'initial_highres2highres.mat ' diro 'example_func2initial_highres.mat'])
    unix(['flirt -ref ' dirs 'clean_uni_reo_brain -in ' dir1 'AP -out ' diro 'example_func2highres -applyxfm -init ' diro 'example_func2highres.mat -interp spline']) %%fMRI2highres (no correction)
    unix(['flirt -in ' dirs 'clean_uni_reo_brain -ref /usr/local/fsl/data/standard/MNI152_T1_1mm_brain -out ' diro 'highres2standard -omat ' diro 'highres2standard.mat -cost corratio -dof 12 -searchrx -180 180 -searchry -180 180 -searchrz -180 180 -interp spline'])
%     unix(['fnirt --iout=' diro 'highres2standard_head --in=' dirs 'clean_uni_reo --aff=' diro 'highres2standard.mat --cout=' diro 'highres2standard_warp --iout=' diro 'highres2standard --jout=' diro 'highres2std_jac --config=/export/data/vmalekian/MT_Alice_dataset_20_revision/T1_2_MNI152_1mm --ref=/usr/local/fsl/data/standard/MNI152_T1_1mm  --warpres=10,10,10 --interp=spline'])
    unix(['fnirt --in=' dirs 'clean_uni_reo.nii --ref=' dir_reg_no 'standard.nii.gz --aff=' diro 'highres2standard.mat --cout=' diro 'highres2standard_warp'])
    unix(['applywarp -i ' dirs 'clean_uni_reo_brain -r /usr/local/fsl/data/standard/MNI152_T1_1mm_brain -o ' diro 'highres2standard -w ' diro 'highres2standard_warp --interp=spline'])
    unix(['convert_xfm -inverse -omat ' diro 'standard2highres.mat ' diro 'highres2standard.mat'])
    unix(['convert_xfm -omat ' diro 'example_func2highres.mat -concat ' diro 'initial_highres2highres.mat ' diro 'example_func2initial_highres.mat'])
    unix(['flirt -ref ' dirs 'clean_uni_reo_brain -in ' dir1 'AP -out ' diro 'example_func2highres -applyxfm -init ' diro 'example_func2highres.mat -interp spline'])
    unix(['convert_xfm -inverse -omat ' diro 'highres2example_func.mat ' diro 'example_func2highres.mat'])
    unix(['convert_xfm -omat ' diro 'example_func2standard.mat -concat ' diro 'highres2standard.mat ' diro 'example_func2highres.mat'])
    unix(['convertwarp --ref=/usr/local/fsl/data/standard/MNI152_T1_1mm_brain --premat=' diro 'example_func2highres.mat --warp1=' diro 'highres2standard_warp --out=' diro 'example_func2standard_warp'])
    unix(['applywarp --ref=/usr/local/fsl/data/standard/MNI152_T1_1mm_brain --in=' dir1 'AP --out=' diro 'example_func2standard --warp=' diro 'example_func2standard_warp --interp=spline'])
    unix(['convert_xfm -inverse -omat ' diro 'standard2example_func.mat ' diro 'example_func2standard.mat'])
end


%%Reversed-PE LS and JAC
for sub=1:22
    
    dirs=  ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/MP2RAGE_WBIC_0pt65_PAT3_PF68_240Hz_UNI/'];
    dir1= ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/'];
    dir2= ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/nc_epi3d_v2d_MTw_WholeBrain_0005/'];
    diro = ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_top_3/'];
    dir_nc = ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_noc_3/'];

    unix(['flirt -in ' dir1 'AP_dc_def -ref ' dir2 'AP_dc_def -out ' diro 'example_func2initial_highres -omat ' diro 'example_func2initial_highres.mat -cost corratio -dof 6 -searchrx -90 90 -searchry -90 90 -searchrz -90 90 -interp spline'])
    unix(['epi_reg --epi=' dir2 'AP_dc_def --t1=' dirs 'clean_uni_reo --t1brain=' dirs 'clean_uni_reo_brain --out=' diro 'initial_highres2highres --wmseg=' dirs 'c2clean_uni_reo_wmseg.nii.gz']) %% initial_highres2highres (LSR)
    unix(['convert_xfm -omat ' diro 'example_func2highres.mat -concat ' diro 'initial_highres2highres.mat ' diro 'example_func2initial_highres.mat'])   
    unix(['flirt -ref ' dirs 'clean_uni_reo_brain -in ' dir1 'AP_dc_def_jac -out ' diro 'example_func2highres -applyxfm -init ' diro 'example_func2highres.mat -interp spline']) %% fMRI2highres (Jac)
    unix(['convert_xfm -omat ' diro 'example_func2highres.mat -concat ' diro 'initial_highres2highres.mat ' diro 'example_func2initial_highres.mat'])
    unix(['flirt -ref ' dirs 'clean_uni_reo_brain -in ' dir2 'AP_dc_def_jac -out ' diro 'initial_highres2highres_jac -applyxfm -init ' diro 'initial_highres2highres.mat -interp spline'])  %% initial_highres2highres (Jac)
    
   
%     unix(['applywarp -r ' dir1 'AP -i ' dir1 'AP --warp=' dir1 'warpfield_01 -o ' dir1 'AP_warpped --interp=spline'])
%     unix(['applywarp -r ' dir2 'AP -i ' dir2 'AP --warp=' dir2 'warpfield_01 -o ' dir2 'AP_warpped --interp=spline']) 
%     unix(['convertwarp --ref=' dirs 'clean_uni_reo_brain --postmat=' diro 'initial_highres2highres.mat --warp1=' dir2 'warpfield_01 --out=' diro 'initial_highres2highres_warp --relout'])
%     unix(['convertwarp --ref=' dir2 'AP_dc_def --postmat=' diro 'example_func2initial_highres.mat --warp1=' dir1 'warpfield_01 --out=' diro 'func2initial_highres_warp1 --relout'])
%     unix(['convertwarp --ref=' dirs 'clean_uni_reo_brain --postmat=' dir_nc 'example_func2initial_highres.mat --warp1=' dir1 'warpfield_01 --out=' diro 'func2initial_highres_warp --relout'])
%     unix(['convertwarp --ref=' dirs 'clean_uni_reo_brain  --warp2=' diro 'initial_highres2highres_warp --warp1=' diro 'func2initial_highres_warp --out=' diro 'func2highres_warp --relout'])    
%     unix(['applywarp -r ' dir2 'AP -i ' dir1 'AP --warp=' diro 'func2initial_highres_warp1 -o ' diro 'example_func2initalhighres_2 --interp=spline'])
%     unix(['applywarp -r ' dirs 'clean_uni_reo_brain -i ' dir2 'AP --warp=' diro 'initial_highres2highres_warp -o ' diro 'initalhighres2highres_1 --interp=spline'])

end


%%B0_field_mapping
for sub=1:22
    
    dirs=  ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/MP2RAGE_WBIC_0pt65_PAT3_PF68_240Hz_UNI/'];
    dir1= ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/'];
    dir2= ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/nc_epi3d_v2d_MTw_WholeBrain_0005/'];
    diro = ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_FM_3/'];
    

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



%% Native space

for sub=1:22
    
    dirs=  ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/MP2RAGE_WBIC_0pt65_PAT3_PF68_240Hz_UNI/'];
    dir2= ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/nc_epi3d_v2d_MTw_WholeBrain_0005/'];
    dir_nc = ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_noc_3/'];
    dir_b0 = ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_FM_3/'];
    dir_bp = ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_top_3/'];

    %%NC
     unix(['convert_xfm -inverse -omat ' dir_nc 'highres2initial_highres.mat ' dir_nc 'initial_highres2highres.mat'])
     unix(['applywarp --ref=' dir2 'AP --in=' dirs 'clean_uni_reo_brain --out=' dir_nc 'clean_uni_reo_brain2EPI --premat=' dir_nc 'highres2initial_highres.mat --interp=spline'])
     unix(['applywarp --ref=' dir2 'AP --in=' dirs 'c1clean_uni_reo --out=' dir_nc 'c1brain2EPI --premat=' dir_nc 'highres2initial_highres.mat --interp=spline'])
     
    %%B0
    unix(['invwarp -w ' dir_b0 'initial_highres2highres_warp -o ' dir_b0 'highres2intial_highres_warp -r ' dir2 'AP'])
    unix(['applywarp --ref=' dir2 'AP --in=' dirs 'clean_uni_reo_brain --out=' dir_b0 'clean_uni_reo_brain2EPI --warp=' dir_b0 'highres2intial_highres_warp --interp=spline'])
    unix(['applywarp --ref=' dir2 'AP --in=' dirs 'c1clean_uni_reo --out=' dir_b0 'c1brain2EPI --warp=' dir_b0 'highres2intial_highres_warp --interp=spline'])
     
    %%BP
    unix(['convertwarp --ref=' dirs 'clean_uni_reo_brain --postmat=' dir_bp 'initial_highres2highres.mat --warp1=' dir2 'warpfield_01 --out=' dir_bp 'initial_highres2highres_warp --relout'])
    unix(['invwarp -w ' dir_bp 'initial_highres2highres_warp -o ' dir_bp 'highres2intial_highres_warp -r ' dir2 'AP'])
    unix(['applywarp --ref=' dir2 'AP --in=' dirs 'clean_uni_reo_brain --out=' dir_bp 'clean_uni_reo_brain2EPI --warp=' dir_bp 'highres2intial_highres_warp --interp=spline'])
    unix(['applywarp --ref=' dir2 'AP --in=' dirs 'c1clean_uni_reo --out=' dir_bp 'c1brain2EPI --warp=' dir_bp 'highres2intial_highres_warp --interp=spline'])
     

end


%% Transform Harvard-Oxford cortical Atlas 2 MP2RAGE space

for sub=1:22
    dir_atlas = '/export/data/vmalekian/MT_Alice_dataset_20_revision/atlas/';
    dirs=  ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/MP2RAGE_WBIC_0pt65_PAT3_PF68_240Hz_UNI/'];
    diro_nc = ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_noc_3/'];
    unix(['invwarp --ref=' dirs 'clean_uni_reo_brain --warp=' diro_nc 'highres2standard_warp --out=' diro_nc 'highres2standard_warp_inv'])
    unix(['applywarp --ref=' dirs 'clean_uni_reo_brain --in=' dir_atlas 'HarvardOxford-cort-maxprob-thr0-1mm.nii.gz --warp=' diro_nc 'highres2standard_warp_inv --out=' diro_nc 'oxf_thr0_2highres_nonlinear --interp=nn'])
%     unix(['applywarp --ref=' dirs 'clean_uni_reo_brain --in=' dir_atlas 'HarvardOxford-sub-maxprob-thr0-1mm.nii.gz --warp=' diro_nc 'highres2standard_warp_inv --out=' diro_nc 'oxf_sub_thr0_2highres_nonlinear --interp=nn'])
end

%% MT-3DEPI (initial_highres2highres) GM segmentation using SPM 


%%SPM format preparation
for k=1:22
    
 
    dir_reg_nc= ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',k) '/fMRI0003/AP_noc_3/'];
    dir_reg_b0= ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',k) '/fMRI0003/AP_FM_3/'];
    dir_reg_bp_jclsr= ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',k) '/fMRI0003/AP_top_3/'];
    
    unix(['fslmaths ' (dir_reg_nc) 'initial_highres2highres -thr 0 ' (dir_reg_nc) 'initial_highres2highres'])
    unix(['fslmaths ' (dir_reg_b0) 'initial_highres2highres -thr 0 ' (dir_reg_b0) 'initial_highres2highres'])
    unix(['fslmaths ' (dir_reg_bp_jclsr) 'initial_highres2highres -thr 0 ' (dir_reg_bp_jclsr) 'initial_highres2highres'])
    unix(['fslmaths ' (dir_reg_bp_jclsr) 'initial_highres2highres_jac -thr 0 ' (dir_reg_bp_jclsr) 'initial_highres2highres_jac'])
    
    unix(['gunzip ' (dir_reg_nc) 'initial_highres2highres.nii.gz'])
    unix(['gunzip ' (dir_reg_b0) 'initial_highres2highres.nii.gz'])
    unix(['gunzip ' (dir_reg_bp_jclsr) 'initial_highres2highres.nii.gz'])
    unix(['gunzip ' (dir_reg_bp_jclsr) 'initial_highres2highres_jac.nii.gz'])
    
 
end

%%SPM segmentation
 for k=1:22

    matlabbatch{1}.spm.spatial.preproc.channel.vols = {

          ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',k) '/fMRI0003/AP_noc_3/initial_highres2highres.nii,1'] 
          ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',k) '/fMRI0003/AP_FM_3/initial_highres2highres.nii,1'] 
          ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',k) '/fMRI0003/AP_top_3/initial_highres2highres.nii,1'] 
          ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',k) '/fMRI0003/AP_top_3/initial_highres2highres_jac.nii,1'] 

}

    matlabbatch{1}.spm.spatial.preproc.channel.biasreg = 0.01;
    matlabbatch{1}.spm.spatial.preproc.channel.biasfwhm = 30;
    matlabbatch{1}.spm.spatial.preproc.channel.write = [1 1];
    matlabbatch{1}.spm.spatial.preproc.tissue(1).tpm = {'/export/data/vmalekian/spm12/tpm/TPM.nii,1'};
    matlabbatch{1}.spm.spatial.preproc.tissue(1).ngaus = 1;
    matlabbatch{1}.spm.spatial.preproc.tissue(1).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(1).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(2).tpm = {'/export/data/vmalekian/spm12/tpm/TPM.nii,2'};
    matlabbatch{1}.spm.spatial.preproc.tissue(2).ngaus = 1;
    matlabbatch{1}.spm.spatial.preproc.tissue(2).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(2).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(3).tpm = {'/export/data/vmalekian/spm12/tpm/TPM.nii,3'};
    matlabbatch{1}.spm.spatial.preproc.tissue(3).ngaus = 2;
    matlabbatch{1}.spm.spatial.preproc.tissue(3).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(3).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).tpm = {'/export/data/vmalekian/spm12/tpm/TPM.nii,4'};
    matlabbatch{1}.spm.spatial.preproc.tissue(4).ngaus = 3;
    matlabbatch{1}.spm.spatial.preproc.tissue(4).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(4).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).tpm = {'/export/data/vmalekian/spm12/tpm/TPM.nii,5'};
    matlabbatch{1}.spm.spatial.preproc.tissue(5).ngaus = 4;
    matlabbatch{1}.spm.spatial.preproc.tissue(5).native = [1 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(5).warped = [0 0];
    matlabbatch{1}.spm.spatial.preproc.tissue(6).tpm = {'/export/data/vmalekian/spm12/tpm/TPM.nii,6'};
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
