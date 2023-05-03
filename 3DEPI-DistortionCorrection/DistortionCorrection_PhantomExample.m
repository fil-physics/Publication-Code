%% 3DEPI distortion correction piplines at 7T for phantom data
% No-correction, B0 field mapping and reversed-PE (fsl_topup) techniques
% V Malekian, FIL Physics

clc;clear;close all;
addpath('/export/home/vmalekian/MI_Distortion/NIfTI_20140122/')
addpath('/export/data/vmalekian/MI_Distortion/mi/')
addpath('/export/data/vmalekian/MI_Distortion/NIfTI_20140122/')
addpath('/export/data/vmalekian/MI_Distortion/InfoTheory/')

dir1=  '/export/data/vmalekian/test/phantom/phantom_20220412/nc_epi3d_v2u_0pt92_WB_Seg8_AP_PA_TE25_160s_0012/';
dir2=  '/export/data/vmalekian/test/phantom/phantom_20220412/nc_epi3d_v2u_0pt92_WB_Seg4_AP_PA_TE25_160s_0013/';

dir_magn=  '/export/data/vmalekian/test/phantom/phantom_20220412/aa_b0_Sag_WE_08mm_0014/';
dir_phase= '/export/data/vmalekian/test/phantom/phantom_20220412/aa_b0_Sag_WE_08mmB0_Map_0015/';

dir_magn_low=  '/export/data/vmalekian/test/phantom/phantom_20220412/aa_b0_Sag_WE_20mm_0016/';
dir_phase_low= '/export/data/vmalekian/test/phantom/phantom_20220412/aa_b0_Sag_WE_20mmB0_Map_0017/';


%% Reversed-PE approach for two-fold and non-segmented data

%%Two-fold segmented EPI
unix(['fslmerge -t ' (dir1) 'AP_PA_roi1 ' (dir1) 'AP_roi1.nii ' (dir1) 'PA_roi1.nii'])
unix(['topup --imain=' (dir1) 'AP_PA_roi1 --datain=param_seg8_inv.txt --config=b02b0.cnf --out=' (dir1) 'topup_AP_PA_b0_def_roi1 --fout=' (dir1) 'my_field_def_roi1 --iout=' (dir1) 'my_unwarp_def_roi1 --jacout=' dir1 'jac --rbmout=' dir1 'xfm --dfout=' dir1 'warpfield --logout=' dir1 'new_topup.log'])
unix(['applytopup --imain=' (dir1) 'AP_roi1 --topup=' (dir1) 'topup_AP_PA_b0_def_roi1 --datain=param_seg8_inv.txt --inindex=1 --out=' (dir1) 'AP_dc_def_roi1_jac --method=jac'])
unix(['applytopup --imain=' (dir1) 'PA_roi1 --topup=' (dir1) 'topup_AP_PA_b0_def_roi1 --datain=param_seg8_inv.txt --inindex=2 --out=' (dir1) 'PA_dc_def_roi1_jac --method=jac'])
unix(['fslmaths ' (dir1) 'AP_dc_def_roi1_jac -thr 0 ' (dir1) 'AP_dc_def_roi1_jac'])
unix(['fslmaths ' (dir1) 'PA_dc_def_roi1_jac -thr 0 ' (dir1) 'PA_dc_def_roi1_jac'])

unix(['applywarp -r ' dir1 'AP_roi1 -i ' dir1 'AP_roi1 --warp=' dir1 'warpfield_01 -o ' dir1 'AP_roi1_warp01'])
unix(['fslmaths ' dir1 'AP_roi1_warp01 -mul ' dir1 'jac_01.nii.gz ' dir1 'AP_roi1_warp01_jac'])
unix(['applywarp -r ' dir1 'PA_roi1 -i ' dir1 'PA_roi1 --warp=' dir1 'warpfield_02 -o ' dir1 'PA_roi1_warp02'])
unix(['fslmaths ' dir1 'PA_roi1_warp02 -mul ' dir1 'jac_02.nii.gz ' dir1 'PA_roi1_warp02_jac'])

%%Non-segmented EPI
unix(['fslmerge -t ' (dir2) 'AP_PA_roi1 ' (dir2) 'AP_roi1.nii ' (dir2) 'PA_roi1.nii'])
unix(['topup --imain=' (dir2) 'AP_PA_roi1 --datain=param_seg4_inv.txt --config=b02b0.cnf --out=' (dir2) 'topup_AP_PA_b0_def_roi1 --fout=' (dir2) 'my_field_def_roi1 --iout=' (dir2) 'my_unwarp_def_roi1 --jacout=' dir2 'jac --rbmout=' dir2 'xfm --dfout=' dir2 'warpfield --logout=' dir2 'new_topup.log'])
unix(['applytopup --imain=' (dir2) 'AP_roi1 --topup=' (dir2) 'topup_AP_PA_b0_def_roi1 --datain=param_seg4_inv.txt --inindex=1 --out=' (dir2) 'AP_dc_def_roi1_jac --method=jac'])
unix(['applytopup --imain=' (dir2) 'PA_roi1 --topup=' (dir2) 'topup_AP_PA_b0_def_roi1 --datain=param_seg4_inv.txt --inindex=2 --out=' (dir2) 'PA_dc_def_roi1_jac --method=jac'])
unix(['fslmaths ' (dir2) 'AP_dc_def_roi1_jac -thr 0 ' (dir2) 'AP_dc_def_roi1_jac'])
unix(['fslmaths ' (dir2) 'PA_dc_def_roi1_jac -thr 0 ' (dir2) 'PA_dc_def_roi1_jac'])

unix(['applywarp -r ' dir2 'AP_roi1 -i ' dir2 'AP_roi1 --warp=' dir2 'warpfield_01 -o ' dir2 'AP_roi1_warp01'])
unix(['fslmaths ' dir2 'AP_roi1_warp01 -mul ' dir2 'jac_01.nii.gz ' dir2 'AP_roi1_warp01_jac'])
unix(['applywarp -r ' dir2 'PA_roi1 -i ' dir2 'PA_roi1 --warp=' dir2 'warpfield_02 -o ' dir2 'PA_roi1_warp02'])
unix(['fslmaths ' dir2 'PA_roi1_warp02 -mul ' dir2 'jac_02.nii.gz ' dir2 'PA_roi1_warp02_jac'])

%%Deformation fields
unix(['fslmaths ' (dir1) 'my_field_def_roi1 -mul ' (dir_magn) 'magn1_roi1_thr_bin '  (dir1) 'my_field_def_roi1_masked'])
unix(['fslmaths ' (dir2) 'my_field_def_roi1 -mul ' (dir_magn) 'magn1_roi1_thr_bin '  (dir2) 'my_field_def_roi1_masked'])

%% High-resolution B0 field map (0.8*0.8*0.8 mm3)

unix(['fslmaths ' (dir_magn) 'magn1_roi1 -thr 200 ' (dir_magn) 'magn1_roi1_thr']) % masking with threshold of 200
unix(['fslmaths ' (dir_magn) 'magn1_roi1_thr -bin ' (dir_magn) 'magn1_roi1_thr_bin'])
unix(['fslmaths ' (dir_magn) 'magn1_roi1 -mul ' (dir_magn) 'magn1_roi1_thr_bin '  (dir_magn) 'magn1_roi1_brain'])
unix(['fslmaths ' (dir_magn) 'magn1_roi1_brain.nii.gz -ero ' (dir_magn) 'magn1_roi1_brain_erd.nii.gz'])

%Edge detector
unix(['fslmaths ' (dir_magn) 'magn1_roi1_thr -bin ' (dir_magn) 'magn1_roi1_thr_bin'])
unix(['fslmaths ' (dir_magn) 'magn1_roi1_thr_bin -edge -bin -mas ' (dir_magn) 'magn1_roi1_thr_bin ' (dir_magn) 'magn1_roi1_thr_bin_edge'])


unix(['fsl_prepare_fieldmap SIEMENS ' (dir_phase) 'phase_roi1.nii.gz ' (dir_magn) 'magn1_roi1_brain.nii.gz ' (dir_phase) 'fmap_roi1_rads.nii.gz 4.3'])
unix(['fslmaths ' (dir_phase) 'fmap_roi1_rads.nii.gz -div 6.28 ' (dir_phase) 'fmap_roi1_hz.nii.gz'])

unix(['fsl_prepare_fieldmap SIEMENS ' (dir_phase) 'phase_roi1.nii.gz ' (dir_magn) 'magn1_roi1_brain_erd.nii.gz ' (dir_phase) 'fmap_roi1_rads_erd.nii.gz 4.3'])
unix(['fslmaths ' (dir_phase) 'fmap_roi1_rads_erd.nii.gz -div 6.28 ' (dir_phase) 'fmap_roi1_hz_erd.nii.gz'])

%%Erosion, regularisation and intensity correction effects
%%Non-eroded
unix(['fugue -i ' (dir1) 'AP_roi1 --loadfmap=' (dir_phase) 'fmap_roi1_rads --dwell=0.0001525 --asym=0.0043 --despike --despikethreshold=2.1 -u  ' (dir1) 'AP_roi1_FM_spike'])
unix(['fugue -i ' (dir2) 'AP_roi1 --loadfmap=' (dir_phase) 'fmap_roi1_rads --dwell=0.0003050 --asym=0.0043 --despike --despikethreshold=2.1 -u  ' (dir2) 'AP_roi1_FM_spike'])

unix(['fugue -i ' (dir1) 'PA_roi1 --loadfmap=' (dir_phase) 'fmap_roi1_rads --dwell=0.0001525 --asym=0.0043 --despike --despikethreshold=2.1 --unwarpdir=y- -u  ' (dir1) 'PA_roi1_FM_spike'])
unix(['fugue -i ' (dir2) 'PA_roi1 --loadfmap=' (dir_phase) 'fmap_roi1_rads --dwell=0.0003050 --asym=0.0043 --despike --despikethreshold=2.1 --unwarpdir=y- -u  ' (dir2) 'PA_roi1_FM_spike'])

%%Eroded
unix(['fugue -i ' (dir1) 'AP_roi1 --loadfmap=' (dir_phase) 'fmap_roi1_rads_erd --dwell=0.0001525 --asym=0.0043 --despike --despikethreshold=2.1 -u  ' (dir1) 'AP_roi1_FM_erd_spike'])
unix(['fugue -i ' (dir2) 'AP_roi1 --loadfmap=' (dir_phase) 'fmap_roi1_rads_erd --dwell=0.0003050 --asym=0.0043 --despike --despikethreshold=2.1 -u  ' (dir2) 'AP_roi1_FM_erd_spike'])

unix(['fugue -i ' (dir1) 'PA_roi1 --loadfmap=' (dir_phase) 'fmap_roi1_rads_erd --dwell=0.0001525 --asym=0.0043 --despike --despikethreshold=2.1 --unwarpdir=y- -u  ' (dir1) 'PA_roi1_FM_erd_spike'])
unix(['fugue -i ' (dir2) 'PA_roi1 --loadfmap=' (dir_phase) 'fmap_roi1_rads_erd --dwell=0.0003050 --asym=0.0043 --despike --despikethreshold=2.1 --unwarpdir=y- -u  ' (dir2) 'PA_roi1_FM_erd_spike'])

%%Using polynomial kernel
unix(['fugue -i ' (dir1) 'AP_roi1 --loadfmap=' (dir_phase) 'fmap_roi1_rads_erd --dwell=0.0001525 --asym=0.0043 --poly=3 -u  ' (dir1) 'AP_roi1_FM_erd_ply'])
unix(['fugue -i ' (dir2) 'AP_roi1 --loadfmap=' (dir_phase) 'fmap_roi1_rads_erd --dwell=0.0003050 --asym=0.0043 --poly=3 -u  ' (dir2) 'AP_roi1_FM_erd_ply'])

%%Fugue with intensity correction (icorr)
unix(['fugue -i ' (dir1) 'PA_roi1 --loadfmap=' (dir_phase) 'fmap_roi1_rads_erd --dwell=0.0001525 --asym=0.0043 --despike --despikethreshold=2.1 --unwarpdir=y- --icorr -u  ' (dir1) 'PA_roi1_FM_erd_spike_corr']) % intensity correction
unix(['fugue -i ' (dir2) 'PA_roi1 --loadfmap=' (dir_phase) 'fmap_roi1_rads_erd --dwell=0.0003050 --asym=0.0043 --despike --despikethreshold=2.1 --unwarpdir=y- --icorr -u  ' (dir2) 'PA_roi1_FM_erd_spike_corr'])
unix(['fugue -i ' (dir1) 'AP_roi1 --loadfmap=' (dir_phase) 'fmap_roi1_rads_erd --dwell=0.0001525 --asym=0.0043 --despike --despikethreshold=2.1 --unwarpdir=y  --icorr -u  ' (dir1) 'AP_roi1_FM_erd_spike_corr']) % intensity correction
unix(['fugue -i ' (dir2) 'AP_roi1 --loadfmap=' (dir_phase) 'fmap_roi1_rads_erd --dwell=0.0003050 --asym=0.0043 --despike --despikethreshold=2.1 --unwarpdir=y  --icorr -u  ' (dir2) 'AP_roi1_FM_erd_spike_corr'])


%%PA-AP subtraction
unix(['fslmaths ' (dir1) 'AP_roi1 -sub ' (dir1) 'PA_roi1 '  (dir1) 'Sub_AP_PA_roi1']) %% no correction

unix(['fslmaths ' (dir1) 'AP_roi1_warp01 -sub ' (dir1) 'PA_roi1_warp02 '  (dir1) 'Sub_AP_PA_roi1_warp']) %% unwarping using topup field map
unix(['fslmaths ' (dir1) 'AP_roi1_warp01_jac -sub ' (dir1) 'PA_roi1_warp02_jac '  (dir1) 'Sub_AP_PA_roi1_warp_jac'])%% unwarping & intensity correction using topup 

unix(['fslmaths ' (dir1) 'AP_roi1_FM_erd_spike -sub ' (dir1) 'PA_roi1_FM_erd_spike '  (dir1) 'Sub_AP_PA_roi1_FM_erd_spike'])%% unwarping using B0 field map
unix(['fslmaths ' (dir1) 'AP_roi1_FM_erd_spike_corr -sub ' (dir1) 'PA_roi1_FM_erd_spike_corr '  (dir1) 'Sub_AP_PA_roi1_FM_erd_spike_corr'])%% unwarping & intensity correction using B0 field map

%% Low-resolution B0 field map (2*2*2 mm3)

%%registeration to higher resolution data
unix(['flirt -ref ' (dir_magn) 'magn1_roi1 -in ' (dir_magn_low) 'magn1_roi1 -omat ' (dir_magn_low) 'my_affine_transf.mat -dof 6 -nosearch'])
unix(['flirt -in ' (dir_magn_low) 'magn1_roi1 -ref ' (dir_magn) 'magn1_roi1 -out ' (dir_magn_low) 'magn1_reg_roi1 -init ' (dir_magn_low) 'my_affine_transf.mat -applyxfm'])
unix(['flirt -in ' (dir_phase_low) 'phase_roi1 -ref ' (dir_phase) 'phase_roi1 -out ' (dir_phase_low) 'phase_reg_roi1 -init ' (dir_magn_low) 'my_affine_transf.mat -applyxfm'])

%%Phantom magnitude mask
unix(['fslmaths ' (dir_magn_low) 'magn1_reg_roi1 -thr 200 ' (dir_magn_low) 'magn1_reg_roi1_thr'])
unix(['fslmaths ' (dir_magn_low) 'magn1_reg_roi1_thr -bin ' (dir_magn_low) 'magn1_reg_roi1_thr_bin'])
unix(['fslmaths ' (dir_magn_low) 'magn1_reg_roi1 -mul ' (dir_magn_low) 'magn1_reg_roi1_thr_bin '  (dir_magn_low) 'magn1_reg_roi1_brain'])
unix(['fslmaths ' (dir_magn_low) 'magn1_reg_roi1_brain.nii.gz -ero ' (dir_magn_low) 'magn1_reg_roi1_brain_erd.nii.gz'])

%%unwrapping phase data
unix(['fsl_prepare_fieldmap SIEMENS ' (dir_phase_low) 'phase_reg_roi1.nii.gz ' (dir_magn_low) 'magn1_reg_roi1_brain.nii.gz ' (dir_phase_low) 'fmap_reg_roi1_rads.nii.gz 3.1'])
unix(['fslmaths ' (dir_phase_low) 'fmap_reg_roi1_rads.nii.gz -div 6.28 ' (dir_phase_low) 'fmap_reg_roi1_hz.nii.gz'])

unix(['fsl_prepare_fieldmap SIEMENS ' (dir_phase_low) 'phase_reg_roi1.nii.gz ' (dir_magn_low) 'magn1_reg_roi1_brain_erd.nii.gz ' (dir_phase_low) 'fmap_reg_roi1_rads_erd.nii.gz 3.1'])
unix(['fslmaths ' (dir_phase_low) 'fmap_reg_roi1_rads_erd.nii.gz -div 6.28 ' (dir_phase_low) 'fmap_reg_roi1_hz_erd.nii.gz'])

%%unwarp EPI data
unix(['fugue -i ' (dir1) 'AP_roi1 --loadfmap=' (dir_phase_low) 'fmap_reg_roi1_rads --dwell=0.0001525 --asym=0.0031 --despike --despikethreshold=2.1  -u ' (dir1) 'AP_roi1_reg_FM_spike'])
unix(['fugue -i ' (dir2) 'AP_roi1 --loadfmap=' (dir_phase_low) 'fmap_reg_roi1_rads --dwell=0.0003050 --asym=0.0031 --despike --despikethreshold=2.1   -u ' (dir2) 'AP_roi1_reg_FM_spike'])

unix(['fugue -i ' (dir1) 'AP_roi1 --loadfmap=' (dir_phase_low) 'fmap_reg_roi1_rads_erd --dwell=0.0001525 --asym=0.0031 --despike --despikethreshold=2.1  -u  ' (dir1) 'AP_roi1_reg_FM_erd_spike'])
unix(['fugue -i ' (dir2) 'AP_roi1 --loadfmap=' (dir_phase_low) 'fmap_reg_roi1_rads_erd --dwell=0.0003050 --asym=0.0031 --despike --despikethreshold=2.1  -u  ' (dir2) 'AP_roi1_reg_FM_erd_spike'])
unix(['fslmaths ' (dir_phase_low) 'fmap_reg_roi1_hz -mul ' (dir_magn) 'magn1_roi1_thr_bin '  (dir_phase_low) 'fmap_reg_roi1_hz_masked'])


%% Computing Dice in two phantom regions

dir1=  '/export/data/vmalekian/test/phantom/phantom_20220412/nc_epi3d_v2u_0pt92_WB_Seg8_AP_PA_TE25_160s_0012/';
dir2=  '/export/data/vmalekian/test/phantom/phantom_20220412/nc_epi3d_v2u_0pt92_WB_Seg4_AP_PA_TE25_160s_0013/';
dir_magn=  '/export/data/vmalekian/test/phantom/phantom_20220412/aa_b0_Sag_WE_08mm_0014/';
dir_phase= '/export/data/vmalekian/test/phantom/phantom_20220412/aa_b0_Sag_WE_08mmB0_Map_0015/';

%%B0 magnitude data
string1 = 'magn1_roi1_thr_bin_dilD.nii.gz';
a0m=dir([dir_magn,string1]);
nf=load_untouch_nii([dir_magn,a0m(1).name]);
mask = imrotate((nf.img)>0,90);

string1 = 'magn1_roi1_brain.nii.gz';
a0m=dir([dir_magn,string1]);
nf=load_untouch_nii([dir_magn,a0m(1).name]);
str = imrotate(single(nf.img),90);
structure = imrotate(single(nf.img),90);

%%create mask from B0 magnitude data
thresh=200;P=1000;conp=6;
str_bin = (bwareaopen(str,P,conp).*mask)>0;

%%Two-fold segmented EPI
string1 = 'AP_roi1.nii.gz';
a0m=dir([dir1,string1]);
nf=load_untouch_nii([dir1,a0m(1).name]);
seg2 = imrotate(single(nf.img),90);
seg2_bin = imbinarize(seg2,thresh).*mask;
seg2_bin = bwareaopen(seg2_bin,P,conp);

string1 = 'AP_dc_def_roi1_jac.nii.gz';
a0m=dir([dir1,string1]);
nf=load_untouch_nii([dir1,a0m(1).name]);
seg2_bp_jac = imrotate(single(nf.img),90);
seg2_bp_jac_bin = imbinarize(seg2_bp_jac,thresh).*mask;
seg2_bp_jac_bin = bwareaopen(seg2_bp_jac_bin,P,conp);

string1 = 'AP_roi1_FM_erd_spike.nii.gz';
a0m=dir([dir1,string1]);
nf=load_untouch_nii([dir1,a0m(1).name]);
seg2_b0 = imrotate(single(nf.img),90);
seg2_b0_bin = imbinarize(seg2_b0,thresh).*mask;
seg2_b0_bin = bwareaopen(seg2_b0_bin,P,conp);

%%Non-segmented EPI
string1 = 'AP_roi1.nii.gz';
a0m=dir([dir2,string1]);
nf=load_untouch_nii([dir2,a0m(1).name]);
seg1 = imrotate(single(nf.img),90);
seg1_bin = imbinarize(seg1,thresh).*mask;
seg1_bin = bwareaopen(seg1_bin,P,conp);

string1 = 'AP_dc_def_roi1_jac.nii.gz';
a0m=dir([dir2,string1]);
nf=load_untouch_nii([dir2,a0m(1).name]);
seg1_bp_jac = imrotate(single(nf.img),90);
seg1_bp_jac_bin = imbinarize(seg1_bp_jac,thresh).*mask;
seg1_bp_jac_bin = bwareaopen(seg1_bp_jac_bin,P,conp);

string1 = 'AP_roi1_FM_erd_spike.nii.gz';
a0m=dir([dir2,string1]);
nf=load_untouch_nii([dir2,a0m(1).name]);
seg1_b0 = imrotate(single(nf.img),90);
seg1_b0_bin = imbinarize(seg1_b0,thresh).*mask;
seg1_b0_bin = bwareaopen(seg1_b0_bin,P,conp);

%% Phantom ROI DC calculation

%%ROI I
x1=70;
str_cut1= str_bin(1:x1,:,1+3:end-3);
figure,imshow(squeeze(str_cut1(:,:,11)),[])

seg1_bin_cut1 = seg1_bin(1:x1,:,1+3:end-3);
seg1_bp_jac_bin_cut1 = seg1_bp_jac_bin(1:x1,:,1+3:end-3);
seg1_b0_bin_cut1 = seg1_b0_bin(1:x1,:,1+3:end-3);
seg2_bin_cut1 = seg2_bin(1:x1,:,1+3:end-3);
seg2_bp_jac_bin_cut1 = seg2_bp_jac_bin(1:x1,:,1+3:end-3);
seg2_b0_bin_cut1 = seg2_b0_bin(1:x1,:,1+3:end-3);

dc_seg1_cut1= dice(seg1_bin_cut1,str_cut1);
dc_bp_jac_seg1_cut1= dice(seg1_bp_jac_bin_cut1,str_cut1);
dc_b0_seg1_cut1= dice(seg1_b0_bin_cut1,str_cut1);
dc_seg2_cut1= dice(seg2_bin_cut1,str_cut1);
dc_bp_jac_seg2_cut1= dice(seg2_bp_jac_bin_cut1,str_cut1);
dc_b0_seg2_cut1= dice(seg2_b0_bin_cut1,str_cut1);

DC_jac=[dc_bp_jac_seg1_cut1,dc_bp_jac_seg2_cut1];

%%ROI II
str_cut2= str_bin(71:end,:,1+3:end-3);
figure,imshow(squeeze(str_cut2(:,:,11)),[])

seg1_bin_cut2 = seg1_bin(71:end,:,1+3:end-3);
seg1_bp_jac_bin_cut2 = seg1_bp_jac_bin(71:end,:,1+3:end-3);
seg1_b0_bin_cut2 = seg1_b0_bin(71:end,:,1+3:end-3);
seg2_bin_cut2 = seg2_bin(71:end,:,1+3:end-3);
seg2_bp_jac_bin_cut2 = seg2_bp_jac_bin(71:end,:,1+3:end-3);
seg2_b0_bin_cut2 = seg2_b0_bin(71:end,:,1+3:end-3);

dc_seg1_cut2= dice(seg1_bin_cut2,str_cut2);
dc_bp_jac_seg1_cut2= dice(seg1_bp_jac_bin_cut2,str_cut2);
dc_b0_seg1_cut2= dice(seg1_b0_bin_cut2,str_cut2);
dc_seg2_cut2= dice(seg2_bin_cut2,str_cut2);
dc_bp_jac_seg2_cut2= dice(seg2_bp_jac_bin_cut2,str_cut2);
dc_b0_seg2_cut2= dice(seg2_b0_bin_cut2,str_cut2);

DC2_jac=[dc_bp_jac_seg1_cut2,dc_bp_jac_seg2_cut2];

