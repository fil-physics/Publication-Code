function m1_spmroi(thr,dilR,res)
% Generate brain mask for QSMbox using SPM tissue probability segments
% Run in SPM Segment directory
%
% INPUT
% thr	c1+c2+c3 probability cut-off threshold
% dilR  dilation (spherical) kernel radius in mm 
%        [dilR=0: no dilation; dilR=-1: dilation by one voxel]
% res   voxel resolution (required if non-zero dilR)
%
% Created by Julio Acosta-Cabronero, 08/12/2018

%% Arg
if nargin<1
    thr = 1e-4;
end

if nargin<2
    dilR = 0;
end

if dilR~=0 && nargin<3
    error('ERROR: m1_spmroi requires voxel resolution as input argument.')
end

%% Load data
disp('Load SPM segments')
% GM
fname = dir('c1*.nii*');
froot = m1_fname_split(fname.name,'gz');
[c1,st] = m1_nifti_load(froot);

% WM
fname = dir('c2*.nii*');
froot = m1_fname_split(fname.name,'gz');
[c2] = m1_nifti_load(froot);

% CSF
fname = dir('c3*.nii*');
froot = m1_fname_split(fname.name,'gz');
[c3] = m1_nifti_load(froot);

%% Processing
%
disp('Generate initial mask')
roi = c1+c2+c3;
roi(roi>thr)=1;
roi(roi<=thr)=0;

%
disp('Fill holes (3D approximation)')
roifillZ=roi*0;
for z=1:size(roi,3)
    roifillZ(:,:,z) = imfill(roi(:,:,z),'holes');
end
 
roifillY=roi*0;
for y=1:size(roi,2)
    roifillY(:,y,:) = imfill(squeeze(roi(:,y,:)),'holes');
end
roifill = roifillZ.*roifillY;
clear roifillZ roifillY

roifillX=roi*0;
for x=1:size(roi,1)
    roifillX(x,:,:) = imfill(squeeze(roi(x,:,:)),'holes');
end
roifill = roifill.*roifillX;
clear roifillX

% Dilate
dil_flag = true;

if dilR==-1
    R = 1;
elseif dilR==0
    dil_flag = false; 
else
    R = round(dilR/mean(res(:)));
end

if dil_flag
    disp(['Dilate ROI mask by ' num2str(R) ' voxel(s)'])
    seo = strel('sphere',R); % struct elem obj
    roifill = imdilate(roifill,seo);
end

%% Save brain mask
m1_nifti_save(roifill,st,'roi')
