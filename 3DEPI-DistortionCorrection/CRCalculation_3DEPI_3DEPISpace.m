%% Correlation ratio (CR) calculation of distorted MP2RAGE data in 3DEPI space 
% No-correction, B0 field mapping and reversed-PE (fsl_topup) techniques
% V Malekian, FIL Physics 

clc;clear;close all;
addpath('/export/home/vmalekian/MI_Distortion/NIfTI_20140122/')

SUBJT=1:22;
p=0;


for sub= SUBJT(~ismember(SUBJT,[4,7])) %% subjects 4 and 7 were excluded due to extensive motion [4,7]
    p=p+1;

    dir_seg_str= ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/'];
    dir_reg_nc= ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_noc_3/'];
    dir_reg_bp= ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_top_3/'];
    dir_reg_b0= ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_FM_3/'];

    string1 = 'AP_brain_mask.nii.gz';
    a0m=dir([dir_seg_str,string1]);
    nf=load_untouch_nii([dir_seg_str,a0m(1).name]);
    mask = double(nf.img);

    string1 = 'c1AP_bin_edge_dil_thr_mask_2fMREPI.nii.gz'; %% GM boundary mask
    a0m=dir([dir_reg_nc,string1]);
    nf=load_untouch_nii([dir_reg_nc,a0m(1).name]);
    gray = im2double(nf.img).*mask;

    
    dir_ats= ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_noc_3/'];
    string_ats= 'oxf_thr0_2func.nii.gz'; %% atlas in fMRI space for parcellation
    a0m=dir([dir_ats,string_ats]);
    nf=load_untouch_nii([dir_ats,a0m(1).name]);
    ats = double(nf.img);
    regions= max(max(max(ats)));


    string1 = 'AP.nii';  %% fMRI data as a reference
    a0m=dir([dir_seg_str,string1]);
    nf=load_untouch_nii([dir_seg_str,a0m(1).name]);
    vol_str = double(nf.img);

    string1 = 'clean_uni_reo_brain2fMREPI.nii.gz';%% distorted MP2RAGE data in fMRI space using no correction
    a0m=dir([dir_reg_nc,string1]);
    nf=load_untouch_nii([dir_reg_nc,a0m(1).name]);
    vol_nc = double(nf.img);
    nc_th = mean(vol_nc(mask==1))/10;

    string1 = 'clean_uni_reo_brain2fMREPI_rev.nii.gz';%% distorted MP2RAGE data in fMRI space using reversed-PE
    %'clean_uni_reo_brain2EPI.nii.gz';
    a0m=dir([dir_reg_bp,string1]);
    nf=load_untouch_nii([dir_reg_bp,a0m(1).name]);
    vol_bp = double(nf.img);
  
    string1 = 'clean_uni_reo_brain2fMREPI.nii.gz';%% distorted MP2RAGE data in fMRI space using B0 field mapping
    a0m=dir([dir_reg_b0,string1]);
    nf=load_untouch_nii([dir_reg_b0,a0m(1).name]);
    vol_b0 = double(nf.img);


    for k1=1:regions
        
        BW = ((ats==k1).*gray);

        BWmask = BW.*gray;
        vol_nc_bw = ((vol_nc)>nc_th).*BWmask;
        mask_th = sum(sum(sum(vol_nc_bw)))/sum(sum(sum(BWmask)));
        mask_rev(sub,k1)=mask_th;
           
        
        vol_str2=round((vol_str(BW==1)));
        vol_nc2=round(vol_nc(BW==1));
        vol_bp2=round(vol_bp(BW==1));
        vol_b02=round(vol_b0(BW==1));

        %%CR calculation
        R1(p,k1) = correlation_ratio((vol_str2),(vol_nc2));
        R2(p,k1) = correlation_ratio((vol_str2),(vol_bp2));
        R3(p,k1) = correlation_ratio((vol_str2),(vol_b02));

    end

end

load('/export/home/vmalekian/MT_Alice_dataset_20_revision/code_backup/pulication_code/fMRI_index_in.mat')

R1 =R1.*index_in';
R2 =R2.*index_in';
R3 =R3.*index_in';
RR123 = cat(3,R1,R2,R3);
save('CR_3DEPI_run1_native.mat','RR123')

%% CR boxplot figure
RW=zeros(size(R1,1),3*regions);
for w1=1:regions
    q=3*(w1-1)+1;
    RW(:,q:q+2)=[R1(:,w1),R2(:,w1),R3(:,w1)];

end

color=['r','g','b'];
colorss = repmat(color,1,48);
figure
hold on,
for i = 1:3*regions
 errorbar(i, mean(RW(:,i)),-1*std(RW(:,i)),1*std(RW(:,i)),'*','Color', colorss(i),"LineWidth",3,"CapSize",5)

end
hold off

xlabel('ROI numbers based on atlas regions'),ylabel('CR')
legend('No Correction','Reversed-PE','B0 field mapping');
colororder({'k','k'})
yyaxis right,ylim([0.1 0.9])
set(gca, 'XTick', 1+1:3:3*regions+1);xticklabels({1:48})
colororder({'k','k'})
set(gca, 'YTick', 0:0.1:1);yyaxis left,ylim([0.1 0.9])
set(gca,'FontSize',16),set(gcf,'Color',[1 1 1])
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'on';

mean_R1 =mean(R1);
mean_R2 =mean(R2);
mean_R3 =mean(R3);
std_R1 = std(R1);
std_R2 = std(R2);
std_R3 = std(R3);
aa1=round(mean_R1,3);
aa2=round(mean_R2,3);
aa3=round(mean_R3,3);
bb1=round(std_R1,3);
bb2=round(std_R2,3);
bb3=round(std_R3,3);

%% ANOVA analysis

for r=1:regions
    if index_in(r)
    a1 = [R1(:,r)',R2(:,r)',R3(:,r)']';
    a2 = [1*ones(1,size(R1,1)),2*ones(1,size(R1,1)),3*ones(1,size(R1,1))]';
    [pval,tbl,stats] = anova1(a1,a2,'off');
    Fstat(r) = tbl{2,5};
    Pstat(r) = tbl{2,6};
    Pstat1(r) = pval;
    [Xc,Xm,Xh,Xnms] = multcompare(stats);
    XcT(:,:,r)= Xc;
    XmT(:,:,r)= Xm;
    end
end

XcP = squeeze(XcT(:,6,:));
XcP_selected = (XcP<0.05)

selected_Pstat= (Pstat<0.05)
sum(selected_Pstat)
XcP = squeeze(XcT(:,6,:));
FF =round(Fstat,1);
FF(selected_Pstat);


