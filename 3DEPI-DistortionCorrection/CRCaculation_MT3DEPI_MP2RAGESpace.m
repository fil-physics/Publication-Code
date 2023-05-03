%% Correlation ratio (CR) calculation of MT-3DEPI data in MP2RAGE space 
% No-correction, B0 field mapping and reversed-PE (fsl_topup) techniques
% V Malekian, FIL Physics 

clc;clear;close all;
addpath('/export/home/vmalekian/MI_Distortion/NIfTI_20140122/')

SUBJT=1:22;
p=0;

for sub= SUBJT(~ismember(SUBJT,[4,7])) %% subjects 4 and 7 were excluded due to extensive motion [4,7]
    p=p+1;

    dir_reg_nc= ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_noc_3/'];
    dir_reg_bp= ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_top_3/'];
    dir_reg_b0= ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_FM_3/'];

    dir_ats = ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_noc_3/'];
    string_ats= 'oxf_thr0_2highres_nonlinear.nii.gz'; %% atlas in MP2RAGE space for parcellation
    a0m=dir([dir_ats,string_ats]);
    nf=load_untouch_nii([dir_ats,a0m(1).name]);
    ats = (nf.img);
    regions= max(max(max(ats)));

    string1 = 'initial_highres2highres.nii'; %% MT-3DEPI data in MP2RAGE space using no correction
    a0m=dir([dir_reg_nc,string1]);
    nf=load_untouch_nii([dir_reg_nc,a0m(1).name]);
    vol_nc = double(nf.img);


    string1 = 'initial_highres2highres_jac.nii'; %% MT-3DEPI data in MP2RAGE space using reversed-PE
    a0m=dir([dir_reg_bp,string1]);
    nf=load_untouch_nii([dir_reg_bp,a0m(1).name]);
    vol_bp = double(nf.img);

    string1 = 'initial_highres2highres.nii'; %% MT-3DEPI data in MP2RAGE space using B0 field mapping
    a0m=dir([dir_reg_b0,string1]);
    nf=load_untouch_nii([dir_reg_b0,a0m(1).name]);
    vol_b0 = double(nf.img);

   


    dir_seg_str= ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/MP2RAGE_WBIC_0pt65_PAT3_PF68_240Hz_UNI/'];
    string1 = 'c1clean_uni_reo_bin_edge_dil_thr_mask.nii.gz'; %% GM boundary mask 
    
    a0m=dir([dir_seg_str,string1]);
    nf=load_untouch_nii([dir_seg_str,a0m(1).name]);
    gray = im2double(nf.img);

    string1 = 'clean_uni_reo.nii';  %% MP2RAGE data as a reference
    a0m=dir([dir_seg_str,string1]);
    nf=load_untouch_nii([dir_seg_str,a0m(1).name]);
    vol_str = double(nf.img);


    for k1=1:regions


        BW = ((ats==k1).*gray);

        vol_str2=round(vol_str(BW==1));
        vol_nc2=round(vol_nc(BW==1));
        vol_bp2=round(vol_bp(BW==1));
        vol_b02=round(vol_b0(BW==1));

        R1(p,k1) = correlation_ratio((vol_nc2),(vol_str2));
        R2(p,k1) = correlation_ratio((vol_bp2),(vol_str2));
        R3(p,k1) = correlation_ratio((vol_b02),(vol_str2));

    end

end

RR123 = cat(3,R1,R2,R3);
save('CR_WB.mat','RR123')

%%Check if CR distrubation is noraml
for r=1:regions
    h_R1(r)= chi2gof(R1(:,r));
    h_R2(r)= chi2gof(R2(:,r));
    h_R3(r)= chi2gof(R3(:,r));
end
figure,plot(h_R1,'--g'),hold on,plot(h_R2,'b-'),plot(h_R3,'r-.'),legend('NC','BP','B0')

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
 errorbar(i, mean(RW(:,i)),-1*std(RW(:,i)),1*std(RW(:,i)),'Marker','*','Color', colorss(i),"LineWidth",3,"CapSize",5)
end
hold off


xlabel('ROI numbers based on atlas regions'),ylabel('CR')
legend('No Correction','Reversed-PE','B0 field mapping');
colororder({'k','k'})
yyaxis right,ylim([0.1 0.9])
set(gca, 'XTick', 1+1:3:3*regions+1);xticklabels({1:48})
set(gca, 'YTick', 0.1:0.1:1);yyaxis left,ylim([0.1 0.9])
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

XcP = squeeze(XcT(:,6,:));
XcP_selected = (XcP<0.05);

selected_Pstat= (Pstat<0.05);
XcP = squeeze(XcT(:,6,:));
FF =round(Fstat,1);
FF(selected_Pstat);

%% Relative CR map

NT1 = ((mean_R2-mean_R1)./mean_R1)*100;
NT2 = ((mean_R3-mean_R1)./mean_R1)*100;

NT=[NT1;NT2];
dir_ats = '/export/home/vmalekian/MT_Alice_dataset_20_revision/atlas/';
string_ats = 'HarvardOxford-cort-maxprob-thr0-1mm.nii.gz';
a0m=dir([dir_ats,string_ats]);
nf=load_untouch_nii([dir_ats,a0m(1).name]);
ats = double(nf.img);
regions= max(max(max(ats)));

dir_std = '/export/home/vmalekian/MT_Alice_dataset_20_revision/standard/'; %%to display on group-level GM transformed into standard space
string_st = 'c1clean_uni_reo2standard_no_total_mean.nii.gz';
a0m=dir([dir_std,string_st]);
nf=load_untouch_nii([dir_std,a0m(1).name]);
st = double(nf.img);
st=st>0.75;

[M1,N1,P1]=size(st);
lbt1=zeros(M1,N1,P1,regions);
lbt2=zeros(M1,N1,P1,regions);

for kk=1:regions
    label = (ats==kk);
    lbt1(:,:,:,kk)= NT1(kk).*label.*st;
    lbt2(:,:,:,kk)= NT2(kk).*label.*st;

end

lbtmean1=sum(lbt1,4);

%%Group-level GM edge extraction for relative CR map

lbtmean1_temp = lbtmean1;
BW_edge = abs(lbtmean1)>0;
CC_edge = bwconncomp(BW_edge);
numPixels = cellfun(@numel,CC_edge.PixelIdxList);
[biggest,idx] = max(numPixels);
lbtmean1_temp(CC_edge.PixelIdxList{idx}) = 1000;
BW_edge= lbtmean1_temp>500;
BW_edge_final = edge3(BW_edge,'sobel',0.5);
string_magn = 'tempt_bound.nii.gz';
save_nii(make_nii(single(BW_edge)),string_magn)
nf1=load_untouch_nii('tempt_bound.nii.gz');
nf1.hdr.hist = nf.hdr.hist;
save_untouch_nii(nf1,string_magn)

string_magn = 'temp_bound.nii.gz';
save_nii(make_nii(single(abs(lbtmean1)>0)),string_magn)
nf1=load_untouch_nii('temp_bound.nii.gz');
nf1.hdr.hist = nf.hdr.hist;
save_untouch_nii(nf1,string_magn)
dir_atst = '/export/home/vmalekian/MT_Alice_dataset_20_revision/code_backup/';
unix(['fslmaths ' (dir_atst) 'temp_bound -bin ' (dir_atst) 'temp_bound_bin'])
unix(['fslmaths ' (dir_atst) 'temp_bound_bin -edge ' (dir_atst) 'temp_bound_bin_edge'])
unix(['fslmaths ' (dir_atst) 'temp_bound_bin_edge -thr 0.4 -bin ' (dir_atst) 'temp_bound_bin_edge_bin'])


%%Final maps for both B0 field mapping and reversed-PE 
lbtmean1=sum(lbt1,4);
string_magn = 'WB_NTbp_th0_oxf.nii.gz';
save_nii(make_nii(single(lbtmean1)),string_magn)
nf1=load_untouch_nii('WB_NTbp_th0_oxf.nii.gz');
nf1.hdr.hist = nf.hdr.hist;
save_untouch_nii(nf1,string_magn)

lbtmean2=sum(lbt2,4);
string_magn = 'WB_NTb0_th0_oxf.nii.gz';
save_nii(make_nii(single(lbtmean2)),string_magn)
nf1=load_untouch_nii('WB_NTb0_th0_oxf.nii.gz');
nf1.hdr.hist = nf.hdr.hist;
save_untouch_nii(nf1,string_magn)


