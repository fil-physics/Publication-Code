%% Correlation ratio (CR) vs Mutual information: distorted MP2RAGE data in MT-3DEPI space 
% No-correction, B0 field mapping and reversed-PE (fsl_topup) techniques
% V Malekian, FIL Physics 

clc;clear;close all;
addpath('/export/home/vmalekian/MI_Distortion/NIfTI_20140122/')
addpath('/export/home/vmalekian/MI_Distortion/InfoTheory/')

SUBJT=1:22;
p=0;

for sub= SUBJT(~ismember(SUBJT,[4,7])) %% subjects 4 and 7 were excluded due to extensive motion [4,7]
    p=p+1;

    dir_seg_str= ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/nc_epi3d_v2d_MTw_WholeBrain_0005/'];
    string1 = 'c1AP_bin_edge_dil_thr_mask.nii.gz'; %% GM boundary mask
    a0m=dir([dir_seg_str,string1]);
    nf=load_untouch_nii([dir_seg_str,a0m(1).name]);
    gray = im2double(nf.img);


    dir_reg_nc= ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_noc_3/'];
    dir_reg_bp= ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_top_3/'];
    dir_reg_b0= ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_FM_3/'];


    dir_ats = ['/export/home/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_noc_3/']; %%2 ROIs: vmPFC and dmPFC
    string_ats1= 'vmPFC_2initialhighres.nii.gz';
    string_ats2= 'dmPFC_2initialhighres.nii.gz';


    a0m=dir([dir_ats,string_ats1]);
    nf=load_untouch_nii([dir_ats,a0m(1).name]);
    ats1 = (nf.img);

    a0m=dir([dir_ats,string_ats2]);
    nf=load_untouch_nii([dir_ats,a0m(1).name]);
    ats2 = (nf.img);
    regions=2;
    ats = cat(4,ats1,ats2);

    string1 = 'clean_uni_reo_brain2EPI.nii'; %% distorted MP2RAGE data in MT-3DEPI space (no correction)
    a0m=dir([dir_reg_nc,string1]);
    nf=load_untouch_nii([dir_reg_nc,a0m(1).name]);
    vol_nc = double(nf.img);

    string1 = 'clean_uni_reo_brain2EPI_rev.nii';%%  distorted MP2RAGE data in MT-3DEPI space (reversed-PE)
    a0m=dir([dir_reg_bp,string1]);
    nf=load_untouch_nii([dir_reg_bp,a0m(1).name]);
    vol_bp = double(nf.img);

    string1 = 'clean_uni_reo_brain2EPI.nii';%%  distorted MP2RAGE data in MT-3DEPI space (B0 field mapping)
    a0m=dir([dir_reg_b0,string1]);
    nf=load_untouch_nii([dir_reg_b0,a0m(1).name]);
    vol_b0 = double(nf.img);

    string1 = 'AP.nii';
    a0m=dir([dir_seg_str,string1]);
    nf=load_untouch_nii([dir_seg_str,a0m(1).name]);
    vol_str = double(nf.img);

    for k1=1:regions% vmpfc & dmpfc

        BW = (squeeze(ats(:,:,:,k1)).*gray);

        vol_str2=round(vol_str(BW==1));
        vol_nc2=round(vol_nc(BW==1));
        vol_bp2=round(vol_bp(BW==1));
        vol_b02=round(vol_b0(BW==1));


        %%MI calculation
        Rm1(p,k1) = mutualinfo((vol_nc2),(vol_str2));
        Rm2(p,k1) = mutualinfo((vol_bp2),(vol_str2));
        Rm3(p,k1) = mutualinfo((vol_b02),(vol_str2));

        %%CR calculation

        R1(p,k1) = correlation_ratio((vol_str2),(vol_nc2));
        R2(p,k1) = correlation_ratio((vol_str2),(vol_bp2));
        R3(p,k1) = correlation_ratio((vol_str2),(vol_b02));

    end

end

%% MI boxplot figure
RWm=zeros(size(Rm1,1),3*regions);
for w1=1:regions
    q=3*(w1-1)+1;
    RWm(:,q:q+2)=[Rm1(:,w1),Rm2(:,w1),Rm3(:,w1)];

end
figure,boxplot(RWm),grid('on')

colors=[0,0,1;0,1,0;1,0,0];
h = findobj (gca, 'Tag' , 'Box' );
qq=1;
for w1=1:3:3*regions
    patch(get(h(w1),'XData'),get(h(w1),'YData'),colors(1,:),'FaceAlpha',.5);
    patch(get(h(w1+1),'XData'),get(h(w1+1),'YData'),colors(2,:),'FaceAlpha',.5);
    patch(get(h(w1+2),'XData'),get(h(w1+2),'YData'),colors(3,:),'FaceAlpha',.5);

end

title('MI in MT-3DEPI space')
ylabel('MI')
legend('B0 field mapping','Reversed-PE','No Correction');

set(gca, 'XTick', 1+1:3:3*regions+1);
ylim([5 8])

set(gca,'FontSize',11),set(gcf,'Color',[1 1 1])
ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'off';
%%ANOVA analysis
figure,
for r=1:regions
    a1 = [Rm1(:,r)',Rm2(:,r)',Rm3(:,r)']';
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
XcP_selected = (XcP<0.05)
selected_Pstat= (Pstat<0.05);
XcP = squeeze(XcT(:,6,:));

%% CR boxplot figure
RW=zeros(size(R1,1),3*regions);
for w1=1:regions
    q=3*(w1-1)+1;
    RW(:,q:q+2)=[R1(:,w1),R2(:,w1),R3(:,w1)];

end
figure,boxplot(RW),grid('on')

colors=[0,0,1;0,1,0;1,0,0];
h = findobj (gca, 'Tag' , 'Box' );
qq=1;
for w1=1:3:3*regions
    patch(get(h(w1),'XData'),get(h(w1),'YData'),colors(1,:),'FaceAlpha',.5);
    patch(get(h(w1+1),'XData'),get(h(w1+1),'YData'),colors(2,:),'FaceAlpha',.5);
    patch(get(h(w1+2),'XData'),get(h(w1+2),'YData'),colors(3,:),'FaceAlpha',.5);

end

title('CR in MT-3DEPI space')
ylabel('CR')
legend('B0 field mapping','Reversed-PE','No Correction');

set(gca, 'XTick', 1+1:3:3*regions+1);
set(gca,'FontSize',11),set(gcf,'Color',[1 1 1])
ylim([0.15 0.85])

ax = gca;
ax.XGrid = 'on';
ax.YGrid = 'off';

%%ANOVA analysis
figure,
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
XcP_selected = (XcP<0.05)
selected_Pstat= (Pstat<0.05);
XcP = squeeze(XcT(:,6,:));
