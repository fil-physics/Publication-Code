%% Normalized MI calculation for fMRI-EPI data run1

clear
close all
clc

addpath('/export/data/vmalekian/MI_Distortion/mi/')
addpath('/export/data/vmalekian/MI_Distortion/NIfTI_20140122/')
addpath('/export/data/vmalekian/MI_Distortion/InfoTheory/')


LL=1:22;
p=0;
thresh=0.8;

for sub= LL(~ismember(LL,[4,7])) %% subjects 4 and 7 were excluded due to extensive motion
        
        p=p+1;

        dir_seg_str= ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/MP2RAGE_WBIC_0pt65_PAT3_PF68_240Hz_UNI/'];

        
        string1 = 'c1clean_uni_reo.nii';
        a0m=dir([dir_seg_str,string1]);
        nf=load_untouch_nii([dir_seg_str,a0m(1).name]);
        gray = im2double(nf.img);
       
        string1 = 'clean_uni_reo_masck.nii.gz';
        a0m=dir([dir_seg_str,string1]);
        nf=load_untouch_nii([dir_seg_str,a0m(1).name]);
        mask = im2double(nf.img);
        
        dir_reg_nc= ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_noc_3/'];
        dir_reg_bp= ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_top_3/'];
        dir_reg_b0= ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_FM_3/'];
      
        
        string1 = 'c1initial_highres2highres.nii';
        a0m=dir([dir_reg_nc,string1]);
        nf=load_untouch_nii([dir_reg_nc,a0m(1).name]);
        c1_0 = im2double(nf.img);
        
        string1 = 'c1initial_highres2highres_jac.nii';
        a0m=dir([dir_reg_bp,string1]);
        nf=load_untouch_nii([dir_reg_bp,a0m(1).name]);
        c1_1 = im2double(nf.img);
        
        string1 = 'c1initial_highres2highres.nii';
        a0m=dir([dir_reg_b0,string1]);
        nf=load_untouch_nii([dir_reg_b0,a0m(1).name]);
        c1_2 = im2double(nf.img);
        
       
        dir_ats= ['/export/data/vmalekian/MT_Alice_dataset_20_revision/' sprintf('%02d',sub) '/fMRI0003/AP_noc_3/'];
        string_ats= 'oxf_thr0_2highres_nonlinear.nii.gz';
        a0m=dir([dir_ats,string_ats]);
        nf=load_untouch_nii([dir_ats,a0m(1).name]);
        ats = double(nf.img);
        regions= max(max(max(ats)));
        
        
        string1 = 'example_func2highres.nii.gz';
        a0m=dir([dir_reg_nc,string1]);
        nf=load_untouch_nii([dir_reg_nc,a0m(1).name]);
        vol_nc1 = double(nf.img);
        nc_th = mean(vol_nc1(mask==1))/10;
        nc_th_temp (sub)=nc_th;
        vol_nc=vol_nc1.*(vol_nc1>nc_th);
        
        string1 = 'example_func2highres.nii.gz';
        a0m=dir([dir_reg_bp,string1]);
        nf=load_untouch_nii([dir_reg_bp,a0m(1).name]);
        vol_bp1 = double(nf.img);    
        vol_bp=vol_bp1.*(vol_nc1>nc_th);


        string1 = 'example_func2highres.nii.gz';
        a0m=dir([dir_reg_b0,string1]);
        nf=load_untouch_nii([dir_reg_b0,a0m(1).name]);
        vol_b01 = double(nf.img);
        vol_b0=vol_b01.*(vol_nc1>nc_th);

        
        string1 = 'clean_uni_reo.nii';
        a0m=dir([dir_seg_str,string1]);
        nf=load_untouch_nii([dir_seg_str,a0m(1).name]);
        vol_str1 = double(nf.img);
        vol_str=vol_str1.*(vol_nc1>nc_th);

        
        
        for k1=1:regions            %% loop for each region
            BW = (ats==k1);
            
            
            gray_masked=gray;
            c0_masked = c1_0;
            c1_masked = c1_1;
            c2_masked = c1_2;
            
            gray_bin = imbinarize(gray_masked,thresh);
            c0_bin = imbinarize(c0_masked,thresh);
            c1_bin = imbinarize(c1_masked,thresh);
            c2_bin = imbinarize(c2_masked,thresh);
            
  
            BWmask = BW.*gray_bin;
            vol_nc_bw = ((vol_nc)>nc_th).*BWmask;
            mask_th = sum(sum(sum(vol_nc_bw)))/sum(sum(sum(BWmask)));
            mask_rev(sub,k1)=mask_th; 
            
            
            vol_str1 =gray_bin.*vol_str;
            vol_nc1 =c0_bin.*vol_nc;
            vol_bp1 =c1_bin.*vol_bp;
            vol_b01 =c2_bin.*vol_b0;
            
            
            vol_str2=vol_str1(BW==1);
            vol_nc2=vol_nc1(BW==1);
            vol_bp2=vol_bp1(BW==1);
            vol_b02=vol_b01(BW==1);
            

            R1(p,k1) = nmi(round(vol_nc2),round(vol_str2));
            R2(p,k1) = nmi(round(vol_bp2),round(vol_str2));
            R3(p,k1) = nmi(round(vol_b02),round(vol_str2));
            
            
        end
end
RR1 = cat(3,R1,R2,R3);
save('R_run1.mat','RR1')

%% NMI box plot figure

index = find(~(median(mask_rev)'> 0.75))
save('fMRI_index','index')
sum(~(median(mask_rev)'> 0.75))
index_in = median(mask_rev)'> 0.75; %% remove regions less than 0.75 overlap due to partial coverage of fMRI data

q=1;
RW=zeros(size(R1,1),3*regions);
for w1=1:regions
    q=3*(w1-1)+1;
    RW(:,q:q+2)=[R1(:,w1),R2(:,w1),R3(:,w1)];
end

RRW=RW;
roi_exc=index;
temp3=zeros(1,length(roi_exc)*3);
for kl=1:length(roi_exc)
    temp= roi_exc(kl);
    temp3(1+(kl-1)*3:kl*3)= 1+(temp-1)*3:temp*3;
end
RRW(:,temp3)=0;
figure,boxplot(RRW)

colors=[0,0,1;0,1,0;1,0,0];
h = findobj (gca, 'Tag' , 'Box' );
qq=1;
for w1=1:3:3*regions
    patch(get(h(w1),'XData'),get(h(w1),'YData'),colors(1,:),'FaceAlpha',.5);
    patch(get(h(w1+1),'XData'),get(h(w1+1),'YData'),colors(2,:),'FaceAlpha',.5);
    patch(get(h(w1+2),'XData'),get(h(w1+2),'YData'),colors(3,:),'FaceAlpha',.5);

end
% title('Normalized MI between fMRI EPI data and MP2RAGE in Harvard-Oxford cortical atlas regions'),
xlabel('ROI numbers based on atlas regions'),ylabel('NMI')
ylim([0 .8]),legend('B0 field mapping','Reversed-PE','No Correction')
set(gca, 'XTick', 1+1:3:3*regions+1);
set(gca, 'YTick', 0:0.1:1);
set(gca,'FontSize',16),set(gcf,'Color',[1 1 1])


med_R1 =median(R1);
med_R2 =median(R2);
med_R3 =median(R3);
iqr_R1 = iqr(R1);
iqr_R2 = iqr(R2);
iqr_R3 = iqr(R3);

aa1=round(med_R1,2)
aa2=round(med_R2,2)
aa3=round(med_R3,2)
bb1=round(iqr_R1,2)
bb2=round(iqr_R2,2)
bb3=round(iqr_R3,2)

%%ANOVA analysis
figure,
for r=1:regions
    a1 = [R1(:,r)',R2(:,r)',R3(:,r)']';
    a2 = [1*ones(1,size(R1,1)),2*ones(1,size(R1,1)),3*ones(1,size(R1,1))]';
    [pval,tbl,stats] = anova1(a1,a2,'off');
    Fstat(r) = tbl{2,5};
    Pstat(r) = tbl{2,6};
   [Xc,Xm,Xh,Xnms] = multcompare(stats);
   XcT(:,:,r)= Xc;
   XmT(:,:,r)= Xm;
end

XcP = squeeze(XcT(:,6,:));
XcP_selected = (XcP<0.05);
sum(XcP_selected')

selected_Pstat= (Pstat<0.05);
sum(selected_Pstat)



%% Improvement map

NT1 = ((med_R2-med_R1)./med_R1)*100;
NT2 = ((med_R3-med_R1)./med_R1)*100;
NT1(isnan(NT1))=0;
NT2(isnan(NT2))=0;
NT1(roi_exc)=0;
NT2(roi_exc)=0;


NT=[NT1;NT2].*[index_in';index_in']
result=cat(1,Fstat,Pstat,selected_Pstat,NT);
xlswrite('result_fMRI_run1.xlsx',result)

dir_ats = '/export/data/vmalekian/MT_Alice_dataset_20_revision/atlas/';
% string_ats = 'MNI-maxprob-thr25-1mm.nii.gz';
string_ats = 'HarvardOxford-cort-maxprob-thr0-1mm.nii.gz';
a0m=dir([dir_ats,string_ats]);
nf=load_untouch_nii([dir_ats,a0m(1).name]);
ats = double(nf.img);
regions= max(max(max(ats)));

dir_std = '/export/data/vmalekian/MT_Alice_dataset_20_revision/standard/';
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
%%final maps for both B0 and reversed-PE maps
lbtmean1=sum(lbt1,4);
string_magn = 'NTbp_th0_oxf_fmri_run1_20s_med.nii.gz';
save_nii(make_nii(single(lbtmean1)),string_magn)
nf1=load_untouch_nii('NTbp_th0_oxf_fmri_run1_20s_med.nii.gz');
nf1.hdr.hist = nf.hdr.hist;
save_untouch_nii(nf1,string_magn)

lbtmean2=sum(lbt2,4);
string_magn = 'NTb0_th0_oxf_fmri_run1_20s_med.nii.gz';
save_nii(make_nii(single(lbtmean2)),string_magn)
nf1=load_untouch_nii('NTb0_th0_oxf_fmri_run1_20s_med.nii.gz');
nf1.hdr.hist = nf.hdr.hist;
save_untouch_nii(nf1,string_magn)

