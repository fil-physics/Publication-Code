function [ph_cmb,unph_diff_cmb,ph_offsets] = m1_shf_mgre_comb(img,vox,te,mask,OS_type)
% Multi-GRE phase combination
%
%   IMG        complex data, 4D/5D: [3D_image, echoes, receivers (optional)]
%   TE         echo times
%   vox        voxel size, e.g. [1 1 1]
%   PH_CMB     combined phase
%   OS_type    3DSRNCP unwrapper precompiled version
%               [2.10] Debian Linux,   Ubuntu 16.04 (default)
%               [2.11] Debian Linux,   Ubuntu 14.04
%               [2.12] Debian Linux,   Debian 7 Wheezy
%               [2.20] Red Hat Linux,  CentOS 6.9
%               [2.30] OSX,            10.13.6
%
% Adapted from code by HF Sun
% Last modified by Julio Acosta-Cabronero

if nargin<5
    OS_type = 2.10; % Debian Linux (Ubuntu 16.04)
end

[~,~,~,ne,nrcvrs] = size(img);
TE1 = te(1);
TE2 = te(2);
imsize = size(img);

img_diff = img(:,:,:,2,:)./img(:,:,:,1,:);
ph_diff = img_diff./abs(img_diff);
ph_diff_cmb = sum(abs(img(:,:,:,1,:)).*ph_diff,5);
ph_diff_cmb(isnan(ph_diff_cmb)) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% best-path unwrapping, locate 3DSRNCP
[pathstr] = fileparts(which('3DSRNCP.m'));
setenv('pathstr',pathstr);
setenv('nv',num2str(imsize(1)));
setenv('np',num2str(imsize(2)));
setenv('ns',num2str(imsize(3)));

fid = fopen('wrapped_phase_diff.dat','w');
fwrite(fid,angle(ph_diff_cmb),'float');
fclose(fid);

mask_unwrp = uint8(abs(mask)*255);
fid = fopen('mask_unwrp.dat','w');
fwrite(fid,mask_unwrp,'uchar');
fclose(fid);

switch OS_type
case 2.10
    unwrapper = '3DSRNCP_deb';
case 2.11
    unwrapper = '3DSRNCP_ubuntu-14.04';
case 2.12
    unwrapper = '3DSRNCP_deb7_wheezy';
case 2.20
    unwrapper = '3DSRNCP_rh';
case 2.30
    unwrapper = '3DSRNCP_osx-10.13.6';
end
bash_script = ['chmod +x ${pathstr}/' unwrapper...
               '; ${pathstr}/' unwrapper ' wrapped_phase_diff.dat mask_unwrp.dat '...
               'unwrapped_phase_diff.dat $nv $np $ns reliability_diff.dat'];
unix(bash_script);

fid = fopen(['unwrapped_phase_diff.dat'],'r');
tmp = fread(fid,'float');
unph_diff_cmb = reshape(tmp - round(mean(tmp(mask==1))/(2*pi))*2*pi ,imsize(1:3)).*mask;
fclose(fid);

% Housekeeping
delete *.dat
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

unph_te1_cmb = unph_diff_cmb*TE1/(TE2-TE1);
offsets = img(:,:,:,1,:)./repmat(exp(1j*unph_te1_cmb),[1,1,1,1,nrcvrs]);
offsets = offsets./abs(offsets);
offsets(isnan(offsets)) = 0;

for chan = 1:nrcvrs
    offsets(:,:,:,:,chan) = smooth3(offsets(:,:,:,:,chan),'box',round(10./vox/2)*2+1); 
    offsets(:,:,:,:,chan) = offsets(:,:,:,:,chan)./abs(offsets(:,:,:,:,chan));
end

ph_offsets = angle(offsets);

% combine phase according to complex summation
offsets = repmat(offsets,[1,1,1,ne,1]);
img = img./offsets;
ph_cmb = angle(sum(img,5));

ph_cmb(isnan(ph_cmb)) = 0;

