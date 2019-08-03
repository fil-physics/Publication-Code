function [unph,reliability_msk] = m1_shf_3DSRNCP(ph_corr,mask,OS_type)
% DESCRIPTION
% 
% Calls the best-path phase unwrapping algorithm described in "Fast 
% three-dimensional phase-unwrapping algorithm based on sorting by 
% reliability following a noncontinuous path". Abdul-Rahman H, et al.
% 
% INPUT ARGUMENTS
% 
%  ph_corr:     wrapped phase
%  mask:        ROI mask
%  OS_type:     3DSRNCP precompiled version
%                [2.1]  Debian Linux,  Ubuntu 16.04 (default)
%                [2.11] Debian Linux,  Ubuntu 14.04
%                [2.12] Debian Linux,  Debian 7 Wheezy
%                [2.2]  Red Hat Linux, CentOS 6.9
%                [2.3]  OSX,           10.13.6
% 
% NOTES
% 
%  The 3DSRNCP unwrapper was originally written by the authors in C++
%  The code was modified by HF Sun to output 3D second differences
%  If you hit issues with the precompiled versions of 3DSRNCP, you can
%  generate a binary executable for your own OS using Intel's C++ compiler 
%  The C++ code can be found in the present directory
%  If you do so, do not forget to add a new entry in m1_shf_3DSRNCP.m 
%  (switch OS_type), and update ptb.opt.OS in ptbs_XXX.m’s CUSTOM SETTINGS
%
%  Install Intel C++ compiler from Parallel Studio XE
%   https://software.intel.com/en-us/qualify-for-free-software
%
%  Compile, e.g.
%   icpc -std=c++11 3DSRNCP_jac.cpp -o 3DSRNCP_deb
% 
% Created by Hongfu Sun
% Last modified by Julio Acosta-Cabronero

if nargin<3
    OS_type = 2.1; % Debian Linux
end

imsize = size(ph_corr);

mask_unwrp = uint8(abs(mask)*255);
fid = fopen('mask_unwrp.dat','w');
fwrite(fid,mask_unwrp,'uchar');
fclose(fid);

pathstr = fileparts(which('3DSRNCP.m'));
setenv('pathstr',pathstr);
setenv('nv',num2str(imsize(1)));
setenv('np',num2str(imsize(2)));
setenv('ns',num2str(imsize(3)));

unph            = zeros(imsize);
reliability_msk = ones(size(ph_corr(:,:,:,1)));

for echo_num = 1:size(ph_corr,4)
    if size(ph_corr,4)>1
        disp(['Echo #' num2str(echo_num)])
    end
    setenv('echo_num',num2str(echo_num));
    fid = fopen(['wrapped_phase' num2str(echo_num) '.dat'],'w');
    fwrite(fid,ph_corr(:,:,:,echo_num),'float');
    fclose(fid);

    switch OS_type
        case 2.1
            unwrapper = '3DSRNCP_deb';
        case 2.11
            unwrapper = '3DSRNCP_ubuntu-14.04';
	    case 2.12
	        unwrapper = '3DSRNCP_deb7_wheezy';
        case 2.2
            unwrapper = '3DSRNCP_rh';
        case 2.3
            unwrapper = '3DSRNCP_osx-10.13.6';
    end
    
    bash_script = ['chmod +x ${pathstr}/' unwrapper...
       '; ${pathstr}/' unwrapper ' wrapped_phase${echo_num}.dat mask_unwrp.dat'...
       ' unwrapped_phase${echo_num}.dat $nv $np $ns reliability${echo_num}.dat'];
    unix(bash_script);

    fid = fopen(['unwrapped_phase' num2str(echo_num) '.dat'],'r');
    tmp = fread(fid,'float');
    unph(:,:,:,echo_num) = reshape(tmp - round(mean(tmp(mask==1))/(2*pi))*2*pi ,imsize(1:3)).*mask;
    fclose(fid);

    fid = fopen(['reliability' num2str(echo_num) '.dat'],'r');
    reliability_raw = fread(fid,'float');
    reliability_raw = reshape(reliability_raw,imsize(1:3));
    fclose(fid);
   
    reliability_msk(:,:,:,echo_num) = mask.*reliability_raw;     
end

% Housekeeping
delete *.dat
