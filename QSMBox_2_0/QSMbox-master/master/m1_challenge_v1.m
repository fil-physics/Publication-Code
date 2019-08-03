function [rmse_all,hfen_all,ssim_all] = m1_challenge_v1(ptbsfile,lambda)
% Created by Juio Costa-Cabronero, 19/01/2019

%%
d0              = pwd;
ptbdir          = fileparts(which(mfilename));
challengedir    = [ptbdir '/ptb/_Challenge_v1'];

addpath([ptbdir '/ptb/_Utils'])
addpath([ptbdir '/ptb/_spm12b_nifti'])

%%
if nargin<1
    % Find PTBS file
    ptbsfile = dir('ptbs_*.m');

    if length(ptbsfile)>1
        disp('More than one PTBS file found - I will use the first one.')
    end

    ptbsfile = ptbsfile(1).name;
    ptbsfile = ptbsfile(1:end-2);
end

disp(['PTBS file set for evaluation: ' ptbsfile])

ptbsdir = [d0 '/' ptbsfile];

mkdir(ptbsdir)

copyfile([d0 '/' ptbsfile '.m'],[ptbsdir '/' ptbsfile '.m'])

%%
cd(challengedir)

% Load data
qsmroot     = 'chi_cosmos'; % 'chi_33';
[qsm_0,stp] = m1_nifti_load(qsmroot); 
[u_phase]   = m1_nifti_load('phs_tissue');
[msk,stm]   = m1_nifti_load('roi');
[magn]      = m1_nifti_load('magn');
sdiff_phase = m1_nifti_load('secdiff_phase');

cd(ptbsdir)

magn = msk.*magn;

% Constants
TE_ref      = 20e-3;
B0_ref      = 3;
CF          = B0_ref*42.57747892*1e6;

TE = 25e-3;
B0 = 2.89;

gyro = 2*pi*42.57747892;
ppm2rad = double(B0*gyro*TE);
rad2ppm = double(1/ppm2rad);

res         = [1.0625,1.0625,1.0714286];
matrix      = size(magn);
B0_dir      = [0 0 1];

u_phase = ppm2rad*msk.*u_phase;
u_phase = u_phase*(TE_ref/TE)*(B0_ref/B0);

m1_nifti_save(u_phase,stp,'u_phase')
m1_nifti_save(magn,stp,'magn')
m1_nifti_save(msk,stm,'roi')
m1_nifti_save(sdiff_phase,stp,'secdiff_phase')
m1_nifti_save(qsm_0,stp,qsmroot)

% Save params
PI = pi;
save('ptb_res.txt','res','-ascii');
save('ptb_B0.txt','B0_ref','-ascii');
save('ptb_B0_dir.txt','B0_dir','-ascii');
save('ptb_TE.txt','TE_ref','-ascii');
save('ptb_scaling_factor.txt','PI','-ascii');

%%
% Save/load ptb.mat
eval([ptbsfile '(''saveptb'')'])
load ptb.mat

if nargin>=2
    ptb.qsm.MSDI.lambda = lambda;
else
    lambda = ptb.qsm.MSDI.lambda;
end
    
ptb.pipeline = {'prep.load_custom1',...
                'prep.load_custom2',...
                'prep.load_custom3',...
                'prep.import_params',...
                'qsm.MSDI'};

% Image #1
ptb.prep.load_custom1.froot       = 'u_phase';
ptb.prep.load_custom1.varname     = 'u_phase';
ptb.prep.load_custom1.structname  = 'stp';

% Image #2
ptb.prep.load_custom2.froot       = 'magn';
ptb.prep.load_custom2.varname     = 'magn';
ptb.prep.load_custom2.structname  = 'stp';

% Image #3
ptb.prep.load_custom3.froot       = 'roi';
ptb.prep.load_custom3.varname     = 'msk';
ptb.prep.load_custom3.structname  = 'stm';

save('ptb.mat','ptb')

%%
eval([ptbsfile '(''runextptb'')'])

%%
reconfiles = dir('qsm_INTEGRAL_*.nii*');

for x = 1:length(reconfiles)
    reconroot = reconfiles(x).name;
    reconroot = m1_fname_split(reconroot,'gz');
    
    saveSSIMmap = true;
    
    [rmse,hfen,ssim] = m1_simil(reconroot,qsmroot,'roi',saveSSIMmap);
    
    simildir = [d0 '/' ptbsfile '/' reconroot];
    mkdir(simildir)
    movefile('simil_*.txt',simildir)
    
    if saveSSIMmap
        movefile('simil_ssim_map.nii*',simildir)
    end
    
    rmse_all(x) = rmse;
    hfen_all(x) = hfen;
    ssim_all(x) = ssim;
end

save('simil_rmse_all.txt','rmse_all','-ascii')
save('simil_hfen_all.txt','hfen_all','-ascii')
save('simil_ssim_all.txt','ssim_all','-ascii')

rmse_all
hfen_all
ssim_all

%%
cd(d0)
