% PHASE TOOLBOX
%
% Collection of processing pipelines for complex GRE data
%
%
% USAGE
%
%  1. Select a pipeline. Uncomment one of the suggested pipelines. If in
%     doubt, use the GENERIC routine suitable for your data
%  
%     If you feel more adventurous, you could create your own pipeline
%     under PROCESSING STEPS. See list of AVAILABLE METHODS inmediately
%     below. Create a 'ptb.pipename' entry in SELECT PIPELINE to enable
%     your new pipeline
%
%  2. Edit CUSTOM SETTINGS accordingly
%
%  3. Save & Run (F5)
%
%
% REQUIREMENTS
%
%  The 'ptb' directory must be added to the MATLAB path
%
%  If 'ptb.roi.bet2.run' is enabled in CUSTOM SETTINGS), FSL's brain 
%  extraction tool v2 (BET2) will be required (UNIX only)
% 
% 
% DEFAULTS
% 
%  Data input could be automated if the following two files existed in the
%  present working directory (pwd):
%
%  + magn_orig.nii(.gz)              magnitude data
% 
%  + phase_orig.nii(.gz)             phase data
%
%
%  The above datasets could be: 
% 
%  + Coil-combined single-TE GRE     3D
%
%  + Uncombined single-TE GRE        4D 
% 
%  + Coil-combined multi GRE         4D
%
% 
%  If 'ptb.roi.custom.run' is enabled in CUSTOM SETTINGS, the following 
%  file in pwd will skip manual input:
%
%  + roi.nii(.gz)                    ROI mask
% 
% 
%  Scan parameter input could also be automated if the following files were
%  available in pwd:
% 
%  + ptb_res.txt                     voxel size in mm, e.g. 1 1 1
% 
%  + ptb_TE.txt                      echo time(s) in s, e.g. 0.02 (20 ms)
% 
%  + ptb_B0.txt                      field strength in T, e.g. 3
%
%  + ptb_B0_dir.txt                  main field direction, e.g. 0 0 1
%
%  + ptb_scaling_factor.txt          phase scaling, e.g. 4096 (Siemens convention)
%
%
%  Alternatively, the above (except scaling factor) could be extracted from
%  a DICOM file. This should be included in a 'dicom' (small caps) dir
%
%
%  Otherwise, you will be prompted for manual input
% 
%
% Created by Julio Acosta-Cabronero, 2017.03.01

clear all

disp(' '); disp('PHASE TOOLBOX by Julio Acosta-Cabronero'); disp(' ')

%% SELECT PIPELINE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ptb.pipename = false;

% GENERIC routines %%%%%%%%%%%%%%%%%%%%%%%
% <<< generic single GRE only >>>
% ptb.pipename = {'gen_uncomb'}; 

% MSDI %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ptb.pipename = {'jac012_MSDI'};
% ptb.pipename = {'jac013_MSDI'};

% QSM INVERSION ONLY %%%%%%%%%%%%%%%%%%%%%
% ptb.pipename = {'qsm_MSDI170706'};
% ptb.pipename = {'qsm_fansi20170330'};
% ptb.pipename = {'qsm_tkd'};
% ptb.pipename = {'qsm_iswim'};
% ptb.pipename = {'qsm_ilsqr_stisuite22'};

% HELLO WORLD %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ptb.pipename = {'hello_world'};

ttt=cell2mat(ptb.pipename); if ~ttt; error('ERROR: Enable a ptb.pipename'); end; clear ttt

%% PROCESSING STEPS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if strcmpi(cell2mat(ptb.pipename),'gen_uncomb')
    disp(' '); disp('>>> Generic pipeline for uncombined GRE data')
    ptb.pipeline = {'prep.all'...
                    'comb.clapl_lbv_vsmv',...
                    'comb.complexfit',...
                    'unwr.clapl_cornell',...
                    'prep.unwr4d.srncp',...
                    'qsm.MSDI_170706'};

elseif strcmpi(cell2mat(ptb.pipename),'gen_comb')
    disp(' '); disp('>>> Generic pipeline for coil-combined GRE data')
    ptb.pipeline = {'prep.all'...
                    'comb.complexfit',...
                    'unwr.clapl_cornell',...
                    'prep.unwr4d.srncp',...
                    'bfe.lbv',...
                    'bfe.vsmv_stisuite22',...
                    'qsm.MSDI_170706'};
                
elseif strcmpi(cell2mat(ptb.pipename),'jac012_MSDI')
    disp(' '); disp('>>> Pipeline for multi-spectral QSM (MSDI v170706)')
    ptb.pipeline = {'prep.all'...
                    'prep.init_offset_corr',...
                    'prep.unwr4d.srncp',...
                    'comb.echofit',...
                    'bfe.lbv',...
                    'bfe.vsmv_stisuite22',...
                    'qsm.MSDI_170706'};

elseif strcmpi(cell2mat(ptb.pipename),'jac013_MSDI')
    disp(' '); disp('>>> Pipeline for multi-spectral QSM (MSDI v170706)')
    ptb.pipeline = {'prep.all'...
                    'prep.init_offset_corr',...
                    'prep.unwr4d.srncp',...
                    'comb.complexfit',...
                    'prep.unwr4d.srncp',...
                    'bfe.lbv',...
                    'qsm.MSDI_170706'};
                
elseif strcmpi(cell2mat(ptb.pipename),'qsm_MSDI170706')
    disp(' '); disp('>>> QSM inversion only (MSDI)')
    ptb.pipeline = {'prep.load_custom1',...
                    'prep.load_custom2',...
                    'prep.load_custom3',...
                    'prep.load_custom4',...
                    'prep.import_params',...
                    'qsm.MSDI_170706'};
                
elseif strcmpi(cell2mat(ptb.pipename),'qsm_fansi20170330')
    disp(' '); disp('>>> QSM inversion only (FANSI v20170330)')
    ptb.pipeline = {'prep.load_custom1',...
                    'prep.load_custom2',...
                    'prep.load_custom3',...
                    'prep.load_custom4',...
                    'prep.import_params',...
                    'qsm.fansi20170330'};

elseif strcmpi(cell2mat(ptb.pipename),'qsm_tkd')
    disp(' '); disp('>>> QSM inversion only (TKD)')
    ptb.pipeline = {'prep.load_custom1',...
                    'prep.load_custom2',...
                    'prep.load_custom3',...
                    'prep.load_custom4',...
                    'prep.import_params',...
                    'qsm.tkd'};

elseif strcmpi(cell2mat(ptb.pipename),'qsm_iswim')
    disp(' '); disp('>>> QSM inversion only (iSWIM)')
    ptb.pipeline = {'prep.load_custom1',...
                    'prep.load_custom2',...
                    'prep.load_custom3',...
                    'prep.load_custom4',...
                    'prep.import_params',...
                    'qsm.iswim'};
            
elseif strcmpi(cell2mat(ptb.pipename),'qsm_ilsqr_stisuite22')
    disp(' '); disp('>>> QSM inversion only (iLSQR)')
    ptb.pipeline = {'prep.load_custom1',...
                    'prep.load_custom2',...
                    'prep.load_custom3',...
                    'prep.load_custom4',...
                    'prep.import_params',...
                    'qsm.ilsqr_stisuite22'};
                
elseif strcmpi(cell2mat(ptb.pipename),'hello_world')
    disp(' '); disp('>>> Test run')
    ptb.pipeline = {'test.hw'};
end

% Print steps
for x=1:length(ptb.pipeline); disp([cell2mat(ptb.pipeline(x))]); end

%% LIST OF AVAILABLE METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prep.all  -------------------------  Run all preparation steps (see below) [recommended]
%  prep.load_default ----------------   Load 3D/4D magnitude & phase sets into magn[phase]_orig MATLAB variables. If magn[phase]_orig.nii* in pwd, they will be loaded automagically
%  prep.import_params  --------------   Import voxel size, echo time(s), field strength, field direction and phase scaling factor. 
%  prep.rescaling  ------------------   Rescale wrapped phases to the range: [-pi,pi] radians
%  prep.zeropad  --------------------   Upsample data to isotropic resolution via zero-padding (third dimension only)
%  prep.magn_comb  ------------------   Combine magnitude image as the root-sum-of-squares along magn_orig's 4th dimension
%  prep.roi  ------------------------   Estimate ROI mask. See CUSTOM SETTINGS below
%  prep.brainmask  ------------------  Import brain mask to constrain inversion consistency. See CUSTOM SETTINGS below
% prep.load_custom1[2..n]   ---------  As an alternative to prep.load_default, each prep.load_customN step will let you load any dataset into any variable of your choice (see CUSTOM SETTINGS below)
% prep.tvdenoise_magn  --------------  Total variation denoising of magnitude image - it could be generalised to take any input, but currently method will just load magn.nii/.gz with default denoising parameters
% prep.init_offset_corr -------------  Initial-phase offset correction (HF Sun). Separate odd/even echo adjustment for bipolar readouts
% prep.unwr4d.srncp  ----------------  Best-path unwrapping for 4D set and 2*pi jump temporal correction. Analogous to 'unwr.srncp' (below) - loop across 4th dim
% prep.bipolar_corr  ----------------  Phase correction for mGRE acquisitions using bipolar readouts (JQ Li, ECNU)
% comb.clapl_lbv  -------------------  Coil combination - magnitude-weighted filtered phase sum. Filters: continuous Laplacian + Laplacian boundary value [recommended]
% comb.clapl_lbv_vsmv  --------------  Coil combination - magnitude-weighted filtered phase sum. Filters: continuous Laplacian + Laplacian boundary value + variable spherical mean value (STI Suite v2.2)
% comb.dlapl_lbv  -------------------  Coil combination - magnitude-weighted filtered phase sum. Filters: discrete Laplacian + Laplacian boundary value [recommended]
% comb.dlapl_lbv_vsmv  --------------  Coil combination - magnitude-weighted filtered phase sum. Filters: discrete Laplacian + Laplacian boundary value + variable spherical mean value (STI Suite v2.2)
% comb.complexfit  ------------------  Phase estimation via magnitude-weighted nonlinear complex fitting (Cornell, JAC mod)
% comb.echofit  ---------------------  Phase estimation via magnitude-weighted LS phase fitting (HF Sun) [recommended]
% unwr.clapl_cornell  ---------------  Continuous Laplacian phase unwrapping (Cornell) [recommended]
% unwr.dlapl_bouwman  ---------------  Discrete Laplacian phase unwrapping (Bouwman)
% unwr.RG_cornell  ------------------  Region growing phase unwrapping (Cornell)
% unwr.PRELUDE  ---------------------  Phase region-expanding labeller for unwrapping discrete estimates (FSL's PRELUDE phase unwrapping)
% unwr.srncp  -----------------------  Fast 3D phase unwrapping based on sorting by reliability following a noncontinuous path (Abdul-Rahman)
% bfe.lbv  --------------------------  Full-multigrid Laplacian boundary value filtering [recommended]
% bfe.vsmv_stisuite22  --------------  Variable spherical mean value filtering (STI Suite v2.2) [recommended]
% bfe.vsmv_bilgic_bouwman  ----------  Variable spherical mean value filtering (Bilgic/Bouwman)
% bfe.esharp  -----------------------  Recover cortical phase through harmonic extension (SHARP edges, Topfer)
% bfe.poly3d  -----------------------  First-to-forth order polynomial fit to e.g. remove residual transmit phase (HF Sun)
% qsm.MSDI_170706  ------------------  Multi-spectral dipole inversion (MSDI) for QSM
% qsm.mwSB_Lsweep_bilgic_bouwman  ---  L1-regularised magnitude-weighted, Split-Bregman inversion with automated parameter selection
% qsm.fansi20170330  ----------------  Magnitude-constrained (vector field/L1-norm/L2-norm based), linear/nonlinear TV/TGV regularised inversion with continuous/discrete/integrated Green dipole kernel implementations (FANSI, C Milovic) 
% qsm.tkd  --------------------------  Truncated k-space division inversion method (K Shmueli)
% qsm.iswim  ------------------------  Susceptibility weighted imaging and mapping (M Haacke)
% qsm.ilsqr_stisuite22  -------------  Iterative LSQR inversion (STI Suite v2.2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CUSTOM SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% GENERAL OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Operating system
ptb.opt.OS = 2.1;
           % [1]   Windows
           % [2.1] Debian Linux
           % [2.2] Red Hat Linux
           
% Save all intermediate outputs
ptb.opt.save_interm = false;

% Save intermediate 4D files
ptb.opt.save_interm_ima4d = false;

%% PREP OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

ptb.iso_reslice = true;
                % if true, and slice thickness > 1.1*mean_inplane_vx_size,
                %  reslice (via zeropadding) to isotropic vx resolution

%% ROI MASK ESTIMATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Custom ROI mask
ptb.roi.custom.run = false;    
                   % default: roi.nii/.gz
                   %  if default mask exists in pwd, it will be autoloaded, 
                   %  otherwise it will prompt for manual input

% BET2, FSL's brain extraction tool
ptb.roi.bet2.run = true;
ptb.roi.bet2.f   = 0.1;     
                 % fractional threshold
                 
%% TRANSMIT PHASE CORRECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prep.init_offset_corr
ptb.readout = 'bipolar';
            % ['monopolar']

%% MULTI-ECHO PHASE FITTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% comb.echofit 
ptb.comb.echofit.intercept = 0; 
                           % zero/non-zero fitting intercept
                           % [0] if e.g. prep.init_offset_corr enabled
                           % [1] corrects for transmit/CF mismatch
                           
%% BACKGROUND FIELD EXTRACTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% LBV filtering
ptb.bfe.lbv.tol    = 0.001;
ptb.bfe.lbv.peel   = 3; 
ptb.bfe.lbv.depth  = -1; 
ptb.bfe.lbv.N1     = 30; 
ptb.bfe.lbv.N2     = 100; 
ptb.bfe.lbv.N3     = 100;

% vSMV filtering (STI Suite v2.2)
ptb.bfe.vsmv_stisuite22.R0 = 25;   
                           % initial convolution kernel radius in mm

% vSMV filtering (Bilgic/Bouwman)
ptb.bfe.vsmv_bilgic_bouwman.Rvector = 25:-1:1;
                                    % convolution kernel radii, vector in mm

% eSHARP (Topfer)
ptb.bfe.esharp.tik_reg = 5e-4;
                       % Tikhonov multiplier for reSHARP
ptb.bfe.esharp.smv_rad = 4;
                       % SMV convolution kernel radius in mm
ptb.bfe.esharp.radius  = [10,10,10];
                       % extent (ellipsoid semiaxes in mm) of background 
                       % field extrapolation

% Polynomial fit (Sun Hong Fu)
ptb.bfe.poly3d.poly_order = 4;
                      % order of polynomial function

%% PHASE RELIABILITY ESTIMATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prep.init_offset_corr | comb.complexfit
ptb.reliab_init_offset  = false;
                        % NOTE: untested feature
% prep.unwr4d.srncp 
ptb.reliab_unwr         = true;

% comb.echofit          | comb.complexfit
ptb.reliab_fitres       = false;

% prep.bipolar_corr 
ptb.reliab_bipolar      = false;
   
%% BRAIN MASK IMPORT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% prep.brainmask - morphological constrain to inversion consistency
ptb.flag_brainmask = false;
%                    [true]  msk_brain.nii/.gz (default), or will prompt
%                             for manual input
%                    [false] it will skip this step (import & masking)

%% QSM INVERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Multi-spectral dipole inversion v170706
ptb.qsm.MSDI_170706.lambda        = 750;
                                  % Lagrange multiplier
ptb.qsm.MSDI_170706.reliab        = 2;
                                  % Scale-wise consistency masking of unreliable phases
                                  %  [1] [R_smv(s)/R_smv(2)]x reliab_pct/100 
                                  %  [2] wC = 2^-log2(R_smv(s) [NO MASKING, eff. lambda scaling]
                                  %  [3] 1-2^-log2(R_smv(s)) [UNTESTED]
ptb.qsm.MSDI_170706.reliab_excl1  = true;
                                  % Do not apply masking for scale #1 (req. reliab=1 only)
ptb.qsm.MSDI_170706.reliab_pct    = 10;
                                  % Init extent, percentage masking (req. reliab=1 only)
ptb.qsm.MSDI_170706.merit         = 1;        
                                  % MEDI's dynamic model error reduction
                                  %  [0] disabled
                                  %  [1] enabled
ptb.qsm.MSDI_170706.smv           = true; % only enabled for version c+
ptb.qsm.MSDI_170706.R_smv         = [2,4,8,16];        
                                  % SMV kernel radii in vx
                                  %  R_smv could be a vector (multi-scale), 
                                  %  or a single scalar
ptb.qsm.MSDI_170706.e                          = 1e-6;
ptb.qsm.MSDI_170706.cg_max_iter                = 1000;
ptb.qsm.MSDI_170706.cg_tol                     = 0.1;
ptb.qsm.MSDI_170706.max_iter                   = 15;
ptb.qsm.MSDI_170706.tol_norm_ratio             = 0.1;
ptb.qsm.MSDI_170706.data_weighting_mode        = 1;
ptb.qsm.MSDI_170706.gradient_weighting_mode    = 3;
ptb.qsm.MSDI_170706.dynedge_mode               = 0;
                                  
% FANSI v20170330
ptb.qsm.fansi20170330.alpha                 = 2e-4;
ptb.qsm.fansi20170330.mu                    = 1e-2;
ptb.qsm.fansi20170330.options.nonlinear     = true;
                    % linear or nonlinear formulation?
ptb.qsm.fansi20170330.options.tgv           = false;
                    % TV or TGV regulariser?
ptb.qsm.fansi20170330.options.kernel_mode   = 2;
                    % [0] continuous kernel (Salomir et al. 2003)
                    % [1] discrete kernel (Milovic et al. 2017)
                    % [2] integrated Green function (Jenkinson et al. 2004)
ptb.qsm.fansi20170330.options.gradient_mode = 0;
                    % [0] vector field
                    % [1] L1 norm 
                    % [2] L2 norm

% TKD
ptb.qsm.tkd.thresh = 0.1;
                   % dipole kernel truncation level
                    
% iSWIM
ptb.qsm.iswim.thresh = 0.1;
                     % dipole kernel truncation level
ptb.qsm.iswim.erode  = 0;
                     % 1 mask erosion 

% iLSQR (STI Suite v2.2)
ptb.qsm.ilsqr_stisuite22.params.niter       = 30;
% ptb.qsm.ilsqr_stisuite22.params.padsize     = [0 0 0];
% ptb.qsm.ilsqr_stisuite22.params.cropsize    = cropsize;
% ptb.qsm.ilsqr_stisuite22.params.tol_step1   = tol_step1;
% ptb.qsm.ilsqr_stisuite22.params.tol_step2   = tol_step2;
% ptb.qsm.ilsqr_stisuite22.params.Kthreshold  = Kthreshold;
                   
%% LOAD CUSTOM IMAGES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Image #1
ptb.prep.load_custom1.froot       = 'unwr_clapl'; %'vsmv_stisuite22_R1-32';
                                  % NIFTI file root
ptb.prep.load_custom1.varname     = 'u_phase';
                                  % MATLAB variable name
ptb.prep.load_custom1.structname  = 'stp';
                                  % NIFTI structure (MATLAB variable name)

% Image #2
ptb.prep.load_custom2.froot       = 'magn';
ptb.prep.load_custom2.varname     = 'magn';
ptb.prep.load_custom2.structname  = 'stp';

% Image #3
ptb.prep.load_custom3.froot       = 'roi_comb'; %'msk_vsmv_stisuite22_R1-32';
ptb.prep.load_custom3.varname     = 'msk';
ptb.prep.load_custom3.structname  = 'stm';

% Image #4
ptb.prep.load_custom4.froot       = 'N_std';
ptb.prep.load_custom4.varname     = 'N_std';
ptb.prep.load_custom4.structname  = 'stp';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% SAVE SETTINGS
save('ptb.mat','ptb')

%% RUN PIPELINE
ptbm_001
