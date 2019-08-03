function ptbs_cse_003_def_nmedi_l400(flag)
% QSMbox
%  Default nMEDI pipeline for coil-combined single-GRE data
%
% USAGE 
%  + Edit CUSTOM SETTINGS (below) if necessary 
%  + Save & Run (F5)
% 
% REQUIREMENTS
%  + If 'ptb.roi.bet2.run' enabled, FSL installation required
% 
% DEFAULTS 
%  + Program will autoload the following files if in pwd:
%    - magn_orig.nii(.gz)              magnitude data           [3D_vol]
%    - phase_orig.nii(.gz)             phase data               [3D_vol]
% 
%  + If 'ptb.roi.custom.run' enabled the following file (in pwd) will skip manual input:
%    - roi.nii(.gz)                    ROI mask
%  
%  + Scan parameter input could also be automated if the following files were available in pwd:
%    - ptb_scaling_factor.txt          phase scaling,           e.g. 4096 (Siemens convention)
%    - ptb_res.txt                     voxel size in mm,        e.g. 1 1 1 
%    - ptb_TE.txt                      echo time(s) in s,       e.g. 0.02 (20 ms)
%    - ptb_B0.txt                      field strength in T,     e.g. 3
%    - ptb_B0_dir.txt                  main field direction,    e.g. 0 0 1
%    NOTE: The above -except scaling factor- can be autoextracted from a single DICOM file in a 'dicom' (small caps) dir (in pwd)
% 
%  + Dir structure example:
%     working_dir/magn_orig.nii(.gz)
%                /phase_orig.nii(.gz)
%                /ptb_scaling_factor.txt
%                /dicom/sample_dicom.dcm
%
%  + If any (or all) of the above missing you will be prompted for manual input
% 
% PIPELINE STEPS
%  1. prep.all  -------------------------  Preparation steps (see below)
%     - prep.load_default ---------------   Load 3D/4D magnitude & phase sets into magn[phase]_orig MATLAB variables. If magn[phase]_orig.nii* in pwd, they will be loaded automagically
%     - prep.import_params  -------------   Import voxel size, echo time(s), field strength, field direction and phase scaling factor 
%     - prep.rescaling  -----------------   Rescale wrapped phases to the range: [-pi,pi] rads
%     - prep.magn_comb  -----------------   Combine magnitude image as the root-sum-of-squares along magn_orig's 4th dimension
%     - prep.roi  -----------------------   Estimate/load ROI mask. See CUSTOM SETTINGS below
%  2. comb.complexfit  ------------------  [Transparent to single-echo data - used to ensure variable name consistency]
%  3. unwr.dlapl  -----------------------  Discrete Laplacian phase unwrapping (MEDI toolbox implementation) 
%  4. prep.unwr4d.srncp  ----------------  Second difference (phase reliability) calculation (based on code from Abdul-Rahman et al. Applied Optics 2007)
%  5. bfe.lbv  --------------------------  Full-multigrid Laplacian boundary value filtering (available from MEDI toolbox)
%  6. bfe.vsmv  -------------------------  Variable spherical mean value filtering (Bilgic/Bouwman implementation)
%  7. qsm.MSDI  -------------------------  Multi-scale dipole inversion (MSDI) for QSM (internally using a modified non-linear MEDI solver with SMV deconvolution)
% 
% Created by Julio Acosta-Cabronero

if nargin<1
    flag = 'run';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Uncomment below if you wish the present settings file to appear in the 
% list of PRESET OPTIONS.
%
% flag = 'preset';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% CUSTOM SETTINGS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% DEFINE PIPELINE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('Default nMEDI pipeline, L=400')
ptb.pipename = {'cse_003_def_nmedi_l400'}; 
ptb.pipeline = {'prep.all'...
                'comb.complexfit',...
                'unwr.dlapl',...
                'prep.unwr4d.srncp',...
                'bfe.lbv',...
                'bfe.vsmv',...
                'qsm.MSDI'};

%% PREP OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zero-padding
ptb.iso_reslice = false;
                % if [true], and slice thickness > 1.1*mean_inplane_vx_size,
                %  reslice (via zeropadding) to isotropic vx resolution

%% ROI MASK ESTIMATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Custom ROI mask
ptb.roi.custom.run = false;    
                   % [true] custom ROI mask (roi.nii/.gz in pwd) will be autoloaded 
                   % NOTE: if [true], but roi.nii/.gz not in pwd, it will prompt for manual input

% BET2, FSL's brain extraction tool (UNIX only)
ptb.roi.bet2.run = true;
ptb.roi.bet2.f   = 0.1;
                 % fractional threshold

%% BACKGROUND FIELD EXTRACTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LBV filtering
ptb.bfe.lbv.tol    = 0.001;
ptb.bfe.lbv.peel   = 2;
ptb.bfe.lbv.depth  = -1; 
ptb.bfe.lbv.N1     = 30; 
ptb.bfe.lbv.N2     = 100; 
ptb.bfe.lbv.N3     = 100;

% vSMV filtering
ptb.bfe.vsmv.R0 = 25;   
                % initial SMV kernel radius in mm

%% QSM INVERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% L1-regularised multi-scale inversion
ptb.qsm.MSDI.lambda        = 400;
                                  % Lagrange multiplier
ptb.qsm.MSDI.reliab        = 1;
                                  % Scale-wise consistency masking of unreliable phases
ptb.qsm.MSDI.reliab_excl1  = true;
                                  % Do not apply masking for scale #1   [if reliab==1]
ptb.qsm.MSDI.reliab_pct    = 10;
                                  % Init extent, percentage masking     [if reliab==1]
ptb.qsm.MSDI.merit         = 1;        
                                  % Dynamic model error reduction
ptb.qsm.MSDI.smv           = false; % <- No deconvolution, standard nMEDI
ptb.qsm.MSDI.R_smv         = [2];        
                                  % SMV kernel radii (in mm) defining each scale
                                  %  R_smv could be a vector (i.e. multi-scale QSM), 
                                  %  or a single scalar (i.e. MEDIN with SMV deconvolution)
ptb.qsm.MSDI.e                          = 1e-6;
                                               % L1-smoothing parameter
ptb.qsm.MSDI.cg_max_iter                = 2000;
                                               % Max iter # for internal conjugate gradient loop
ptb.qsm.MSDI.cg_tol                     = 0.1;
                                               % Stopping tolerance for CG loop
ptb.qsm.MSDI.max_iter                   = 25;
                                               % Max iter # for outer quasi-Newton loop
ptb.qsm.MSDI.tol_norm_ratio             = 0.1;
                                               % Stopping tolerance for outer qN loop
ptb.qsm.MSDI.data_weighting_mode        = 1;
                                               % Magnitude weighting (consistency)
ptb.qsm.MSDI.gradient_weighting_mode    = 3;
                                               % Edge masking (regulariser)
                                               %  [1] apply to all scales
                                               %  [2] apply to all scales - if reliab==1, prioritise consistency masking
                                               %                             i.e. M_edge=1, if M_reliab==0
                                               %  [3] apply to scale #1 only
ptb.qsm.MSDI.dynedge_mode               = 0;
                                               % Edge masking type
                                               % [0] Static:  magnitude driven (all scales)
                                               % [1] Dynamic: magnitude (scale #1) + QSM driven (subsequently) 

%% SAVE SETTINGS
if strcmpi(flag,'preset')
    [~,ptb.settings_file]   = fileparts(which(mfilename));
    ptb.settings_file       = [ptb.settings_file '.m'];
    save('ptb.mat','ptb')
end

if strcmpi(flag,'run')
    save('ptb.mat','ptb')
end

%% RUN PIPELINE
ptbm_001(flag)
