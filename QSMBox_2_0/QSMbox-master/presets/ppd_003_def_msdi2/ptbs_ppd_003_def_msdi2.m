function ptbs_ppd_003_def_msdi2(flag)
% QSMbox
%  Default MSDI v2 pipeline for pre-processed data
%
% USAGE 
%  + Edit CUSTOM SETTINGS (below) if necessary 
%  + Save & Run (F5)
% 
% REQUIREMENTS
%  + Data pre-processed for QSM
% 
% DEFAULTS 
%  + See LOAD CUSTOM IMAGES for details (edit if necessary)
%  
%  + Scan parameter input could also be automated if the following files were available in pwd:
%    - ptb_scaling_factor.txt          phase scaling,           e.g. 4096 (Siemens convention)
%    - ptb_res.txt                     voxel size in mm,        e.g. 1 1 1 
%    - ptb_TE.txt                      echo times in s,         e.g. 2.34e-3 4.64e-3 6.94e-3 9.24e-3 11.54e-3 13.84e-3 16.14e-3 18.44e-3
%    - ptb_B0.txt                      field strength in T,     e.g. 3
%    - ptb_B0_dir.txt                  main field direction,    e.g. 0 0 1 (slices perp to B0), or 0 0.5 0.866 (30deg rot)
%    NOTE: The above -except scaling factor & TEs- can be autoextracted from a single DICOM file in a 'dicom' (small caps) dir (in pwd)
% 
%  + Dir structure example:
%     working_dir/custom1_fname.nii(.gz)
%                /custom2_fname.nii(.gz)
%                /custom3_fname.nii(.gz)
%                /custom4_fname.nii(.gz)
%                /ptb_scaling_factor.txt
%                /ptb_TE.txt
%                /dicom/sample_dicom.dcm
%
%  + If any of the parameters are missing you will be prompted for manual input - not images though
% 
% PIPELINE STEPS
%  1. prep.load_custom1  -----------------  Load custom image #1 (pre-processed phase)
%  2. prep.load_custom2  -----------------  Load custom image #2 (magnitude)
%  3. prep.load_custom3  -----------------  Load custom image #3 (ROI mask)
%  4. prep.load_custom4  -----------------  Load custom image #4 (1/SNR)
%  5. prep.import_params  ----------------  Import voxel size, echo time(s), field strength, field direction and phase scaling factor 
%  6. qsm.MSDI  --------------------------  Multi-scale dipole inversion (MSDI) for QSM (see Acosta-Cabronero et al. Neuroimage 2018)
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
disp('Default MSDI v2 pipeline for pre-processed data')
ptb.pipename = {'ppd_003_def_msdi2'};
ptb.pipeline = {'prep.load_custom1',...
                'prep.load_custom2',...
                'prep.load_custom3',...
                'prep.load_custom4',...
                'prep.import_params',...
                'qsm.MSDI'};
            
%% LOAD CUSTOM IMAGES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Image #1
ptb.prep.load_custom1.froot       = 'unwr_3DSRNCP';
                                  % NIFTI file root
                                  % e.g. 'vsmv_R50-1vx';
ptb.prep.load_custom1.varname     = 'u_phase';
                                  % MATLAB variable name
ptb.prep.load_custom1.structname  = 'stp';
                                  % NIFTI structure (MATLAB variable name)

% Image #2
ptb.prep.load_custom2.froot       = 'magn';
ptb.prep.load_custom2.varname     = 'magn';
ptb.prep.load_custom2.structname  = 'stp';

% Image #3
ptb.prep.load_custom3.froot       = 'roi_comb';
                                  % 'msk_vsmv_R50-1vx';
                                  % 'Mask';
ptb.prep.load_custom3.varname     = 'msk';
ptb.prep.load_custom3.structname  = 'stm';

% Image #4
ptb.prep.load_custom4.froot       = 'N_std';
ptb.prep.load_custom4.varname     = 'N_std';
ptb.prep.load_custom4.structname  = 'stp';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% QSM INVERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Multi-scale dipole inversion v2 implementation
ptb.qsm.MSDI.smv            = true;
                            % If true, enable SMV deconvolution
ptb.qsm.MSDI.R_smv          = [-1,5];        
                            % SMV kernel radii (in mm) defining each scale
                            %  R_smv could be a vector (i.e. multi-scale QSM), 
                            %  or a single scalar (i.e. MEDIN with SMV deconvolution)
                            % N.B. If R_smv(1)==-1 & length(R_smv)>1: 
                            %  (i)    R_smv(1)=1 mm
                            %  (ii)   Optimise lambda for local scale
                            %  (iii)  Fine inversion with optimal lambda 
                            % This feature results in preservation of phase noise texture 
                            %  and low susceptibility contrast features
                            
% Self-optimised local scale            [if R_smv(1)==-1 & length(R_smv)>1] 
ptb.qsm.MSDI.search_L                   = logspace(4,5,20);
                                        % Sweep vector for lambda search
ptb.qsm.MSDI.search_full_sweep          = 0;
                                        % Sweep command
                                        %  [0] Stop search when optimisation criterion has been met
                                        %  [1] Full grid search
ptb.qsm.MSDI.search_cg_max_iter         = 5;
                                        % # of CG iterations for lambda search
ptb.qsm.MSDI.fine_cg_tol                = 1e-3;
                                        % Stopping tolerance for local scale's fine inversion

% Non-local solver                      [if R_smv(n)>0]
ptb.qsm.MSDI.lambda                     = round(10^2.6); % 10^2.6~398
                                        % Lagrange (consistency) multiplier
ptb.qsm.MSDI.e                          = 1e-6;
                                        % L1-smoothing parameter
ptb.qsm.MSDI.cg_max_iter                = 1000;
                                        % Max iter # for internal conjugate gradient loop
ptb.qsm.MSDI.cg_tol                     = 0.1;
                                        % Stopping tolerance for CG loop
ptb.qsm.MSDI.max_iter                   = 40;
                                        % Max iter # for outer quasi-Newton loop
ptb.qsm.MSDI.tol_norm_ratio             = 0.01;
                                        % Stopping tolerance for outer qN loop
                                        % N.B. If R_smv(n)>6 & tol_norm_ratio<0.1, the latter will be reset to 0.1
                                        
% Model adjustments
ptb.qsm.MSDI.data_weighting_mode        = 1;
                                        % Base noise model, magnitude weighting (consistency)
ptb.qsm.MSDI.smv_N_std_update           = true;
                                        % Adjust noise model as a function of R_smv

ptb.qsm.MSDI.reliab                     = 1;
                                        % Scale-wise consistency masking of unreliable phases
ptb.qsm.MSDI.reliab_excl1               = true;
                                        % Do not apply masking for scale #1   [if reliab==1]
ptb.qsm.MSDI.reliab_pct                 = 10;
                                        % Init extent, percentage masking     [if reliab==1]

ptb.qsm.MSDI.merit                      = 1;        
                                        % Dynamic model error reduction
ptb.qsm.MSDI.merit_f                    = 6;        
                                        % MERIT's threshold
ptb.qsm.MSDI.merit_p                    = 2;        
                                        % MERIT's noise model corrective power
                                                                                
ptb.qsm.MSDI.gradient_weighting_mode    = 2;
                                        % Edge masking (regulariser)
                                        %  [1] apply to all scales
                                        %  [2] apply to all scales   [if reliab==1, prioritise consistency masking]
                                        %                             i.e. M_edge=1, if M_reliab==0
                                        %  [3] apply to scale #1 only
ptb.qsm.MSDI.dynedge_mode               = 3;
                                        % Edge masking type
                                        %  [0] Static I:     magnitude, magn.nii(.gz)
                                        %  [1] Dynamic I:    magnitude (scale #1) + QSM (subsequently)
                                        %  [2] Dynamic II:   scale-dependent SMV filtered field map [if smv]
                                        %  [3] Static II:    input field map, nfm_0.nii(.gz)

%% SAVE SETTINGS
flag = m1_save_ptb(ptb,flag);

%% RUN PIPELINE
ptbm_001(flag)
