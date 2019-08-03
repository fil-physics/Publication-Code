function ptbs_ume_001_step3(flag)
% QSMbox
%  Default pipeline for coil-uncombined multi-GRE data (STEP #3)
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
%    - magn_orig.nii(.gz)              magnitude data           4D, [3D_vols,coils]
%    - phase_orig.nii(.gz)             phase data               "Idem"
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
%  2. comb.dlapl_lbv_vsmv  --------------  Coil combination - magnitude-square weighted filtered phase sum. Filters: discrete Laplacian + Laplacian boundary value (MEDI toolbox) + variable spherical mean value (Bilgic/Bouwman implementation)
%  3. comb.complexfit  ------------------  Pass-on data (rename only step)
%  4. unwr.dlapl  -----------------------  Discrete Laplacian phase unwrapping
%  5. prep.unwr4d.srncp  ----------------  Second difference (phase reliability) calculation (based on code from Abdul-Rahman et al. Applied Optics 2007)
%  6. qsm.MSDI  -------------------------  Multi-scale dipole inversion (MSDI) for QSM (see Acosta-Cabronero et al. Neuroimage 2018)
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
disp('Default pipeline for uncombined multi-echo data (STEP #3)')
ptb.pipename = {'ume_001_step3'}; 
ptb.pipeline = {'prep.all'...
                'comb.dlapl_lbv_vsmv',...
                'comb.complexfit',...
                'unwr.dlapl',...
                'prep.unwr4d.srncp',...
                'qsm.MSDI'};
            
%% PREP OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Zero-padding
ptb.iso_reslice = false; % N.B. IF PADDING REQUIRED, IT MUST BE ENABLED IN STEPS 1 & 2
                % if [true], and slice thickness > 1.1*mean_inplane_vx_size,
                %  reslice (via zeropadding) to isotropic vx resolution

%% ROI MASK ESTIMATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Custom ROI mask
ptb.roi.custom.run = true;    
                   % [true] custom ROI mask (roi.nii/.gz in pwd) will be autoloaded 
                   % NOTE: if [true], but roi.nii/.gz not in pwd, it will prompt for manual input

%% BACKGROUND FIELD EXTRACTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LBV filtering
ptb.bfe.lbv.tol    = 0.001;
ptb.bfe.lbv.peel   = 0;
ptb.bfe.lbv.depth  = -1; 
ptb.bfe.lbv.N1     = 30; 
ptb.bfe.lbv.N2     = 100; 
ptb.bfe.lbv.N3     = 100;

% vSMV filtering
ptb.bfe.vsmv.R0 = 40;   
                % initial SMV kernel radius in mm

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
ptb.qsm.MSDI.tol_norm_ratio             = 0.1;
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
