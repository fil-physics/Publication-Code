function ptbs_ppd_002_def_cme_unwr(flag)
% QSMbox
%  Default phase preprocessing pipeline for coil-combined multi-GRE data
%
% USAGE 
%  + Edit CUSTOM SETTINGS (below) if necessary 
%  + Save & Run (F5)
% 
% REQUIREMENTS
%  + Select readout type, i.e. monopolar or bipolar
%  + If 'ptb.roi.bet2.run' enabled, FSL installation required
% 
% DEFAULTS 
%  + Program will autoload the following files if in pwd:
%    - magn_orig.nii(.gz)              magnitude data           4D, [3D_vols,echo]
%    - phase_orig.nii(.gz)             phase data               "Idem"
% 
%  + If 'ptb.roi.custom.run' enabled the following file (in pwd) will skip manual input:
%    - roi.nii(.gz)                    ROI mask
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
%     working_dir/magn_orig.nii(.gz)
%                /phase_orig.nii(.gz)
%                /ptb_scaling_factor.txt
%                /ptb_TE.txt
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
%  2. prep.init_offset_corr -------------  Transmit-phase offset correction and (optional) odd/even echo adjustment for bipolar readouts (based on public-domain code from Hongfu Sun)
%  3. prep.unwr4d.srncp  ----------------  Echo-wise best-path unwrapping with second difference (phase reliability) calculation and 2*pi jump temporal correction (based on code from Abdul-Rahman et al. Applied Optics 2007 available from HF Sun's GitHub repo)
%  4. comb.echofit  ---------------------  Phase estimation via WLS phase fitting (adapted from HF Sun implementation)
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
disp('Default phase preprocessing pipeline for coil-combined multi-echo data')
disp('(initial phase-offset correction + best-path unwrapping + echofit)')
ptb.pipename = {'ppd_002_def_cme_unwr'};
ptb.pipeline = {'prep.all'...
                'prep.init_offset_corr',...
                'prep.unwr4d.srncp',...
                'comb.echofit'};

%% PREP OPTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Intermediate outputs
ptb.opt.save_interm       = true;
ptb.opt.save_interm_ima4d = true;

% Zero-padding
ptb.iso_reslice = true;
                % if [true], and slice thickness > 1.1*mean_inplane_vx_size,
                %  reslice (via zeropadding) to isotropic vx resolution

%% ROI MASK ESTIMATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Custom ROI mask
ptb.roi.custom.run = false;    
                   % [true] custom ROI mask (roi.nii/.gz in pwd) will be autoloaded 
                   % NOTE: if [true], but roi.nii/.gz not in pwd, it will prompt for manual input

% SPM, https://www.fil.ion.ucl.ac.uk/spm
ptb.roi.spm.run       = true;
ptb.roi.spm.settings  = '3T_default';
                      % alternative preset: '7T_default';
ptb.roi.spm.thr       = 1e-4;    
                      % c1+c2+c3 probability cut-off threshold
ptb.roi.spm.dilR      = -1;    
                      % ROI dilation (spherical) kernel radius in mm 
                      %  [0]  no dilation
                      %  [-1] dilation by one voxel
                                            
% BET2, FSL's brain extraction tool (UNIX only)
ptb.roi.bet2.run = false;
ptb.roi.bet2.f   = 0.1;
                 % fractional threshold

%% PHASE CORRECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% prep.init_offset_corr
ptb.readout = 'bipolar';
            % Multi-echo readout type
            % ['monopolar']
            % ['bipolar']

% comb.echofit 
ptb.comb.echofit.intercept = 0; 
                           % zero/non-zero fitting intercept
                           % [0] if e.g. prep.init_offset_corr enabled
                           % [1] corrects for transmit/CF mismatch
            
%% SAVE SETTINGS
flag = m1_save_ptb(ptb,flag);

%% RUN PIPELINE
ptbm_001(flag)
