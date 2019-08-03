function ptbm_001(flag)
% QSMbox MASTER PROGRAM, DO NOT EDIT
% Created by Julio Acosta-Cabronero

%% INITIALISATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(flag,'run')
    load ptb.mat
    
% Save new settings as preset option
elseif strcmpi(flag,'preset')
    load ptb.mat    
    ptbdir      = fileparts(which(mfilename));
    presetsdir  = [ptbdir '/../presets'];
    pipename    = cell2mat(ptb.pipename);
    
    if exist([presetsdir '/' pipename],'dir')
        disp(' ')
        disp(['Warning: ''' pipename ''' is an existing preset'])
        disp(' ')
        presetow = input('Would you like to overwrite it? [y/n] ','s');
        switch presetow
            case 'y'
                disp(' ')
                disp(['Export to: presets/' pipename '/'])
                disp(' ')
                new_settings_file = ['ptbs_' pipename '.m'];
                fid = fopen(new_settings_file,'r'); f = fread(fid,'*char')'; fclose(fid);
                f   = strrep(f,'flag = ''preset''','% flag = ''preset'''); 
                fid = fopen(new_settings_file,'w'); fprintf(fid,'%s',f); fclose(fid);
                copyfile(new_settings_file,[presetsdir '/' pipename '/' new_settings_file]);
                edit(['ptbs_' pipename])               
            case 'n'
                return
            otherwise
                error('ERROR: Wrong selection. Start again.')
        end        
    else
        disp(' ')
        disp(['Export to presets/' pipename '/'])
        disp(' ')
        new_settings_file = ['ptbs_' pipename '.m'];
        copyfile(ptb.settings_file,new_settings_file)
        fid = fopen(new_settings_file,'r'); f = fread(fid,'*char')'; fclose(fid);
        f   = strrep(f,'flag = ''preset''','% flag = ''preset'''); 
        fid = fopen(new_settings_file,'w'); fprintf(fid,'%s',f); fclose(fid);
        mkdir([presetsdir '/' pipename])
        copyfile(new_settings_file,[presetsdir '/' pipename '/' new_settings_file]);
        edit(['ptbs_' pipename])               
    end
    disp(['Press F5 to run new pipeline'])
    disp(' ')
    return

else
    return
end

% Print steps
disp(' ')
for x=1:length(ptb.pipeline); disp([cell2mat(ptb.pipeline(x))]); end

% Start timer
ptb.clock = clock;

%% GLOBAL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Operating system
if exist('~/.QSMbox/ptb_OS.m','file')
    run('~/.QSMbox/ptb_OS.m');
    ptb.opt.OS = ans;
else
    if ~exist('ptb.opt.OS','var')
        disp(' ')
        disp('The Operating System setting file ''~/.QSMbox/ptb_OS.m'' could not be found.')
        disp('This is a requirement for best-path unwrapping.')
        disp('If required, QSMbox will blindly use an Ubuntu 16.04 binary.')
        disp('If this is incorrect, please edit and run:')
        disp(' ''QSMbox/bashutils/QSMbox_OS.sh''')
        ptb.opt.OS = 2.10;
    end
end
% Intermediate outputs
if ~isfield(ptb,'opt')
    ptb.opt = [];
end
if ~isfield(ptb.opt,'save_interm')
    % Save all intermediate outputs
    ptb.opt.save_interm       = false;
end
if ~isfield(ptb.opt,'save_interm_ima4d')
    % Save intermediate 4D files
    ptb.opt.save_interm_ima4d = false;
end
if ~isfield(ptb.opt,'gadmpm_flag')
    % MPM-Gadgetron input
    ptb.opt.gadmpm_flag       = false;
end
% prep.brainmask - morphological constrain to inversion consistency
ptb.flag_brainmask          = false;
%                           [true]  msk_brain.nii/.gz (default), or will prompt
%                                    for manual input
%                           [false] it will skip this step (import & masking)
% prep.init_offset_corr | comb.complexfit
ptb.reliab_init_offset      = false;
                            % N.B.: untested feature
% comb.echofit | comb.complexfit
ptb.reliab_fitres           = false;
% prep.bipolar_corr 
ptb.reliab_bipolar          = false;
% prep.unwr4d.srncp 
ptb.reliab_unwr             = true;
% comb
if ~isfield(ptb,'comb')
    ptb.comb = [];
end
if ~isfield(ptb.comb,'keep_uncombined')
    ptb.comb.keep_uncombined = false;
end
% MSDI
if ~isfield(ptb,'qsm')
    ptb.qsm                     = [];
end
if ~isfield(ptb.qsm,'MSDI')
    ptb.qsm.MSDI                = [];
end
% Noise model adjustment due to SMV prefiltering 
if isfield(ptb.qsm.MSDI,'smv_N_std_0')
    ptb.qsm.MSDI.smv_N_std_update = 1-ptb.qsm.MSDI.smv_N_std_0;
end
if ~isfield(ptb.qsm.MSDI,'smv_N_std_update')
    ptb.qsm.MSDI.smv_N_std_update = true;
end
% MERIT arguments
if ~isfield(ptb.qsm.MSDI,'merit_f')
    ptb.qsm.MSDI.merit_f        = 6;
end
if ~isfield(ptb.qsm.MSDI,'merit_p')
    ptb.qsm.MSDI.merit_p        = 2;
end
% Phase reliability weighting
if ~isfield(ptb.qsm.MSDI,'reliab')
    ptb.qsm.MSDI.reliab         = 1;
end
if ~isfield(ptb.qsm.MSDI,'reliab_excl1')
    ptb.qsm.MSDI.reliab_excl1   = true;
end
if ~isfield(ptb.qsm.MSDI,'reliab_pct')
    ptb.qsm.MSDI.reliab_pct     = 10;
end
% Fidelity boost
if ~isfield(ptb.qsm.MSDI,'search_L')
    ptb.qsm.MSDI.search_L               = logspace(4,5,20);
end
if ~isfield(ptb.qsm.MSDI,'search_full_sweep')
    ptb.qsm.MSDI.search_full_sweep      = 0;
end
if ~isfield(ptb.qsm.MSDI,'search_cg_max_iter')
    ptb.qsm.MSDI.search_cg_max_iter     = 5;
end
if ~isfield(ptb.qsm.MSDI,'fine_cg_tol')
    ptb.qsm.MSDI.fine_cg_tol            = 1e-3;
end
if ~isfield(ptb.qsm.MSDI,'fidelity_boost')
    ptb.qsm.MSDI.fidelity_boost         = false;
end
if ~isfield(ptb.qsm.MSDI,'fidelity_boost_factor')
    ptb.qsm.MSDI.fidelity_boost_factor  = 10^2.4;
end
if ~isfield(ptb.qsm.MSDI,'fidelity_boost_maxiter')
    ptb.qsm.MSDI.fidelity_boost_maxiter = 1;
end
if ~isfield(ptb,'roi')
    ptb.roi                 = [];
end
if ~isfield(ptb.roi,'spm')
    ptb.roi.spm.run         = false;
end
if ~isfield(ptb.roi.spm,'settings')
    ptb.roi.spm.settings    = '3T_default';
end
if ~isfield(ptb.roi.spm,'thr')
    ptb.roi.spm.thr         = 1e-4;
end
if ~isfield(ptb.roi.spm,'dilR')
    ptb.roi.spm.dilR        = -1; % ROI dilation by 1 vx
end
% Other flags
reslice_flag        = false;
downsampl_factor    = 1;
gadmpm_flag         = ptb.opt.gadmpm_flag;

%% INIT PIPELINE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add generic toolboxes to MATLAB path
ptbdir = fileparts(which(mfilename));
addpath([ptbdir '/ptb/_PhaseTools'])
addpath([ptbdir '/ptb/_Utils'])
addpath([ptbdir '/ptb/_spm12b_nifti'])

% Run through steps
for stepnumber = 1:length(ptb.pipeline)
    stepname   = cell2mat(ptb.pipeline(stepnumber));

%% HELLO WORLD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(stepname,'test.hw')
    disp(' '); disp('>>> Hello World')
end
    
%% DATA (DEFAULT) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(stepname,'prep.all')...
|| strcmpi(stepname,'prep.load_default')

    disp(' '); disp('>>> Load data')
    
    % Magnitude
    dir0     = pwd;
    datadirM = false;
    if exist('magn_orig.nii') || exist('magn_orig.nii.gz')
        disp('Magnitude data found')
        [magn_froot] = 'magn_orig';
        
    elseif exist('echo1_Magnitude.nii') || exist('echo1_Magnitude.nii.gz')
        disp('Magnitude data found')
        [magn_froot] = 'echo1_Magnitude';
        
    else
        [fname,datadirM] = uigetfile({'*.nii.gz';'*.nii'},...
                                      'Select magnitude data');
        [magn_froot] = m1_fname_split(fname,'gz');
        cd(datadirM);
        [~,ext] = m1_fname_split(fname);
        if strcmpi(ext,'gz') 
            copyfile([datadirM '/' fname],[dir0 '/magn_orig.nii.gz'])
        else
            copyfile([datadirM '/' fname],[dir0 '/magn_orig.nii'])
        end
    end

    % Phase
    datadirP = false;
    if exist('phase_orig.nii') || exist('phase_orig.nii.gz')
        disp('Phase data found')
        [phase_froot] = 'phase_orig';
        if datadirM
            if exist('phase_orig.nii')
                disp('Import phase_orig.nii from the same directory')
                copyfile([datadirM '/phase_orig.nii'],[dir0 '/phase_orig.nii'])
            elseif exist('phase_orig.nii.gz')
                disp('Import phase_orig.nii.gz from the same directory')
                copyfile([datadirM '/phase_orig.nii.gz'],[dir0 '/phase_orig.nii.gz'])
            end
        end
        
    elseif exist('echo1_Phase.nii') || exist('echo1_Phase.nii.gz')
        disp('Phase data found')
        [phase_froot] = 'echo1_Phase';
        
    else
        [fname,datadirP] = uigetfile({'*.nii.gz';'*.nii'},...
                                      'Select phase data' );
        [phase_froot] = m1_fname_split(fname,'gz');
        cd(datadirP);
        [~,ext] = m1_fname_split(fname);
        if strcmpi(ext,'gz')
            copyfile([datadirP '/' fname],[dir0 '/phase_orig.nii.gz'])
        else
            copyfile([datadirP '/' fname],[dir0 '/phase_orig.nii'])
        end
    end

    % Load data
    tic %
    cd(dir0);
    if strcmpi(magn_froot,'echo1_Magnitude') && strcmpi(phase_froot,'echo1_Phase')
        gadmpm_flag = true;
        % Magnitude
        filelist    = dir('echo*_Magnitude.nii*');
        NE          = length(filelist);
        for x=1:NE
           eval(['magn_froot_temp = ''echo' num2str(x) '_Magnitude'';'])
            [magn_orig_temp,stm] = m1_nifti_load(magn_froot_temp);
           if gadmpm_flag; ttt = magn_orig_temp; ttt=permute(ttt,[3 1 2]); ttt=ttt(:,end:-1:1,end:-1:1); magn_orig_temp = ttt; clear ttt; end
           magn_orig(:,:,:,x) = magn_orig_temp;
        end
         clear magn_orig_temp
        
        % Phase 
        filelist    = dir('echo*_Phase.nii*');
        NEM         = NE;
        NE          = length(filelist);
        if NE~=NEM
            error('ERROR: Number of Magnitude and Phase images do not match. Exiting.')
        end
        for x=1:NE
           eval(['phase_froot_temp = ''echo' num2str(x) '_Phase'';'])
            [phase_orig_temp,stp] = m1_nifti_load(phase_froot_temp);
           if gadmpm_flag; ttt = phase_orig_temp; ttt=permute(ttt,[3 1 2]); ttt=ttt(:,end:-1:1,end:-1:1); phase_orig_temp = ttt; clear ttt; end
           phase_orig(:,:,:,x) = phase_orig_temp;
        end
         clear phase_orig_temp
        
    else
        if datadirM; cd(datadirM); end
         [magn_orig,stm]  = m1_nifti_load(magn_froot);
        if datadirP; cd(datadirP); end
         [phase_orig,stp] = m1_nifti_load(phase_froot);
    end
    stm(1).dt(1) = 4; % NIFTI structure: force integer (used for saving binary masks)
    cd(dir0);
    toc %
end

%% DATA (CUSTOM) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(stepname,'prep.load_custom1')
    disp(' '); disp('>>> Custom data import #1')
    eval(['[ttt,' ptb.prep.load_custom1.structname...
          '] = m1_nifti_load(''' ptb.prep.load_custom1.froot ''');'])
    if gadmpm_flag; ttt=permute(ttt,[3 1 2]); ttt=ttt(:,end:-1:1,end:-1:1); end
    eval([ptb.prep.load_custom1.varname ' = ttt;'])
end

if strcmpi(stepname,'prep.load_custom2')
    disp(' '); disp('>>> Custom data import #2')
    eval(['[ttt,' ptb.prep.load_custom2.structname...
          '] = m1_nifti_load(''' ptb.prep.load_custom2.froot ''');'])
    if gadmpm_flag; ttt=permute(ttt,[3 1 2]); ttt=ttt(:,end:-1:1,end:-1:1); end
    eval([ptb.prep.load_custom2.varname ' = ttt;'])
end

if strcmpi(stepname,'prep.load_custom3')
    disp(' '); disp('>>> Custom data import #3')
    eval(['[ttt,' ptb.prep.load_custom3.structname...
          '] = m1_nifti_load(''' ptb.prep.load_custom3.froot ''');'])
    if gadmpm_flag; ttt=permute(ttt,[3 1 2]); ttt=ttt(:,end:-1:1,end:-1:1); end
    eval([ptb.prep.load_custom3.varname ' = ttt;'])
end

if strcmpi(stepname,'prep.load_custom4')
    disp(' '); disp('>>> Custom data import #4')
    eval(['[ttt,' ptb.prep.load_custom4.structname...
          '] = m1_nifti_load(''' ptb.prep.load_custom4.froot ''');'])
    if gadmpm_flag; ttt=permute(ttt,[3 1 2]); ttt=ttt(:,end:-1:1,end:-1:1); end
    eval([ptb.prep.load_custom4.varname ' = ttt;'])
end

clear ttt

%% SCAN PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(stepname,'prep.all')... 
|| strcmpi(stepname,'prep.import_params')

    disp(' '); disp('>>> Import scan parameters')

    % INIT: voxel size, field strength, echo time and phase scaling factor
    res=false; B0=false; TE=false; phase_scaling_factor=false;

    % If they exist, read from text files
    if exist('ptb_res.txt')
        disp('Reading voxel size from ptb_res.txt')
        res = double(load('ptb_res.txt'));
    end
    if strcmpi(stepname,'prep.import_params')
        if exist('ptb_res_new.txt')
            disp('It appears your data was resliced - reading voxel size from ptb_res_new.txt instead')
            res = double(load('ptb_res_new.txt'));
        end
    end
    if exist('ptb_B0.txt')
        disp('Reading field strength from ptb_B0.txt')
        B0 = double(load('ptb_B0.txt'));
    end
    if exist('ptb_TE.txt')
        disp('Reading echo time(s) from ptb_TE.txt')
        TE = double(load('ptb_TE.txt'));
    end
    if exist('ptb_scaling_factor.txt')
        disp('Reading phase scaling factor from ptb_scaling_factor.txt')
        phase_scaling_factor = double(load('ptb_scaling_factor.txt'));
    end

    % Alternatively, read from DICOM file in 'dicom' dir
    if exist('dicom')
        if ~res
            disp('Reading voxel size from DICOM')
            [res,~,~] = m1_dicom_extract_res_B0_TE('dicom');
        end
        if ~B0
            disp('Reading field strength from DICOM')
            [~,B0,~] = m1_dicom_extract_res_B0_TE('dicom');
        end
        if ~TE
            if size(magn_orig,4)>1 && strcmpi(cell2mat(ptb.pipename),'gen_comb')
                disp('4D coil-combined data - all TEs cannot be extracted from a single DICOM file')
                TE = input('Input echo times (vector in s, e.g. [5:5:20]*1e-3): ');
            end
            disp('Reading echo time(s) from DICOM')
            [~,~,TE] = m1_dicom_extract_res_B0_TE('dicom');     
        end

    % Otherwise prompt for manual input
    else
        if ~res
            disp('Voxel size information not found')
            res = input('Input voxel resolution (vector in mm, e.g. [1 1 2]): ');
        end
        if ~B0
            disp('Field strength information not found')
            B0 = input('Input field strength (in T): ');
        end
        if ~TE
            disp('Echo time information not found')
            TE = input('Input echo time (if multi-GRE, vector in s, e.g. [5:5:20]*1e-3): ');
        end
    end
    
    % Sanity check
    if strcmpi(cell2mat(ptb.pipename),'gen_comb')
        if size(magn_orig,4) ~= length(TE)
            error('ERROR: Number of echo times do not match number of input volumes. Exiting.')
        end
    end

    % Read B0_dir from text, otherwise manual input
    if ~exist('ptb_B0_dir.txt')
        disp('B0 direction information not found')
        B0_dir = input('Input main field direction (normalised vector, e.g. [0 0 1] if axial slices perpendicular to main field): ');
        B0_dir = double(B0_dir); 
         save('ptb_B0_dir.txt','B0_dir','-ascii')
    else
        disp('Reading B0 direction from ptb_B0_dir.txt')
        B0_dir = double(load('ptb_B0_dir.txt'));
    end    

    % If scaling factor file not provided, prompt for manual input
    if ~phase_scaling_factor
        disp(' ')
        disp('Phase must be scaled to [-pi,pi] radians')
        disp('Siemens for example writes out phase data in the [-4096,4096] range')
        disp('Thus, phase_scaling_factor = 4096 would be required')
        disp(' ')
        phase_scaling_factor = input(['Input phase scaling factor (INFO: native phase range: [' num2str(min(phase_orig(:))) ',' num2str(max(phase_orig(:))) ']). If phase = [-pi,pi] rad, type pi: ']);
        phase_scaling_factor = double(phase_scaling_factor); 
         save('ptb_scaling_factor.txt','phase_scaling_factor','-ascii')
    end

    % Sanity check
    if ~res
        error('ERROR: Program is missing some scan parameters. Try again.')
    elseif ~B0
        error('ERROR: Program is missing some scan parameters. Try again.')
    elseif ~TE
        error('ERROR: Program is missing some scan parameters. Try again.')

    % Print and save scan parameters
    else
        disp(' '); disp('>>> Scan parameters successfully imported:')
        disp(['Voxel size, mm     ' num2str(res)])
        disp(['Echo time(s), s    ' num2str(TE)])
        disp(['Field strength, T  ' num2str(B0)])
        disp(['B0 direction       ' num2str(B0_dir)])
        disp(['Phase scaling      ' num2str(phase_scaling_factor)])
        if strcmpi(cell2mat(ptb.pipename),'gen_uncomb')
            disp(['Coil elements      ' num2str(size(magn_orig,4))])
        end
        res = double(res); 
         save('ptb_res.txt','res','-ascii')
        TE  = double(TE);  
         save('ptb_TE.txt','TE','-ascii')
        B0  = double(B0);  
         save('ptb_B0.txt','B0','-ascii')
    end
end

%% PHASE RESCALING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(stepname,'prep.all')...
|| strcmpi(stepname,'prep.rescaling')

    if phase_scaling_factor<3.14159 || phase_scaling_factor>3.1416
        disp(' '); disp('>>> Rescale phase to [-pi,pi] rad')
        [phase_orig] = double((phase_orig*pi)/phase_scaling_factor);
        % NIFTI structure: int => float
        stp(1).dt(1)=16;
    end
end

%% ZEROPAD TO ISOTROPIC RESOLUTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(stepname,'prep.all')...
|| strcmpi(stepname,'prep.zeropad')

    if (res(3)/mean([res(1),res(2)])>1.1 && ptb.iso_reslice)
        tic %
        reslice_flag = 1;
        disp(' '); disp('>>> Upsample to isotropic resolution (cosine filtering + zero-padding)')
        disp('N.B.: Depending on matrix size, this might require a large amount of memory')
        disp('      If you hit issues, disable ptb.iso_reslice in the ptbs file')
        res0                    = res;
        full_size               = size(magn_orig);
        matrix_size             = full_size(1:3);
        [dim1,dim2,dim3,dim4]   = size(magn_orig);
        res_inplane_mean        = mean([res0(1),res0(2)]);
        res_slice_factor        = res0(3)/res_inplane_mean;
        pad_slice               = dim3*res_slice_factor-dim3;
        pad                     = [0,0,round(pad_slice/2)];
        kfdoff                  = 10; % # of voxels where the filter will drop off from 1 to 0.
        [kf1,kf2,kf3]           = ndgrid(ones([dim1 1]),ones([dim2 1]),m1_cosfilt(dim3,kfdoff));
        kf                      = min(min(kf1,kf2),kf3); clear kf1 kf2 kf3
        for x=1:dim4
            S = magn_orig(:,:,:,x).*exp(1i*phase_orig(:,:,:,x));
            A=ifftshift(S);                         clear S
            B=fftn(A);                              clear A
            C=fftshift(B);                          clear B
            K = 1/sqrt(length(C(:)))*C;             clear C
            K_real = kf.*real(K);
            K_real_pad = padarray(K_real, pad);     clear K_real
            K_imag = kf.*imag(K);                   clear K
            K_imag_pad = padarray(K_imag, pad);     clear K_imag
            K_pad = complex(K_real_pad,K_imag_pad); clear K_*_pad
            A=ifftshift(K_pad);                     clear K_pad
            B=ifftn(A);                             clear A
            C=fftshift(B);                          clear B
            S_pad = sqrt(length(C(:)))*C;           clear C
            magn_orig_pad(:,:,:,x)  = abs(S_pad);
            phase_orig_pad(:,:,:,x) = angle(S_pad); clear S_pad
            if x==1
                matrix_size     = size(magn_orig_pad);
                res(3)          = (res0(3)*dim3)/(dim3+2*pad(3));
                save('ptb_res_new.txt','res','-ascii')
            end
        end
        magn_orig   =  magn_orig_pad;               clear magn_orig_pad kf
        phase_orig  =  phase_orig_pad;              clear phase_orig_pad
        toc %
    end
end

%% ENSURE EVEN NUMBER OF SLICES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(stepname,'prep.all')
    
    [~,~,dim3,dim4] = size(magn_orig);
    if mod(dim3,2)~=0
        tic %
        reslice_flag = 2;
        disp(' '); disp('>>> Odd number of slices detected - zeropad to make it even')
        res0                    = res;
        pad                     = [0,0,1];
        for x=1:dim4
            S = magn_orig(:,:,:,x).*exp(1i*phase_orig(:,:,:,x));
            A=ifftshift(S);                                 clear S
            B=fftn(A);                                      clear A
            C=fftshift(B);                                  clear B
            K = 1/sqrt(length(C(:)))*C;                     clear C
            K_real = real(K);
            K_real_pad = padarray(K_real, pad, 'post');     clear K_real
            K_imag = imag(K);                               clear K
            K_imag_pad = padarray(K_imag, pad, 'post');     clear K_imag
            K_pad = complex(K_real_pad,K_imag_pad);         clear K_*_pad
            A=ifftshift(K_pad);                             clear K_pad
            B=ifftn(A);                                     clear A
            C=fftshift(B);                                  clear B
            S_pad = sqrt(length(C(:)))*C;                   clear C
            magn_orig_pad(:,:,:,x)  = abs(S_pad);
            phase_orig_pad(:,:,:,x) = angle(S_pad);         clear S_pad
            if x==1
                res(3) = (res0(3)*dim3)/(dim3+1);
                save('ptb_res_new.txt','res','-ascii')
            end
        end
        magn_orig   =  magn_orig_pad;               clear magn_orig_pad
        phase_orig  =  phase_orig_pad;              clear phase_orig_pad
        toc %
    end
end

%% UPDATE NIFTI HEADER %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(stepname,'prep.all')
    
    if reslice_flag==1
        disp(' '); disp('>>> Update NIFTI header to reflect resolution change')
        stp = stp(1);
        % matrix dimensions
        full_size               = size(magn_orig);
        matrix_size             = full_size(1:3);
        stp.dim                 = matrix_size;
        stp.private.dat.dim     = matrix_size;
        % work out signs and remove rotations
        for y1 = 1:3
            res_hdr             = stp.mat(y1,y1);
            sgn_res_hdr(y1)     = abs(res_hdr)/res_hdr;
            for y2 = 1:3
                stp.mat(y1,y2)  = 0;
            end
        end
        for y3 = 1:3
            % new voxel size
            stp.mat(y3,y3)      = sgn_res_hdr(y3)*res(y3);
            % % reset origin to centre of FOV
            % stp.mat(y3,4)       = -sgn_res_hdr(y3)*((matrix_size(y3)*res(y3))/2);
        end
        stp.private.mat         = stp.mat;
        stp.private.mat0        = stp.mat;
        stm                     = stp;
        stm.dt(1)               = 4;
    end
    if reslice_flag
        disp(' ')
        disp(['New matrix dimensions:   ' num2str(matrix_size)])
        disp(['New voxel size, mm:      ' num2str(res)])
    end
end
    
%% MAGNITUDE COMBINATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(stepname,'prep.all')...
|| strcmpi(stepname,'prep.magn_comb')

    if size(magn_orig,4)
        if length(TE)==1
            disp(' '); disp('>>> 4-D magnitude image (single-TE) detected')
            disp('RSoS combination')
            % RSoS combination
            [magn] = sqrt(sum(magn_orig.^2,4));
            ttt = magn; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
             m1_nifti_save(ttt,stp,'magn'); clear ttt
        else
            disp(' '); disp('>>> 4-D magnitude image (multi-TE) detected')
            disp('Synthesise: magn=S0.*exp(-mean(TE).*R2star)')
            % Synthesise magnitude image for mean TE
            [magn] = m1_meanmagn(magn_orig,TE);
            ttt = magn; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
             m1_nifti_save(ttt,stp,'magn'); clear ttt
        end
    else
        [magn] = magn_orig;
    end
end
    
%% ROI MASK ESTIMATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(stepname,'prep.all')...
|| strcmpi(stepname,'prep.roi')
    tic %
    disp(' '); disp('>>> Define ROI')

    % EXTERNAL INPUT
    if ptb.roi.custom.run
        disp('Do not estimate, import predefined mask')
        dir0=pwd; 
        datadirR=false;
        if ~exist('roi.nii') && ~exist('roi.nii.gz')
            [fname,datadirR] = uigetfile({'*.nii.gz';'*.nii'},...
                                          'Select ROI mask');
            [roi_froot] = m1_fname_split(fname,'gz');
            cd(datadirR);
            [~,ext] = m1_fname_split(fname);
            if strcmpi(ext,'gz') 
                copyfile([datadirR '/' fname],[dir0 '/roi.nii.gz'])
            else
                copyfile([datadirR '/' fname],[dir0 '/roi.nii'])
            end
        else
            [roi_froot] = 'roi';
        end
        if datadirR; cd(datadirR); end
        [msk] = m1_nifti_load(roi_froot);
        if gadmpm_flag; ttt = magn; ttt=permute(ttt,[3 1 2]); ttt=ttt(:,end:-1:1,end:-1:1); magn = ttt; clear ttt; end
        cd(dir0);
        if ptb.roi.bet2.run
            disp('BET2 brain mask calculation also enabled. Prioritise external input.') 
            ptb.roi.bet2.run = false;
        end
        if ptb.roi.spm.run
            disp('SPM-based brain mask calculation also enabled. Prioritise external input.') 
            ptb.roi.spm.run = false;
        end
    end
    
    % BET2
    if ptb.roi.bet2.run
        F = num2str(ptb.roi.bet2.f);
        disp(['Brain mask estimation with BET2, fractional threshold: ' F])
        disp(['N.B.: You need a UNIX-based OS and FSL to run this.']) 
        disp(['      If the program throws an error at this point, it means BET2 is not available.']) 
        disp(['      You can also create a brain mask using SPM.'])
        disp(['      Set SPM path in ~/.QSMbox/ptb_spmpath.m.'])
        disp(['      See QSMbox/bashutils/QSMbox_local_settings.sh.'])
        system(['bet2 magn bet2-f' F '_magn -f ' F]);
        [roi_froot] =    ['bet2-f' F '_magn'];
        fname0 = dir([roi_froot '.nii*']);
        fname  = fname0(1).name;
        [~,ext] = m1_fname_split(fname);
        if strcmpi(ext,'gz') 
            copyfile(fname,'roi.nii.gz')
        else
            copyfile(fname,'roi.nii')
        end        
         [msk] = m1_nifti_load(roi_froot);
        if gadmpm_flag; ttt = msk; ttt=permute(ttt,[3 1 2]); ttt=ttt(:,end:-1:1,end:-1:1); msk = ttt; clear ttt; end
    end

    % SPM
    if ptb.roi.spm.run
        disp(['Brain mask estimation with SPM (settings: ' ptb.roi.spm.settings ', probability threshold: ' num2str(ptb.roi.spm.thr) ')'])
        if exist('~/.QSMbox/ptb_spmpath.m','file')
            run('~/.QSMbox/ptb_spmpath.m');
        else
            disp(' ')
            disp('ERROR: The SPM path file ''~/.QSMbox/ptb_spmpath.m'' could not be found.')
            disp('This is a requirement if you wish to use SPM for brain mask calculation.')
            disp('Set SPM path in ~/.QSMbox/ptb_spmpath.m.')
            disp('See QSMbox/bashutils/QSMbox_local_settings.sh.')                
            disp('N.B. Alternatively, you can use BET2 (see ptbs file).')
            disp('Exiting...')
            return %%%
        end
        
        d0=pwd;
        d1='SPM';
        mkdir(d1)
        
        fname0 = dir('magn.nii*');
        fname  = fname0(1).name;
        [~,ext] = m1_fname_split(fname);
        if strcmpi(ext,'gz') 
            copyfile(fname,[d1 '/magn.nii.gz'])
            gunzip([d1 '/magn.nii.gz']);
            delete([d1 '/magn.nii.gz']);
        else
            copyfile(fname,[d1 '/magn.nii'])
        end             
        
        cd(d1)
        
        m1_spmsegment('magn.nii',ptb.roi.spm.settings);     % SPM segmentation
        
        disp(' '); disp('ROI calculation')
        m1_spmroi(ptb.roi.spm.thr,ptb.roi.spm.dilR,res);    % ROI calculation
        
        cd(d0)
        
        if exist([d1 '/roi.nii.gz'])
            copyfile([d1 '/roi.nii.gz'],'roi.nii.gz')
        elseif exist([d1 '/roi.nii'])
            copyfile([d1 '/roi.nii'],'roi.nii')
        else
            error('ERROR: Something went wrong - roi.nii(.gz) was not created. Exiting.')
        end
         [msk] = m1_nifti_load('roi');
        if gadmpm_flag; ttt = msk; ttt=permute(ttt,[3 1 2]); ttt=ttt(:,end:-1:1,end:-1:1); msk = ttt; clear ttt; end
    end
    
    % Force binary
    msk(msk>0) = 1;
    msk_orig   = msk;
    toc %
end

%% BRAIN MASK IMPORT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(stepname,'prep.all')...
|| strcmpi(stepname,'prep.brainmask')
    
    if ptb.flag_brainmask
        % brain mask for morphological constrain to inversion consistency
        disp(' '); disp('>>> Import brainmask')
        if exist('msk_brain.nii') || exist('msk_brain.nii.gz')
            [msk_brain_froot] = 'msk_brain';
             [msk_brain] = m1_nifti_load(msk_brain_froot);
            if gadmpm_flag; ttt = msk_brain; ttt=permute(ttt,[3 1 2]); ttt=ttt(:,end:-1:1,end:-1:1); msk_brain = ttt; clear ttt; end
            msk_brain(msk_brain>0) = 1;
        else
            dir0=pwd;
            datadirB=false;
            [fname,datadirB] = uigetfile({'*.nii.gz';'*.nii'},...
                                      'Select brain mask');
            [msk_brain_froot] = m1_fname_split(fname,'gz');
            if datadirB; cd(datadirB); end
             [msk_brain] = m1_nifti_load(msk_brain_froot);
            if gadmpm_flag; ttt = msk_brain; ttt=permute(ttt,[3 1 2]); ttt=ttt(:,end:-1:1,end:-1:1); msk_brain = ttt; clear ttt; end
            msk_brain(msk_brain>0) = 1;
            cd(dir0);
        end
    end    
end
    
%% COIL COMBINATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(stepname,'comb.dlapl_lbv')...
|| strcmpi(stepname,'comb.dlapl_lbv_vsmv')...
|| strcmpi(stepname,'comb.dlapl_lbv_poly3d')
    tic %
    disp(' ')
    if ~ptb.comb.keep_uncombined
        disp('>>> Laplacian-based coil combination')
    else
        disp('>>> Echo-wise prefiltering')
    end
    addpath([ptbdir '/ptb/_LBV'])

    % Init loop, memory efficient
    msk_comb = ones(size(msk)); 
    if ~ptb.comb.keep_uncombined
        S_comb = zeros(size(msk));
    end
       
    for coil = 1:size(magn_orig,4)
        disp(' ')
        if ~ptb.comb.keep_uncombined
            disp(['>> Channel #' num2str(coil)])
        else
            disp(['>> Echo #' num2str(coil)])
        end
        
        [magn_ttt]    = magn_orig(:,:,:,coil);
        [u_phase_ttt] = phase_orig(:,:,:,coil);

        disp(['DL unwrapping'])
        [u_phase_ttt] = unwrapLaplacian(u_phase_ttt,size(u_phase_ttt),res);

        disp(['LBV filtering'])
        [u_phase_ttt] = LBV(u_phase_ttt,msk,size(u_phase_ttt),res,...
                     ptb.bfe.lbv.tol,ptb.bfe.lbv.depth,ptb.bfe.lbv.peel,...
                     ptb.bfe.lbv.N1,ptb.bfe.lbv.N2,ptb.bfe.lbv.N3);
        if ptb.bfe.lbv.peel>0
            maskID  = fopen(['mask_p' num2str(ptb.bfe.lbv.peel) '.bin']);           
            msk = fread(maskID,'single'); fclose(maskID);            
            msk = reshape(msk,size(u_phase_ttt));  
            msk(msk>0)=1;
        end
        toc %
        
        if strcmpi(stepname,'comb.dlapl_lbv_vsmv')
           disp(['vSMV filtering'])
           if coil==1; disp('N.B.: If this step throws an error, try reducing ptb.bfe.vsmv.R0 in CUSTOM SETTINGS.'); end
           R0 = round(ptb.bfe.vsmv.R0/mean(res));
           Rvector = R0:-1:1;
           [~,u_phase_ttt,msk] = m1_bb_vSHARP(u_phase_ttt,msk,Rvector);
        elseif strcmpi(stepname,'comb.dlapl_lbv_poly3d')
           [u_phase_ttt] = m1_shf_poly3d(u_phase_ttt,msk,ptb.bfe.poly3d.poly_order); 
        end
        
        msk_comb = msk_comb.*msk;
        if ~ptb.comb.keep_uncombined
            S_comb = S_comb + msk.*(magn_ttt.*exp(1i*u_phase_ttt)).^2;
        else
            phase_orig(:,:,:,coil) = u_phase_ttt;
        end
        msk = msk_orig;
    end
    
    clear *_ttt
    if ~ptb.comb.keep_uncombined 
        clear magn_orig phase_orig 
    end
     
    msk = msk_comb;
    roi_froot = 'roi_comb';
    ttt = msk; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
     m1_nifti_save(ttt,stm,roi_froot); clear ttt
    
    if ~ptb.comb.keep_uncombined 
        magn_orig  = msk.*magn;
        phase_orig = msk.*angle(sqrt(S_comb));
         clear *_comb 
    end
end

%% INITIAL PHASE and EVEN/ODD OFFSET CORRECTION for MULTI-ECHO DATA %%%%%%%
if strcmpi(stepname,'prep.init_offset_corr')
    tic %
    addpath([ptbdir '/ptb/_3DSRNCP'])

    % monopolar readout(s)
    if strcmpi(ptb.readout,'monopolar')
        disp(' '); disp('>>> Init-offset correction (monopolar readouts)')
        [phase_orig,unph_diff,phase_offsets] = m1_shf_mgre_comb(magn_orig.*exp(1i*phase_orig),res,TE,msk,ptb.opt.OS);
        ttt = unph_diff; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
         m1_nifti_save(ttt,    stp,'phase_diff');
        ttt = phase_offsets; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
         m1_nifti_save(ttt,stp,'phase_offsets'); clear ttt
         clear unph_diff
         
    % bipolar readouts
    elseif strcmpi(ptb.readout,'bipolar')
        disp(' '); disp('>>> Init-offset and bipolar readout correction')
        disp('Adjust odd echoes')
        phase_ttt                               = phase_orig(:,:,:,1:2:end);
        [phase_ttt,unph_diff_odd,offsets_odd]   = m1_shf_mgre_comb(magn_orig(:,:,:,1:2:end).*exp(1i*phase_orig(:,:,:,1:2:end)),res,TE(1:2:end),msk,ptb.opt.OS);
        phase_orig(:,:,:,1:2:end)               = phase_ttt;
        disp('Adjust even echoes')
        phase_ttt                               = phase_orig(:,:,:,2:2:end);
        [phase_ttt,unph_diff_even,offsets_even] = m1_shf_mgre_comb(magn_orig(:,:,:,2:2:end).*exp(1i*phase_orig(:,:,:,2:2:end)),res,TE(2:2:end),msk,ptb.opt.OS);
        phase_orig(:,:,:,2:2:end)               = phase_ttt; clear phase_ttt 
        ttt = unph_diff_odd; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
         m1_nifti_save(ttt,  stp, 'phase_diff_odd');
        ttt = offsets_odd; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
         m1_nifti_save(ttt,    stp, 'phase_offsets_odd');
        ttt = unph_diff_even; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
         m1_nifti_save(ttt, stp, 'phase_diff_even');
        ttt = offsets_even; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
         m1_nifti_save(ttt,   stp, 'phase_offsets_even'); clear ttt
        phase_offsets = (abs(offsets_odd)+abs(offsets_even))/2;
         clear *_odd *_even
    
    else
        error('ERROR: Incorrect setting for readout polarity. Set ptb.readout in ptbs file to ''monopolar'' or ''bipolar''.')
    end
    
    % Save dependencies
    Debug = false;
    if Debug && ~exist('ptb_temp01.mat')
        disp('Debug mode: saving temporary .mat')
        save ptb_temp01.mat
    end
    
    if ptb.reliab_init_offset
        % Untested
        disp('Reliability estimation (transmit phase offset)')
        [~,phase_offsets_edge] = m1_SMVfilt(phase_offsets,res,1);
        R_initOffset = msk./phase_offsets_edge; R_initOffset(isnan(R_initOffset))=0; R_initOffset(isinf(R_initOffset))=0;
        ttt = R_initOffset; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
         m1_nifti_save(ttt,stp,'R_initOffset');
         clear R_initOffset ttt
    end
     clear phase_offsets*
     
    if ptb.opt.save_interm
        disp('Save 3D NIFTIs')
        outroot = 'initCorr_phase';
        for x=1:size(phase_orig,4)
            ttt = phase_orig(:,:,:,x); if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
             m1_nifti_save(ttt,stp,[outroot '_vol' num2str(10+x)]);
        end
         clear ttt
        if ptb.opt.save_interm_ima4d && ptb.opt.OS>2 % UNIX only
            disp('Merge 4D output')
            system(['fslmerge -t all_' outroot ' ' outroot '_vol*']);
            delete([outroot '_vol*'])
        end
    end
    toc %
end

%% SPATIAL PHASE UNWRAPPING ACROSS ECHOES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(stepname,'prep.unwr4d.srncp')
    tic %
    disp(' '); disp('>>> 3D best-path phase unwrapping')
    disp(['N.B.: If phase unwrapping fails, take a look at: ' ptbdir '/ptb/_3DSRNCP/3DSRNCP_jac_README.txt'])
    addpath([ptbdir '/ptb/_3DSRNCP'])
    if exist('u_phase','var'); phase_orig = u_phase; flag_u_phase = true; else flag_u_phase = false; end
    [phase_orig,uncert_unwr] = m1_shf_3DSRNCP(phase_orig,msk,ptb.opt.OS);
    if ptb.reliab_unwr && ptb.qsm.MSDI.reliab
        disp('Reliability estimation (3D second difference quality map)')
        mean_uncert_unwr = sum(uncert_unwr,4)/size(uncert_unwr,4);
        R_unwr = sqrt(mean_uncert_unwr); R_unwr(isnan(R_unwr))=0; R_unwr(isinf(R_unwr))=0;
        ttt = R_unwr; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
         m1_nifti_save(ttt,stp,'secdiff_phase');
         clear R_unwr ttt
    end
     clear *uncert_unwr
    
    if size(phase_orig,4)>1
        disp('Search for 2*pi jumps across echoes')
        [~,unph_diff,~] = m1_shf_mgre_comb(magn_orig.*exp(1i*phase_orig),res,TE,msk,ptb.opt.OS);
        flag_nojump = true;
        for x = 2:size(phase_orig,4)
            meandiff = phase_orig(:,:,:,x)-phase_orig(:,:,:,1)-double(x-1)*unph_diff;
            meandiff = meandiff(msk==1);
            meandiff = mean(meandiff(:));
            njump = round(meandiff/(2*pi));
            if njump>0; disp([num2str(njump) ' 2*pi jumps for echo #' num2str(x)]); flag_nojump = false; end
            phase_orig(:,:,:,x) = phase_orig(:,:,:,x)-njump*2*pi;
            phase_orig(:,:,:,x) = phase_orig(:,:,:,x).*msk;
        end
         clear meandiff njump unph_diff offsets
        if flag_nojump; disp('No jumps found'); end
    else
        disp('Save unwrapped phase')
        outroot  = 'unwr_3DSRNCP';
        ttt = phase_orig; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
         m1_nifti_save(ttt,stp,outroot); clear ttt
    end
    
    if ptb.opt.save_interm
        disp('Save 3D NIFTIs')
        outroot  = 'unwr_3DSRNCP';
        for x=1:size(phase_orig,4)
            ttt = phase_orig(:,:,:,x); if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
             m1_nifti_save(ttt,stp,[outroot  '_vol' num2str(10+x)]);
        end
         clear ttt
        if ptb.opt.save_interm_ima4d && ptb.opt.OS>2 % UNIX only
            disp('Merge 4D output')
            system(['fslmerge -t all_' outroot ' ' outroot '_vol*']);   % u_phase
            delete([outroot '_vol*'])
        end
    end
    
    if flag_u_phase; u_phase = phase_orig; clear phase_orig; end
    toc %
end

%% PHASE ESTIMATION (LINEAR FIT) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(stepname,'comb.echofit')
    stepname = 'comb.echofit';
    tic %
    disp(' '); disp('>>> Temporal phase fitting')
    if ptb.comb.echofit.intercept==0
        disp('Set to null fitting intercept; it assumes previously corrected for transmit phase offset.')
        [u_phase,fitres] = m1_shf_echofit(phase_orig,magn_orig,TE,0);
        phase_offsets = 1;
    else
        [u_phase,fitres,phase_offsets] = m1_shf_echofit(phase_orig,magn_orig,TE,ptb.comb.echofit.intercept);
        ttt = phase_offsets; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
         m1_nifti_save(ttt,stp,'echofit_offset'); clear ttt
    end
    
    if ptb.reliab_init_offset && ptb.comb.echofit.intercept>0
        disp('Reliability estimation (transmit phase offset)')
        [~,phase_offsets_edge] = m1_SMVfilt(phase_offsets,res,1);
        R_initOffset = msk./phase_offsets_edge; R_initOffset(isnan(R_initOffset))=0; R_initOffset(isinf(R_initOffset))=0;
        ttt = R_initOffset; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
         m1_nifti_save(ttt,stp,'R_initOffset');
         clear R_initOffset ttt
    end
    
     clear magn_orig phase_orig phase_offsets*     
    
    % Save NIFTIs
    ttt = u_phase; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
     m1_nifti_save(ttt,stp,'echofit_phase');
    ttt = fitres; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
     m1_nifti_save(ttt, stp,'echofit_residual'); clear ttt
    
    N_std = msk.*(1./magn); N_std(isnan(N_std))=0; N_std(isinf(N_std))=0;
    ttt = N_std; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
     m1_nifti_save(ttt,stp,'N_std'); clear ttt
    
    if ptb.reliab_fitres
        disp('Reliability estimation (phase fitting residual)')
        R_fitres = sqrt(fitres); R_fitres(isnan(R_fitres))=0; R_fitres(isinf(R_fitres))=0;
        ttt = R_fitres; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
         m1_nifti_save(ttt,stp,'R_fitres');  
         clear R_fitres ttt
    end
     clear fitres
    toc %
end

%% ALTERNATIVE BIPOLAR READOUT CORRECTION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(stepname,'prep.bipolar_corr')
    % Based on Li et al. Phase-corrected bipolar gradients in multi-echo 
    % gradient-echo sequences for quantitative susceptibility mapping. 
    % MAGMA,2014 DOI:10.1007/s10334-014-0470-3
    tic %
    disp(' '); disp('>>> Bipolar-readout phase correction')
    if size(phase_orig,4)>=3
        disp('> 3D best-path phase unwrap first three echoes')
        addpath([ptbdir '/ptb/_3DSRNCP'])
        phase_ttt = phase_orig(:,:,:,1:3);
        [phase_ttt] = m1_shf_3DSRNCP(phase_ttt,msk,ptb.opt.OS);
        
        disp('> Calculate errors')
        [phase_offsets] = 1.0*phase_ttt(:,:,:,2)-...
                          0.5*phase_ttt(:,:,:,1)-...
                          0.5*phase_ttt(:,:,:,3);
        ttt = phase_offsets; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end             
         m1_nifti_save(ttt, stp, 'bipolar_phase_error');
         clear phase_ttt ttt

        if ptb.reliab_bipolar
            disp('Reliability estimation (bipolar offset)')
            [~,phase_offsets_edge] = m1_SMVfilt(phase_offsets,res,1);
            R_bipolarOffset = msk./phase_offsets_edge; R_bipolarOffset(isnan(R_bipolarOffset))=0; R_bipolarOffset(isinf(R_bipolarOffset))=0;
            ttt = R_bipolarOffset; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
             m1_nifti_save(ttt,stp,'R_bipolarOffset');
             clear R_bipolarOffset ttt
        end
         clear phase_offsets    

        disp('> Fit polynomial')
         poly_order = 1;
        [error_nopoly] = m1_shf_poly3d(phase_offsets,msk,poly_order);
        [error_poly]   = phase_offsets-error_nopoly;
        ttt = error_nopoly; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
         m1_nifti_save(ttt, stp, 'bipolar_error_nopoly');
        ttt = error_poly; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
         m1_nifti_save(ttt,   stp, 'bipolar_error_poly'); clear ttt
        
        disp('> Apply phase adjustments to even echoes')
        c=0;
        for x=1:size(phase_orig,4)
            c=c+1;
            if c==2
                phase_orig(:,:,:,x) = phase_orig(:,:,:,x)-error_poly;
                c=0;
            end
        end
         clear *error*
        
    else
        error('ERROR: Multi-GRE data with at least three echoes required. Try again.')
    end

    if ptb.opt.save_interm
        disp('Save 3D NIFTIs')
        outroot = 'bipolarCorr_phase';
        for x=1:size(phase_orig,4)
            ttt = phase_orig(:,:,:,x); if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
             m1_nifti_save(ttt,stp,[outroot '_vol' num2str(10+x)]);
        end
         clear ttt
        if ptb.opt.save_interm.ima4d && ptb.opt.OS>2 % UNIX only
            disp('Merge 4D output')
            system(['fslmerge -t all_' outroot ' ' outroot '_vol*']);
            delete([outroot '_vol*'])
        end
    end
    toc %    
end

%% PHASE ESTIMATION (NONLINEAR FIT) %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(stepname,'comb.complexfit') 
    tic %
    disp(' '); disp('>>> Complex fitting')
    single_echo = false;
    
    if length(TE)==1
        disp('Single GRE data. No fitting will be carried out - only variable preparation.')
        single_echo = true;
%         [u_phase,N_std,relres,phase_0] = ptb_Fit_ppm_complex(S_complex,TE);
        [u_phase]   = phase_orig;
        [delta_TE]  = TE;
    else
        [S_complex] = magn_orig.*exp(-1i*phase_orig);
    end
     clear magn_orig phase_orig
    
    if length(TE)==2
        disp('Dual GRE data')
        [u_phase,N_std,relres,phase_0] = ptb_Fit_ppm_complex(S_complex,TE);
        [delta_TE] = TE(2)-TE(1);

    elseif length(TE)>2
        if abs((TE(2)-TE(1))-(TE(3)-TE(2)))<0.2e-3
            disp('Multi GRE (evenly spaced) data')
            [u_phase,N_std,relres,phase_0] = ptb_Fit_ppm_complex(S_complex,TE);
        else
            disp('Multi GRE (unevenly spaced) data')
            [u_phase,N_std,relres,phase_0] = ptb_Fit_ppm_complex_TE(S_complex,TE);
        end
        [delta_TE] = TE(2)-TE(1);
    end
     clear S_complex
    
    if ptb.reliab_init_offset && strcmpi(ptb.readout,'monopolar')
        disp('Reliability estimation (initial offset)')
        [~,phase_offsets_edge] = m1_SMVfilt(phase_0,res,1);
        R_initOffset = msk./phase_offsets_edge; R_initOffset(isnan(R_initOffset))=0; R_initOffset(isinf(R_initOffset))=0;
        ttt = R_initOffset; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
         m1_nifti_save(ttt,stp,'R_initOffset');
         clear R_initOffset ttt
    end
     clear phase_offsets*    
    
    if ptb.reliab_fitres
        disp('Reliability estimation (fitting residual)')
        R_fitres = 1./relres; R_fitres(isnan(R_fitres))=0; R_fitres(isinf(R_fitres))=0;
        ttt = R_fitres; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
         m1_nifti_save(ttt,stp,'R_fitres');  
         clear R_fitres ttt
    end
         
    if ~single_echo
        ttt = u_phase; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
         m1_nifti_save(ttt,  stp, 'phase_fit');
        ttt = phase_0; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
         m1_nifti_save(ttt,  stp, 'phase_0');
        ttt = N_std; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
         m1_nifti_save(ttt,    stp, 'N_std');
        ttt = 1./N_std; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
         m1_nifti_save(ttt, stp, 'snr_magn');
        ttt = relres; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
         m1_nifti_save(ttt,   stp, 'relres_fit');
         clear relres phase_0 ttt
    end
    toc %
end

%% PHASE UNWRAPPING METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(stepname,'unwr.dlapl')
    tic %
    disp(' '); disp('>>> Discrete Laplacian phase unwrapping')
    [u_phase] = unwrapLaplacian(u_phase,size(u_phase),res);
    ttt = u_phase; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
     m1_nifti_save(ttt,stp,'unwr_dlapl'); clear ttt
    toc %
end
    
if strcmpi(stepname,'unwr.dlapl_bouwman')
    tic %
    disp(' '); disp('>>> Discrete Laplacian phase unwrapping') 
    [u_phase] = m1_bb_discretelaplacianUnwrap(u_phase,ones(size(u_phase)));
    ttt = u_phase; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
     m1_nifti_save(ttt,stp,'unwr_dlapl'); clear ttt
    toc %
end

if strcmpi(stepname,'unwr.RG_cornell')
    tic %
    disp(' '); disp('>>> Region growing phase unwrapping') 
    [u_phase] = unwrapPhase(magn,u_phase,size(u_phase));
    ttt = u_phase; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
     m1_nifti_save(ttt,stp,'unwr_RG'); clear ttt
    toc %
end

if strcmpi(stepname,'unwr.prelude')
    tic %
    disp(' '); disp('>>> PRELUDE phase unwrapping')
    ttt = msk.*u_phase; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
     m1_nifti_save(ttt,stp,'roi_wrapped_phase');
    ttt = msk.*magn; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
     m1_nifti_save(ttt,   stp,'roi_magn'); clear ttt
     if ptb.opt.OS>2
         system(['prelude -a roi_magn'...
                       ' -p roi_wrapped_phase'...
                       ' -u unwr_PRELUDE'...
                       ' -m ' roi_froot...
                       ' -n 12']);
         [u_phase] = m1_nifti_load('unwr_PRELUDE');
        if gadmpm_flag; ttt = u_phase; ttt=permute(ttt,[3 1 2]); ttt=ttt(:,end:-1:1,end:-1:1); u_phase = ttt; clear ttt; end
        toc %
     else
         error('ERROR: FSL PRELUDE not available for your operating system. Exiting.')
     end
end

if strcmpi(stepname,'unwr.srncp')
    tic %
    disp(' '); disp('>>> Best-path phase unwrapping')
    disp(['N.B.: If phase unwrapping fails, take a look at: ' ptbdir '/ptb/_3DSRNCP/3DSRNCP_jac_README.txt'])
    addpath([ptbdir '/ptb/_3DSRNCP'])
    [u_phase,msk] = m1_shf_3DSRNCP(u_phase,msk,ptb.opt.OS);
    ttt = u_phase; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
     m1_nifti_save(ttt,stp,'unwr_3DSRNCP');
    ttt = msk; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
     m1_nifti_save(ttt,stm,'msk_reliability_3DSRNCP'); clear ttt
    toc %
end

%% BACKGROUND FIELD EXTRACTION METHODS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmpi(stepname,'bfe.lbv')  
    stepname = 'bfe.lbv';
    tic %
    disp(' '); disp('>>> Full-multigrid Laplacian boundary value filter')
    addpath([ptbdir '/ptb/_LBV'])
    [u_phase] = LBV(u_phase,msk,size(u_phase),res,...
                     ptb.bfe.lbv.tol,ptb.bfe.lbv.depth,ptb.bfe.lbv.peel,...
                      ptb.bfe.lbv.N1,ptb.bfe.lbv.N2,ptb.bfe.lbv.N3);
    if ptb.bfe.lbv.peel>0
        maskID  = fopen(['mask_p' num2str(ptb.bfe.lbv.peel) '.bin']);           
        msk = fread(maskID,'single'); fclose(maskID);       
        msk = reshape(msk,size(u_phase));  
        msk(msk>0)=1;
        ttt = msk; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
         m1_nifti_save(ttt,stm,'msk_lbv');
    end
    ttt = u_phase; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
     m1_nifti_save(ttt,stp,'lbv'); clear ttt
    toc %
end

if strcmpi(stepname,'bfe.vsmv_stisuite22')
    tic %
    disp(' '); disp('>>> Variable spherical mean value filter')
    addpath([ptbdir '/ptb/_STISUITE22'])
    disp('N.B.: If this step throws an error, try reducing ptb.bfe.vsmv_stisuite22.R0 in CUSTOM SETTINGS')
    R0 = round(ptb.bfe.vsmv.R0/mean(res));
    [u_phase,msk]=V_SHARP(u_phase,msk,'voxelsize',res,'smvsize',R0);
    ttt = u_phase; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
     m1_nifti_save(ttt,stp,['vsmv_stisuite22_R1-' num2str(R0)]);
    ttt = msk; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
     m1_nifti_save(ttt,stm,['msk_vsmv_stisuite22_R1-' num2str(R0)]); clear ttt
    toc %
end

if strcmpi(stepname,'bfe.vsmv_stisuite142')
    tic %
    disp(' '); disp('>>> Variable spherical mean value filter')
    addpath([ptbdir '/ptb/_STISUITE22'])
    R0 = round(ptb.bfe.vsmv.R0/mean(res));
    padsize=[100 100 100];
    [u_phase,msk]=V_SHARP_142(u_phase,msk,R0,padsize,res);
    ttt = u_phase; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
     m1_nifti_save(ttt,stp,['vsmv_stisuite142_R1-' num2str(R0)]);
    ttt = msk; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
     m1_nifti_save(ttt,stm,['msk_vsmv_stisuite142_R1-' num2str(R0)]); clear ttt
    toc %
end

if strcmpi(stepname,'bfe.vsmv')
    tic %
    disp(' '); disp('>>> Variable spherical mean value filter')
    R0 = round(ptb.bfe.vsmv.R0/mean(res));
    Rvector = R0:-1:1;
    [deconvfm,u_phase,msk] = m1_bb_vSHARP(u_phase,msk,Rvector);
    ttt = u_phase; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
     m1_nifti_save(ttt,stp,['vsmv_R' num2str(max(Rvector))...
                                     '-' num2str(min(Rvector)) 'vx']);
    ttt = deconvfm; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
     m1_nifti_save(ttt,stp,'deconvfm');
    ttt = msk; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
     m1_nifti_save(ttt,stm,['msk_vsmv_R' num2str(max(Rvector))...
                                     '-' num2str(min(Rvector)) 'vx']); clear ttt
    toc %
end

if strcmpi(stepname,'bfe.poly3d')
    tic %
    disp(' '); disp('>>> 3D polynomial fit')
    [u_phase] = m1_shf_poly3d(u_phase,msk,ptb.bfe.poly3d.poly_order);
    ttt = u_phase; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
     m1_nifti_save(ttt,stp,'poly3d'); clear ttt
    toc %
end

%% QSM INVERSION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Add library dir to MATLAB path
if strcmpi(stepname,'qsm.MSDI')
    stepname = 'qsm.MSDI';
    addpath([ptbdir '/ptb/_medin_170706'])
    addpath([ptbdir '/ptb/_3DSRNCP'])
end

%% Unreliable phase mask
if strcmpi(stepname,'qsm.MSDI')
    wC = 1; flag_wC = false;
    nR = 0;
    if ptb.qsm.MSDI.reliab
        disp(' '); disp('>>> Compute phase reliability mask')
        if ptb.reliab_init_offset
            disp('Init offset')
             R_initOffset = m1_nifti_load('R_initOffset');
            if gadmpm_flag; ttt = R_initOffset; ttt=permute(ttt,[3 1 2]); ttt=ttt(:,end:-1:1,end:-1:1); R_initOffset = ttt; clear ttt; end
            wC = wC.*R_initOffset; flag_wC = true;
            nR = nR+1;
        end
        if ptb.reliab_unwr
            disp('Quality metric')
             R_unwr = m1_nifti_load('secdiff_phase');
            if gadmpm_flag; ttt = R_unwr; ttt=permute(ttt,[3 1 2]); ttt=ttt(:,end:-1:1,end:-1:1); R_unwr = ttt; clear ttt; end
            wC = wC.*R_unwr; flag_wC = true;
            nR = nR+1;
        end
        if ptb.reliab_bipolar
            disp('Bipolar offset')
             R_bipolarOffset = m1_nifti_load('R_bipolarOffset');
            if gadmpm_flag; ttt = R_bipolarOffset; ttt=permute(ttt,[3 1 2]); ttt=ttt(:,end:-1:1,end:-1:1); R_bipolarOffset = ttt; clear ttt; end
            wC = wC.*R_bipolarOffset; flag_wC = true;
            nR = nR+1;
        end
        if ptb.reliab_fitres
            disp('Fitting residual')
             R_fitres = m1_nifti_load('R_fitres');
            if gadmpm_flag; ttt = R_fitres; ttt=permute(ttt,[3 1 2]); ttt=ttt(:,end:-1:1,end:-1:1); R_fitres = ttt; clear ttt; end
            wC = wC.*R_fitres; flag_wC = true;
            nR = nR+1;
        end
        if flag_wC
            wC = msk.*(wC.^(1/nR));
            ttt = wC; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
             m1_nifti_save(ttt,stp,'wC');
        end
    end
     clear R_* ttt
    
    % Phase noise sigma
    if ~exist('N_std'); disp('N_std does not exist => N_std=1./magn'); N_std=msk./magn; end
end
    
%%
MSDI_Debug = false;
if MSDI_Debug
    stepname ='qsm.MSDI';
end

if strcmpi(stepname,'qsm.MSDI')
    tic %
    disp(' '); disp('>>> Inversion preamble')

    % Save dependencies
    SaveDep = false;
    if ~MSDI_Debug && SaveDep
        save ptb_preinv.mat magn u_phase wC N_std msk TE B0 B0_dir res ptb stm stp stepname MSDI_Debug '-v7.3'
        clear all; load ptb_preinv
    end
        
    % Debugging/development mode, use CTRL+RETURN to run MSDI only 
    if MSDI_Debug
        Debug_level = 2;
        disp(['Entering debug/devel mode ' num2str(Debug_level)])
        switch Debug_level
            case 1
                clear all; load ptb_preinv
            case 2
                addpath ~/matlab/QSMbox/master/ptb/_PhaseTools
                addpath ~/matlab/QSMbox/master/ptb/_Utils
                addpath ~/matlab/QSMbox/master/ptb/_spm12b_nifti
                addpath ~/matlab/QSMbox/master/ptb/_medin_170706
                addpath ~/matlab/QSMbox/master/ptb/_3DSRNCP
                if exist('ptb_preinv.mat'); load ptb_preinv; else load preinversion; end
                ptb.qsm.MSDI.lambda      = 500;
                ptb.qsm.MSDI.reliab_pct  = 10;
                ptb.qsm.MSDI.smv         = true;
                ptb.qsm.MSDI.R_smv       = [2 4 8 16];
                ptb.qsm.MSDI.e           = 1e-6;
                ptb.qsm.MSDI.cg_max_iter = 2000;
                save ptb_preinv.mat magn u_phase wC N_std msk TE B0 B0_dir res ptb stm stp stepname MSDI_Debug '-v7.3'
        end
    end
   
    % MSDI working dir
    date1 = datestr(now,'yymmdd');
    date2 = datestr(now,'HHMMSS');
    d0    = pwd; 
    wdir  = [d0 '/' cell2mat(ptb.pipename) '_' date1 '_' date2];
    mkdir(wdir)
    cd(wdir)
 
    % Constants
    TE_ref      = 20e-3;
    B0_ref      = 3;
    if length(TE)>1; delta_TE=TE(2)-TE(1); else delta_TE=TE; end
    matrix_size = size(u_phase);
    voxel_size  = res;
    
    % Calculate rad <-> ppm conversion factors
    gyro = 2*pi*42.57747892;
    ppm2rad = double(B0*gyro*delta_TE);
    rad2ppm = double(1/ppm2rad);
    save('ptb_ppm2rad.txt','ppm2rad','-ascii')
    save('ptb_rad2ppm.txt','rad2ppm','-ascii')
    ppm2rad = double(B0_ref*gyro*TE_ref);
    rad2ppm = double(1/ppm2rad);
    save('ptb_ppm2rad_ref.txt','ppm2rad','-ascii')
    save('ptb_rad2ppm_ref.txt','rad2ppm','-ascii')

    % Normalise phase (rad) to TE_ref*B0_ref
    u_phase(isnan(u_phase)) = 0;    u_phase(isinf(u_phase)) = 0;
    nfm         = u_phase*(TE_ref/delta_TE)*(B0_ref/B0);      clear u_phase
    delta_TE    = TE_ref;
    B0          = B0_ref;
    CF          = B0_ref*42.57747892*1e6;
    ttt = nfm; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
     m1_nifti_save(ttt,stp,'nfm_0'); clear ttt
    nfm_0       = nfm;
    
    % Mask
    msk(    isnan(msk))     = 0;    msk(    isinf(msk))     = 0;
    msk(msk>0)              = 1;
    Mask                    = msk;                            clear msk
    ttt = Mask; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
     m1_nifti_save(ttt,stp,'Mask'); clear ttt
    
    % Magnitude
    magn(   isnan(magn))    = 0;    magn(   isinf(magn))    = 0;
    ttt = magn; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
     m1_nifti_save(ttt,stp,'magn');
    clear magn ttt
    
    % Phase noise sigma
    m                       = Mask./N_std; m(isnan(m))=0; m(isinf(m))=0; m=m/mean(m(Mask>0)); %<<<<<<<<<<<<< NORMALISE
    ttt = m; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
     m1_nifti_save(ttt,stp,'W_0'); clear ttt
    N_std_0                 = 1./m; N_std_0(isnan(N_std_0))=0; N_std_0(isinf(N_std_0))=0; 
     clear m
    ttt = N_std_0; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
     m1_nifti_save(ttt,stp,'N_std_0'); clear ttt
    N_std                   = N_std_0;
       
    % Initialise loop    
    smv             = ptb.qsm.MSDI.smv;
    R_smv           = ptb.qsm.MSDI.R_smv;
    if R_smv==-1
        disp(' '); disp('Preliminary Poisson-MSDI implementation')
        MSDI_flag = -1;
        SMV_local = mean(res(:));
        R_smv = [SMV_local 4 8 12 16]; % R_smv = SMV_step:SMV_step:16;
    elseif length(R_smv)>1 && R_smv(1)==-1
        disp(' '); disp('MSDI v2 implementation')
        MSDI_flag = -1.1; % This flag will be reset to 0 (original MSDI implementation) after solving for the local scale
        SMV_local = 1; %mean(res(:));
        R_smv = [SMV_local R_smv(2:end)];
    else
        MSDI_flag = 0;
    end
    D               = m1_tl_medin_dipole_kernel(matrix_size,voxel_size,B0_dir);
    qsm_unmasked    = zeros(size(Mask));    
    cost_reg_sum    = 0;
    cost_fid_sum    = 0;
    C               = 0; % scale
    C2              = 0; % reliab
    SMV_radius      = R_smv(1);
    idx_optim       = 1;
    lambda_search   = true; % Original MSDI: MSDI_flag=0 & lambda_search=true throughout 
    fine_inversion  = true; % Exhaustive inversion with lambda_optim (resulting from lambda_search)
    
    while SMV_radius<=max(R_smv(:))
        tic %
        
        C=C+1; C2=C2+1; 
        if ~lambda_search && round(MSDI_flag)==-1
            C=C-1; C2=C2-1;
        end
        
        SMV_radius = R_smv(C);
            
        if smv
            if length(R_smv)>1
                disp(' '); disp(['>>> MSDI (scale #' num2str(C) ' of ' num2str(length(R_smv)) ')'])
            else
                disp(' '); disp(['>>> nMEDI with SMV sharpening (single step)'])
            end
        else
            if length(R_smv)>1
                disp(' '); disp(['>>> Iterative nMEDI (iter #' num2str(C) ' of ' num2str(length(R_smv)) ')'])
            else
                disp(' '); disp(['>>> Standard nMEDI'])
            end
        end
        
        if MSDI_flag==0              
            if SMV_radius==0
                disp(' '); disp('Local scale enabled')
                SMV_radius  = mean(res(:));
                L           = 1;
            elseif SMV_radius==-0.01
                disp(' '); disp('Local scale enabled')
                SMV_radius  = mean(res(:));
                disp('Local scale parameter sweep enabled')
                L           = ptb.qsm.MSDI.search_L; %logspace(1,4,10);
            else
                L           = 1;
            end
        elseif round(MSDI_flag)==-1
            disp(' ')
            if lambda_search
                disp('Parameter sweep enabled')
            else
                disp('Fine inversion with optimum lambda')
            end
            if C==1 && lambda_search
                L   = ptb.qsm.MSDI.search_L; %logspace(4,5,20);
                L   = L(end:-1:1);
            else
                L = L_bg;
            end
            % full_sweep = 0:   minimal grid search [discard lambdas + stop when min(diff_poisson) reached]
            % full_sweep = 0.5: do not discard lambdas
            % full_sweep = 1:   full grid search
            full_sweep = ptb.qsm.MSDI.search_full_sweep;
            if (full_sweep==0 && MSDI_flag==-1) || MSDI_flag==-1.1 
                if idx_optim>1
                    disp(['Discard L = [' num2str(round(L(1:idx_optim-1))) ']']) 
                end
                L   = L(idx_optim:end);
                L_bg = L;
            end
            if lambda_search
                disp(['Lambda grid, L = [' num2str(round(L)) ']'])
            else
                idx_optim = 1;
                L       = L(1);
                disp(['Optim lambda, L = [' num2str(round(L)) ']'])
            end
        end
        
        if ptb.qsm.MSDI.reliab==1
            if ptb.qsm.MSDI.reliab_excl1 && C2==1
                disp(' '); disp('Phase reliability adjustment disabled')
                wC           = 1;
            else
                disp(' '); disp('Phase reliability adjustment enabled')
                F            = R_smv(C)/4; % R_smv(2); <===== fixed to 4 mm normalisation 
                reliab_pct   = F*ptb.qsm.MSDI.reliab_pct; if reliab_pct>100; reliab_pct=100; end
                disp(['Exclude the ' num2str(reliab_pct) '% most unreliable measurements from the consistency cost'])
                copyfile([d0 '/wC.nii'],[wdir '/wC.nii'])
                 [wC]         = m1_nifti_load('wC');
                if gadmpm_flag; ttt = wC; ttt=permute(ttt,[3 1 2]); ttt=ttt(:,end:-1:1,end:-1:1); wC = ttt; clear ttt; end
                [wC]         = Mask.*(1-m1_pctmask(wC,Mask,reliab_pct));
                ttt = wC; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
                 m1_nifti_save(ttt,stm,['qsm_msk_wC_' num2str(C)]); clear ttt
            end
        elseif ptb.qsm.MSDI.reliab==2
            wC           = 2^-log2(R_smv(C));
        elseif ptb.qsm.MSDI.reliab==3
            wC           = 1./sqrt(R_smv(C));
        elseif ptb.qsm.MSDI.reliab==4
            reliab_pct   = 100*(1-1./(sqrt(R_smv(C))));
            copyfile([d0 '/wC.nii'],[wdir '/wC.nii'])
             [wC]         = m1_nifti_load('wC');
            if gadmpm_flag; ttt = wC; ttt=permute(ttt,[3 1 2]); ttt=ttt(:,end:-1:1,end:-1:1); wC = ttt; clear ttt; end
            [wC]         = Mask.*(1-m1_pctmask(wC,Mask,reliab_pct));
            ttt = wC; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
             m1_nifti_save(ttt,stm,['qsm_msk_wC_' num2str(C)]); clear ttt
        end
         
        downsampl = false;
        if smv && SMV_radius/mean(res(:))>=2 && downsampl 
            disp(' '); disp('Downsample data')
            downsampl_factor    = 2; %floor(SMV_radius/mean(res(:)))
            nfm                 = m1_downsampl(nfm,downsampl_factor);
            N_std_hires         = N_std;
            N_std               = m1_downsampl(N_std,downsampl_factor);
            N_std_0_hires       = N_std_0;
            N_std_0             = m1_downsampl(N_std_0,downsampl_factor);
            Mask_hires          = Mask;
            Mask                = m1_downsampl(Mask,downsampl_factor);
            Mask(Mask>=0.5)     = 1;
            Mask(Mask<0.5)      = 0;
            nfm                 = nfm.*Mask;
            N_std               = N_std.*Mask;
            N_std_0             = N_std_0.*Mask;
            if exist('wC') && length(wC(:))>1
                wC              = m1_downsampl(wC,downsampl_factor);
                wC              = wC.*Mask;
            end
            matrix_size_orig    = matrix_size;
            matrix_size         = size(nfm);
            downsampl_vec       = matrix_size_orig./matrix_size
            voxel_size_orig     = voxel_size;
            voxel_size          = voxel_size_orig.*downsampl_vec;
            SaveNifti = false;
            if SaveNifti
                ttt = nfm; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
                 m1_nifti_save(ttt,stp,'nfm_lowres');
                ttt = Mask; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
                 m1_nifti_save(ttt,stp,'Mask_lowres');
                ttt = N_std; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
                 m1_nifti_save(ttt,stp,'N_std_lowres');
                if exist('wC') && length(wC(:))>1
                    ttt = wC; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
                     m1_nifti_save(ttt,stp,'wC_lowres');
                end
                 clear ttt
            end
        end
        
        if smv
            disp(' '); disp(['SMV filtering (kernel radius: ' num2str(SMV_radius) 'mm)'])
            nfm_prev    = nfm;
            [nfm]       = Mask.*m1_SMVfilt(nfm,voxel_size,SMV_radius);
            nfm         = nfm_prev-nfm; clear nfm_prev
            ttt = nfm; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
             m1_nifti_save(ttt,stp,['SMVII' num2str(SMV_radius) '_nfm']); clear ttt
            if ptb.qsm.MSDI.smv_N_std_update
                [N_std]     = Mask.*m1_SMVfilt(N_std,voxel_size,SMV_radius);
                [Mask_SMV]  = Mask.*m1_SMVfilt(Mask,voxel_size,SMV_radius);
                m           = Mask_SMV./N_std; m(isnan(m))=0; m(isinf(m))=0; m = m./mean(m(Mask>0)); %<<<<<<<<<<<<< NORMALISE  
                 clear Mask_SMV
                N_std       = Mask./m; N_std(isnan(N_std))=0; N_std(isinf(N_std))=0;
                 clear m
                N_std       = sqrt(N_std_0.^2+N_std.^2);
            else
                N_std       = N_std_0; %<<<<<<<<<<<<< RESET N_STD_0
            end
        else
            SMV_radius = 0;
        end
        
        if smv
            disp(' '); disp('Initialise dipole inversion with SMV deconvolution')
        else
            disp(' '); disp('Initialise nMEDI')
        end
        RDF = nfm; clear nfm
        lambda                             = round(ptb.qsm.MSDI.lambda);
        e                                  = ptb.qsm.MSDI.e;
        cg_max_iter                        = ptb.qsm.MSDI.cg_max_iter;
        cg_tol                             = ptb.qsm.MSDI.cg_tol;
        max_iter                           = ptb.qsm.MSDI.max_iter;
        tol_norm_ratio                     = ptb.qsm.MSDI.tol_norm_ratio;
        data_weighting_mode                = ptb.qsm.MSDI.data_weighting_mode;
        if ptb.qsm.MSDI.gradient_weighting_mode==3
            if C==1
                gradient_weighting_mode    = 1;
            else
                gradient_weighting_mode    = 0;
            end
        else
                gradient_weighting_mode    = ptb.qsm.MSDI.gradient_weighting_mode;
        end    
        dynedge_mode                       = ptb.qsm.MSDI.dynedge_mode;
        merit                              = ptb.qsm.MSDI.merit;
        merit_f                            = ptb.qsm.MSDI.merit_f;
        merit_p                            = ptb.qsm.MSDI.merit_p;

        for sweep = 1:length(L)    
            if ptb.qsm.MSDI.fidelity_boost
                if MSDI_flag==0
                    if C==1
                        if length(L)>1 
                            ptb.qsm.MSDI.fidelity_boost_factor = L(sweep); 
                        end
                        lambda    = round(ptb.qsm.MSDI.fidelity_boost_factor*ptb.qsm.MSDI.lambda);
                        max_iter  = ptb.qsm.MSDI.fidelity_boost_maxiter;
                        C2        = 0;
                    else
                        lambda    = round(ptb.qsm.MSDI.lambda);
                        max_iter  = ptb.qsm.MSDI.max_iter;
                    end
                end
            end
            if round(MSDI_flag)==-1
                if lambda_search
                    lambda      = round(L(sweep));
                    cg_max_iter = ptb.qsm.MSDI.search_cg_max_iter;
                    cg_tol      = 0;
                    max_iter    = 1;
                else
                    lambda      = round(L);
                    cg_max_iter = 1000;
                    cg_tol      = ptb.qsm.MSDI.fine_cg_tol;
                    max_iter    = 1;
                end
            end
            if SMV_radius>6
                if tol_norm_ratio < 0.1
                    tol_norm_ratio = 0.1;
                end
            end

            % Edge masking type
            % [0] Static I:     magnitude (all scales)
            % [1] Dynamic I:    magnitude (scale #1) + QSM (subsequently)
            % [2] Dynamic II:   scale-dependent SMV filtered field map [if smv]
            % [3] Static II:    input field map, i.e. nfm_0 (all scales)
            if sweep==1
                if (dynedge_mode==0 || dynedge_mode==1) && C==1
                    if exist('magn.nii.gz')
                        copyfile('magn.nii.gz','edgeima.nii.gz');
                    elseif exist('magn.nii')
                        copyfile('magn.nii','edgeima.nii');
                    end
                elseif dynedge_mode==2
                    if smv
                        R_smvkernel = SMV_radius; % dynamic
                        if exist(['SMVII' num2str(R_smvkernel) '_nfm.nii.gz'])
                            copyfile(['SMVII' num2str(R_smvkernel) '_nfm.nii.gz'],'edgeima.nii.gz')
                        elseif exist(['SMVII' num2str(R_smvkernel) '_nfm.nii'])
                            copyfile(['SMVII' num2str(R_smvkernel) '_nfm.nii'],'edgeima.nii')
                        end
                    else
                        error('ERROR: You set dynedge_mode=2, but smv=false (smv=true required). Edit your ptbs file and start again.')
                    end
                elseif dynedge_mode==3 && C==1
                    disp('Static edge masking across scales. Edge mask inferred from input field map, i.e. nfm_0')
                    if exist('nfm_0.nii.gz')
                        copyfile('nfm_0.nii.gz','edgeima.nii.gz');
                    elseif exist('nfm_0.nii')
                        copyfile('nfm_0.nii','edgeima.nii');
                    end
                elseif dynedge_mode==4
                    if smv && C==1
                        R_smvkernel = R_smv(1); % static
                        if exist(['SMVII' num2str(R_smvkernel) '_nfm.nii.gz'])
                            copyfile(['SMVII' num2str(R_smvkernel) '_nfm.nii.gz'],'edgeima.nii.gz')
                        elseif exist(['SMVII' num2str(R_smvkernel) '_nfm.nii'])
                            copyfile(['SMVII' num2str(R_smvkernel) '_nfm.nii'],'edgeima.nii')
                        end
                    else
                        error('ERROR: You set dynedge_mode=4, but smv=false (smv=true required). Edit your ptbs file and start again.')
                    end
                end
            end
                
            save RDF.mat RDF N_std Mask matrix_size voxel_size delta_TE CF B0_dir '-v7.3'
            save genDI.mat stp ptb e cg_max_iter cg_tol max_iter tol_norm_ratio...
                 data_weighting_mode gradient_weighting_mode wC SMV_radius merit_f merit_p...
                 MSDI_flag downsampl_factor gadmpm_flag '-v7.3'
            
            [qsm,cost_reg,cost_fid,N_CG_iter,CG_residual] = medin_170706b('lambda',lambda,'smv',smv,'merit',merit); %,'pad',[100 100 100]);
            
            N_CG_iter   = N_CG_iter(end);
            CG_residual = CG_residual(end);
            
            if round(MSDI_flag)==-1
                cost_reg_all(C,sweep)   = cost_reg;
                cost_fid_all(C,sweep)   = cost_fid;
                N_iter(C,sweep)         = N_CG_iter;
                CG_rsd(sweep)           = CG_residual;
            else
                cost_reg_all(C)         = cost_reg;
                cost_fid_all(C)         = cost_fid;
            end
               
            CopywGsum = false;
            if CopywGsum
                if (gradient_weighting_mode && C==1) || (gradient_weighting_mode && dynedge_mode)
                    if exist('wGsum.nii.gz')
                        copyfile('wGsum.nii.gz',['qsm_wGsum_' num2str(C) '.nii.gz'])
                    elseif exist('wGsum.nii')
                        copyfile('wGsum.nii',['qsm_wGsum_' num2str(C) '.nii'])
                    end
                end
            end
            
            if downsampl_factor>1
                % Upsample QSM
                upsampl_vec = downsampl_vec
                qsm = m1_upsampl(qsm,upsampl_vec);  
                Mask_lowres = Mask;
                Mask        = Mask_hires;
            end
               
            qsm_unmasked_temp = qsm_unmasked+qsm; clear qsm

            if smv
                if length(L)>1 || (MSDI_flag==-1 && lambda_search)
                    qsm_outroot = ['qsm_INTEGRAL_' num2str(C) '_MSDI_l' num2str(lambda) '_sweep' num2str(sweep)];
                else
                    qsm_outroot = ['qsm_INTEGRAL_' num2str(C) '_MSDI_l' num2str(lambda)];
                end
            else
                if length(L)>1 || MSDI_flag==-1
                    qsm_outroot = ['qsm_INTEGRAL_' num2str(C) '_inMEDI_l' num2str(lambda) '_sweep' num2str(sweep)];
                else
                    qsm_outroot = ['qsm_INTEGRAL_' num2str(C) '_inMEDI_l' num2str(lambda)];
                end
            end
            
            disp(' ')
       
            if length(L)>1
                ttt = qsm_unmasked_temp; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
                 m1_nifti_save(ttt,stp,qsm_outroot); clear ttt
            else
                ttt = Mask.*qsm_unmasked_temp; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
                 m1_nifti_save(ttt,stp,qsm_outroot); clear ttt
            end
                        
            if length(L)>1 || round(MSDI_flag)==-1
                poisson_flag = 1;
                if poisson_flag==1
                    [diff_poisson,left] = m1_phase_poisson(Mask.*nfm_0*rad2ppm,Mask.*qsm_unmasked_temp,res);
                elseif poisson_flag==2
                    [diff_poisson] = m1_phase_poisson(Mask.*RDF*rad2ppm,Mask.*(qsm_unmasked_temp-qsm_unmasked),res);
                elseif poisson_flag==3
                    if C==1 && sweep==1
                        lfm = RDF;
                    end
                    [diff_poisson] = m1_phase_poisson(Mask.*lfm*rad2ppm,Mask.*(qsm_unmasked_temp-qsm_unmasked),res);    
                end
                diff_poisson_norm(sweep)  = norm(diff_poisson(Mask>0),2) / numel(diff_poisson(Mask>0));
                diff_poisson_norm4(sweep) = norm(diff_poisson(Mask>0)./left(Mask>0),2) / numel(diff_poisson(Mask>0));
                diff_poisson_norm3(sweep) = diff_poisson_norm(sweep)*CG_residual;
                 clear diff_poisson left
                disp(' ')
                disp(['Lambda               = [' num2str(round(L(1:sweep))) ']'])
                disp(['Poisson residual #1  = [' num2str(diff_poisson_norm) ']'])
                disp(['Poisson residual #2  = [' num2str(diff_poisson_norm4) ']'])
                disp(['CG residual          = [' num2str(CG_rsd) ']'])
                disp(['Composite residual   = [' num2str(diff_poisson_norm3) ']'])
                % Stopping criterion
%                 diff_poisson_norm2(sweep) = diff_poisson_norm3(sweep);
%                 diff_poisson_norm2(sweep) = diff_poisson_norm(sweep);
%                 diff_poisson_norm2(sweep) = diff_poisson_norm4(sweep);
                diff_poisson_norm2(sweep) = CG_rsd(sweep);
            end
            
            disp(' ')
            
            if round(MSDI_flag)==-1 && sweep>1 && diff_poisson_norm2(sweep)-diff_poisson_norm2(sweep-1)>0 && full_sweep<1
                break
            end
            
            if downsampl_factor>1
                Mask = Mask_lowres;
            end            
        end % sweep
        
        if downsampl_factor>1
            Mask = Mask_hires;
        end    
        
        if length(L)>1 || round(MSDI_flag)==-1
            if lambda_search
                disp(['Optimal reconstruction (scale #' num2str(C) '):'])
            else
                disp(['Summary of progress (scale #' num2str(C) '):'])
            end
            
            solution_cost = diff_poisson_norm2(:);
            min_solution_cost(C) = min(solution_cost);
            disp(['Poisson cost = [' num2str(min_solution_cost) ']'])
            
            if lambda_search
                idx_optim = find(solution_cost==min_solution_cost(C));           
            end
            
            if lambda_search && round(MSDI_flag)==-1
                save 'lambda_search_L_sweep.txt' 'L' '-ascii'                
                save 'lambda_search_residual_poisson_1.txt' 'diff_poisson_norm' '-ascii'
                save 'lambda_search_residual_poisson_2.txt' 'diff_poisson_norm4' '-ascii'
                save 'lambda_search_residual_CG.txt' 'CG_rsd' '-ascii'
                save 'lambda_search_residual_composite.txt' 'diff_poisson_norm3' '-ascii'
            end
             clear diff_poisson_norm diff_poisson_norm2 CG_rsd
             
            if lambda_search
                if round(MSDI_flag)==-1
                    idx_optim = idx_optim(1);
                else
                    idx_optim = idx_optim(end);
                end
                
                if ptb.qsm.MSDI.fidelity_boost
                    ptb.qsm.MSDI.fidelity_boost_factor = L(idx_optim); 
                    lambda = round(ptb.qsm.MSDI.fidelity_boost_factor*ptb.qsm.MSDI.lambda);
                else
                    lambda = round(L(idx_optim));
                end
                lambda_optim(C) = lambda;
            end
            
            disp(['Lambda optim = [' num2str(lambda_optim) ']']); disp(' ')
            
            if lambda_search
                if smv
                    qsm_optimroot = ['qsm_INTEGRAL_' num2str(C) '_MSDI_l' num2str(lambda) '_sweep' num2str(idx_optim)];
                else
                    qsm_optimroot = ['qsm_INTEGRAL_' num2str(C) '_inMEDI_l' num2str(lambda) '_sweep' num2str(idx_optim)];
                end
                 clear qsm_unmasked_temp
                qsm_unmasked_old = qsm_unmasked;
                 qsm_unmasked = m1_nifti_load(qsm_optimroot);
                if gadmpm_flag; ttt = qsm_unmasked; ttt=permute(ttt,[3 1 2]); ttt=ttt(:,end:-1:1,end:-1:1); qsm_unmasked = ttt; clear ttt; end
                 delete('qsm_INTEGRAL*sweep*.nii*')
                disp(' ')
            else
                qsm_unmasked_old = qsm_unmasked;
                qsm_unmasked = qsm_unmasked_temp;
                 clear qsm_unmasked_temp
            end
        else
            qsm_unmasked_old = qsm_unmasked;
            qsm_unmasked = qsm_unmasked_temp;
             clear qsm_unmasked_temp
        end
        
        scale_poisson_monitor = true;
        if scale_poisson_monitor
            [diff_poisson] = m1_phase_poisson(Mask.*nfm_0*rad2ppm,Mask.*(qsm_unmasked-qsm_unmasked_old),res);
%             diff_poisson_norm_scale(C) = norm(diff_poisson(Mask>0)./left(Mask>0),2) / numel(diff_poisson(Mask>0));
            diff_poisson_norm_scale(C) = norm(diff_poisson(Mask>0)) / numel(diff_poisson(Mask>0));
            clear diff_poisson left
            disp(['Poisson residual (scale differential) = [' num2str(diff_poisson_norm_scale) ']']); disp(' ')
%             if C>1 && diff_poisson_norm_scale(C)-diff_poisson_norm_scale(C-1)>0
%                 break
%             end
        end
         clear RDF

        if length(L)>1 || round(MSDI_flag)==-1
            if smv
                qsm_outroot = ['qsm_INTEGRAL_' num2str(C) '_MSDI_l' num2str(lambda)];
            else
                qsm_outroot = ['qsm_INTEGRAL_' num2str(C) '_inMEDI_l' num2str(lambda)];
            end
            ttt = Mask.*qsm_unmasked; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
             m1_nifti_save(ttt,stp,qsm_outroot); clear ttt
        end
             
        if C>1
            if smv
                qsm_outroot = ['qsmII_' num2str(C) '_MSDI_l' num2str(lambda)];
            else
                qsm_outroot = ['qsmII_' num2str(C) '_inMEDI_l' num2str(lambda)];
            end
            qsm = qsm_unmasked-qsm_unmasked_old;
            ttt = Mask.*qsm; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
             m1_nifti_save(ttt,stp,qsm_outroot);
             clear qsm qsm_unmasked_old ttt
        end
        
        % Edge masking type
        % [1] Dynamic I:    magnitude (scale #1) + QSM (subsequently)
        if dynedge_mode==1
            if exist([qsm_outroot '.nii.gz'])
                copyfile([qsm_outroot '.nii.gz'],'edgeima.nii.gz');
            elseif exist([qsm_outroot '.nii'])
                copyfile([qsm_outroot '.nii'],'edgeima.nii');
            end
        end
        
        if length(L)==1 || C>1
            if smv
                qsm_outroot = ['qsm_UNMASKED_MSDI'];
            else
                qsm_outroot = ['qsm_UNMASKED_inMEDI'];
            end
            ttt = qsm_unmasked; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
             m1_nifti_save(ttt,stp,qsm_outroot); clear ttt
        end
         
        disp(' '); disp('Compute remnant field')
        qsm_ttt     = qsm_unmasked*(2*pi*delta_TE*CF)*1e-6;
        if downsampl_factor>1
            N_std       = N_std_hires;
            N_std_0     = N_std_0_hires;
        end
        m           = Mask./N_std; m(isnan(m))=0; m(isinf(m))=0;
        S1          = Mask.*m.*exp(1i*real(ifftn(D.*fftn(qsm_ttt)))); clear qsm_ttt m
        nfm1        = angle(S1); clear S1
        [nfm1]      = m1_shf_3DSRNCP(nfm1,Mask,ptb.opt.OS);
        nfm         = Mask.*(nfm_0-nfm1); clear nfm1
        ttt = nfm; if gadmpm_flag; ttt=ttt(:,end:-1:1,end:-1:1); ttt=permute(ttt,[2 3 1]); end
         m1_nifti_save(ttt,stp,['qsm_nfmdiff_' num2str(C) '_l' num2str(lambda)]); clear ttt
        N_std       = N_std_0; % this assumes N_std_nfm1=0
         
        disp(' '); toc %
        
        if SMV_radius==max(R_smv(:)) && ( ( MSDI_flag==-1 && ( (lambda_search && ~fine_inversion) || (~lambda_search && fine_inversion) ) )  || MSDI_flag==0 )
            break
        end

        if lambda_search && fine_inversion && round(MSDI_flag)==-1
            lambda_search = false;
        else
            lambda_search = true;
            if MSDI_flag==-1.1
                MSDI_flag = 0;
            end
        end        
    end % for SMV_radius, C (scale counter)
    
    
    cost_reg_all = double(cost_reg_all); cost_fid_all = double(cost_fid_all);
    save(['cost_reg_l' num2str(lambda) '.txt'],'cost_reg_all','-ascii');
    save(['cost_fid_l' num2str(lambda) '.txt'],'cost_fid_all','-ascii');
    msdiout = 1;
    save('msdiout.txt','msdiout','-ascii');
    cd(d0); save('msdiout.txt','msdiout','-ascii'); cd(wdir)
    copyfile([d0 '/ptbs_' cell2mat(ptb.pipename) '.m'],[wdir '/ptbs_' cell2mat(ptb.pipename) '.m'])

    % Housekeeping
    if ~MSDI_Debug
        delete('RDF.mat')
        delete('genDI.mat'); 
%         delete('x.nii*')
    end

    cd(d0)

end % MSDI

end % for stepname=

%% REPORTING %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(' '); disp('>>> Save ptb structure')
save('ptb.mat','ptb')
disp(' ')
disp(['>>> Finished pipeline: ' cell2mat(ptb.pipename)]);
disp(' ')
for x=1:length(ptb.pipeline); disp([cell2mat(ptb.pipeline(x))]); end
disp(' ')
disp(['Total processing time: ' num2str(round(etime(clock,ptb.clock)/60)) ' minutes.'])
disp(' ');
