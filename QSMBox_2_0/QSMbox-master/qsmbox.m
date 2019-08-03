function qsmbox(pwdiswdir,importdata,scantype,pipeID,pipeoper)
% QSMbox v2.0 (release candidate)
%  MATLAB software for quantitative susceptibility mapping (QSM) reconstruction
% 
% SINTAX  
%  qsmbox
%   Start dialogue
% 
%  qsmbox(pwdiswdir)
%   Is pwd a suitable working directory? ['y'/'n']
% 
%  qsmbox(pwdiswdir,importdata)
%   Would you like to import data from another directory? ['y'/'n']
% 
%  qsmbox(pwdiswdir,importdata,scantype)
%   Select datatype: 
%    [1] Coil-combined single-GRE
%    [2] Coil-combined multi-GRE
%    [3] Coil-uncombined single-GRE
%    [4] Pre-processed data
% 
%  qsmbox(pwdiswdir,importdata,scantype,pipeID)
%   Select preset pipeline
% 
%  qsmbox(pwdiswdir,importdata,scantype,pipeID,pipeoper)
%   [1] Edit settings, or [2] Run?
% 
% EXAMPLE
%  qsmbox('y','n',1,1,1)
% 
% DOWNLOAD
%  https://gitlab.com/acostaj/QSMbox.git
% 
% If you use this software for research, please cite: 
%  Acosta-Cabronero et al. Neuroimage 2018


disp(' ')
disp('QSMbox v2.0 (release candidate) is experimental software developed by Julio Acosta-Cabronero')

% Verbose?
boxverbose = false;

% Establish working directory
d0=pwd;
disp(' ')
disp(['You are here: ' d0])
disp(' ')
if ~exist('pwdiswdir','var')
    pwdiswdir = input('Is this your working directory? [y/n]: ','s');
end
switch pwdiswdir
    case 'y'
    case 'n'
        newwdir = uigetdir('','Select working directory');
        cd(newwdir)
    otherwise
        disp(' ')
        error('Wrong selection. Start again.')
end

% Import data (*_orig.nii*, ptb_*.txt and dicom/) from a remote directory
if ~exist('importdata','var')
    disp(' ')
    importdata = input('Would you like to import data from another directory? [y/n]: ','s');
    disp(' ')
end
switch importdata
    case 'y'   
        datadir = uigetdir('','Select data directory');
        
        % Copy NIFTIs
        ldata = '*_orig.nii*';
        nifti_flag = false;
        if isempty(dir([datadir '/' ldata]))
             disp('Warning: NIFTI data, magn_orig.nii(.gz) and phase_orig.nii(.gz), NOT found')
             disp(' ')
        else
            if isempty(dir(ldata))                
                nifti_flag = true;
            else
                confirm = input('NIFTI files will be overwritten, would you like to continue? [y/n]: ','s'); 
                disp(' ')
                switch confirm
                    case 'y'
                        nifti_flag = true;
                end
            end
        end
        if nifti_flag
            disp('Copy NIFTI files:')
            disp(' ')
            copyfile([datadir '/' ldata],pwd,'f');
            ls('*_orig.nii*')
        end
        
        % Copy parameter files
        ldata = 'ptb_*.txt';
        param_flag = false;
        if isempty(dir([datadir '/' ldata]))
             disp('Warning: Parameter files, ptb_*.txt, NOT found')
             disp(' ')
        else
            if isempty(dir(ldata))
                param_flag = true;
            else
                confirm = input('Parameter files (i.e. ptb_*.txt) will be overwritten, would you like to continue? [y/n]: ','s'); 
                disp(' ')
                switch confirm
                    case 'y'
                        param_flag = true;
                end
            end
        end
        if param_flag
            disp('Copy parameter files:')
            disp(' ')
            copyfile([datadir '/' ldata],pwd,'f');
            ls('ptb_*.txt')
        end
        
        % Copy 'dicom' dir
        dicom_flag = false;
        if exist([datadir '/dicom'],'dir')
            if exist('dicom','dir')
                confirm = input('Dicom directory will be overwritten, would you like to continue? [y/n]: ','s');
                disp(' ')
                switch confirm
                    case 'y'
                        rmdir('dicom','s')
                        dicom_flag = true;
                end
            else
                dicom_flag = true;
            end
        else
            disp('Warning: ''dicom'' directory NOT found.')
            if param_flag
                disp(' ')
                disp('Though if all parameter files are available, DICOM header information might not be necessary.')
            end
            disp(' ')
        end
        if dicom_flag
            disp('Copy ''dicom'' directory:')
            disp(' ')
            copyfile([datadir '/dicom'],'dicom');
            ls('dicom/*')
        end
        
    case 'n'
    otherwise
        error('Wrong selection. Start again.')
end

% Select datatype
boxdir = fileparts(which(mfilename));
masterdir=[boxdir '/master'];
addpath(masterdir)
presetsdir=[boxdir '/presets'];
disp('DATATYPE')
disp(' ')
disp('->[1] Coil-combined single-GRE   |  full pipeline | ''cse'' |')
disp('->[2] Coil-combined multi-GRE    |  full pipeline | ''cme'' |')
disp('->[3] Coil-uncombined single-GRE |  full pipeline | ''use'' |')
disp('->[4] Add-ons                    |  partial proc. | ''ppd'' |')
disp(' ')
if ~exist('scantype','var')
    scantype = input('Select: ');
else
    disp(['Select: ' num2str(scantype)])
end
    
% Select pipeline
disp(' ')
disp('PRESET OPTIONS')
disp(' ')
disp('N.B. The MSDI algorithm is described in Acosta-Cabronero et al. Neuroimage 2018.')
disp('     MSDI v2 is the latest implementation with preservation of phase noise texture')
disp('     and low susceptibility contrast features.')

disp(' ')
switch scantype
    case 1
        scantyperoot = 'cse_';
        Npresets = print_pipelines(presetsdir,scantyperoot,boxverbose);
    case 2
        scantyperoot = 'cme_';
        Npresets = print_pipelines(presetsdir,scantyperoot,boxverbose);
    case 3
        scantyperoot = 'use_';
        Npresets = print_pipelines(presetsdir,scantyperoot,boxverbose);
    case 4
        scantyperoot = 'ppd_';
        Npresets = print_pipelines(presetsdir,scantyperoot,boxverbose);
    otherwise
        error('Wrong selection. Start again.')
end
if Npresets==1
    pipeID = 1;
else
    if ~exist('pipeID','var')
        pipeID = input('Select: ');
    else
        disp(['Select: ' num2str(pipeID)])
    end
end
    
% Edit or run pipeline?
disp(' ')
disp('SELECT OPERATION')
disp(' ')
if ~exist('pipeoper','var')
    pipeoper = input('->[1] Edit settings, or ->[2] Run? ');
else
    disp(['->[1] Edit settings, or ->[2] Run? ' num2str(pipeoper)])
end
switch pipeoper
    case 1
        pipeline_manager(presetsdir,scantyperoot,pipeID,'edit')
    case 2
        pipeline_manager(presetsdir,scantyperoot,pipeID,'run')
    otherwise
        disp(' ')
        error('Wrong selection. Start again.')
end

%% FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function Npresets = print_pipelines(presetsdir,scantyperoot,boxverbose)
    presets  = lsdirs(presetsdir,scantyperoot);
    Npresets = length(presets);
    for x=1:Npresets
        presetroot = cell2mat(presets(x));
        addpath([presetsdir '/' presetroot])
        flag_arch = false;
        if exist([pwd '/ptbs_' presetroot '.m'],'file')
            flag_arch = true;
            arch_uid  = round(now*1e6);
            movefile([pwd '/ptbs_' presetroot '.m'],[pwd '/ptbs_' presetroot '_' num2str(arch_uid) '.m']); 
        end
        eval(['ptbs_' presetroot '(''dummy'')'])
        disp([' ->[' num2str(x) '] ' presetroot])
        if boxverbose
            if flag_arch
                disp(['./ptbs_' presetroot '.m -> ./ptbs_' presetroot '_' num2str(arch_uid) '.m'])
            end 
        end
        disp(' ')
    end
end

function pipeline_manager(presetsdir,scantyperoot,pipeID,pipeopername)
    presets = lsdirs(presetsdir,scantyperoot);
    presetroot = cell2mat(presets(pipeID));
    settings_file = ['ptbs_' presetroot '.m'];
    copyfile([presetsdir '/' presetroot '/' settings_file],['./' settings_file])
    if strcmpi(pipeopername,'edit')
        disp(' ')
        edit(settings_file)
    elseif strcmpi(pipeopername,'run')
        disp(' ')
        eval(['ptbs_' presetroot '(''run'')'])
    end
end

function [presets] = lsdirs(presetsdir,scantyperoot)
    if ~exist('presetsdir','var')
        presetsdir = pwd;
    end
    if ~exist('scantyperoot','var')
        d = dir(presetsdir);
    else
        d = dir([presetsdir '/' scantyperoot '*']);
    end
    isub = [d(:).isdir]; %# returns logical vector
    presets = {d(isub).name}';
    presets(ismember(presets,{'.', '..'})) = [];
end

end % qsmbox
