function m1_spmsegment(fname,settings)
% SPM12 segmentation and bias correction
%
% m1_spmsegment(filename,settings)
%
% INPUTS
%  NIFTI(.gz) file in pwd
%  SETTINGS
%   + '3T_default'
%   + '7T_default'
%
% NOTE
%  Best to run each image in a dedicated directory
%
% Created by Julio Acosta-Cabronero

% Default SPM Segment settings
if nargin < 2
    settings = '3T_default';
end

% Gunzip if necessary
[froot,ext] = m1_fname_split(fname);
if strcmpi(ext,'gz')
    gunzip(fname);
    delete(fname)
    fname = froot;
end

% Locate job file
jspmdir = fileparts(which(mfilename));
jobfile = {[jspmdir '/m1_spmsegment_job.m']};

% Vars needed by job file
save('fname.mat','fname');
save('settings.mat','settings');

% Run job file
spm('defaults', 'PET');
spm_jobman('initcfg');
spm_jobman('run',jobfile);

% % Gzip
% m1_gzip;
