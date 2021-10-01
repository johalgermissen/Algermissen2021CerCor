function [dirs] = set_dirs(rootdir)

% Set directories for EEG pre-processing and analyses based on root 
% directory. Add paths for Fieldtrip.
% 
% INPUT:
% rootdir           =  root directory relative to which paths are initialized
% (default '/project/3017042.02')
% 
% OUTPUT:
% dirs              = structure with directories
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

%% Set directories:

% Check if root directory provided or not:
if nargin < 1
    dirs.root = '/project/3017042.02'; % root directory--each user needs to adapt this to their own folder structure
    fprintf('No root directory specified, assume %s',dirs.root)
else
    dirs.root = rootdir;
end

% ----------------------------------------------------------------------- %
fprintf('Set directories \n')

% Central directory with EEG data:
dirs.data           = fullfile(dirs.root,'Log/EEG');

% Directory with helper files for any analyses:
dirs.helpers        = fullfile(dirs.root,'Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Helpers');

% Polhemus raw data:
dirs.polhemus       = fullfile(dirs.data, 'EEG_electrode_pos');

% Behavioral raw data (.mat files):
dirs.behavior       = fullfile(dirs.root,'Behavior/Data_beh_mat');

% EEG raw data:
dirs.raw            = fullfile(dirs.data, 'MRcorrected');

% Directory where all pre-processed files go:
dirs.results        = fullfile(dirs.data, 'CueLockedResults');
if ~exist(dirs.results,'dir'); mkdir(dirs.results); end

% Preprocessing directories:
dirs.preICA         = fullfile(dirs.results, 'preICA');
if ~exist(dirs.preICA,'dir'); mkdir(dirs.preICA); end

dirs.ICA            = fullfile(dirs.results, 'ICA');
if ~exist(dirs.ICA,'dir'); mkdir(dirs.ICA); end

dirs.ICAplots       = fullfile(dirs.results, 'ICAplots');
if ~exist(dirs.ICAplots,'dir'); mkdir(dirs.ICAplots); end

% Code as logical index with either included (0) or excluded (1)
dirs.rejTrials      = fullfile(dirs.results, 'rejTrials');
if ~exist(dirs.rejTrials,'dir'); mkdir(dirs.rejTrials); end

dirs.finalPP        = fullfile(dirs.results, 'finalPP');
if ~exist(dirs.finalPP,'dir'); mkdir(dirs.finalPP); end

% TF directories set adaptively based on input settings

% Analysis directories:

% dirs.plot           = fullfile(dirs.task,'plots');
% if ~exist(dirs.plot,'dir'); mkdir(dirs.plot); end

%% Set paths to Fieldtrip and helpers directory:
% (for every script; adapt input and targetdir per script):

% ----------------------------------------------------------------------- %
fprintf('Add path to Fieldtrip\n')
addpath /home/common/matlab/fieldtrip; % add your own fieldtrip directory
ft_defaults;

fprintf('Add path to helper files directory\n')
addpath(dirs.helpers);

% END