function EEGfMRIPav_CueLocked_5_TR(rootdir,iSub)

% EEGfMRIPav_CueLocked_5_TR(rootdir,iSub)
% 
% Interactive script.
% Visualize trials with EEGlab, identify trials to reject via visual
% inspection.
% 
% INPUTS:
% rootdir           = string, root directory of project
% iSub              = scalar integer, subject number of data to inspect
%
% OUTPUTS:
% saves vector with logical indices of which trials to reject under
% dirs.rejTrials.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% By J.C.Swart, 2016/ J. Algermissen, 2018-2021.
% -------------------------------------------------------------------------
% NOTE: RUNS BEST IN MATLAB 2013B!!! NEWER VERSION MIGHT NOT SUPPORT THIS
% EEGLAB CODE!!!

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Preprocessing/

%% Complete input settings:

% clear all; close all; clc
% dbstop if error

if ~exist('rootdir','var')
    rootdir = preprocessing_set_rootdir(); % '/project/3017042.02';
    fprintf('rootdir unspecified, assume %s\n',rootdir)
end

if ~exist('iSub','var')
    iSub = 1; 
    fprintf('iSub unspecified, assume %d\n',iSub)
end

%% Set directories:

% Add helper files to path:
fprintf('Add helper files to path\n');
addpath(fullfile(rootdir,'Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Helpers'));

% Set directories and parameters:
dirs    = set_dirs(rootdir);
par     = set_par();

%% Add EEGlab:

fprintf('Add EEGlab to path\n');
dirs.eeglab         = fullfile('~/matlabScripts/eeglab12_0_1_0b');
addpath(dirs.eeglab);

%% Detect input files:

% Detect all available input folders (all subjects):
tmp         = dir(fullfile(dirs.postICA,'*.mat'));
fileList    = {tmp.name};
nInput      = length(fileList);
fprintf('Found data files from %d subjects\n',nInput);

% Extract subject numbers from folders available:
subNames    = str2double(extractBetween(fileList,8,10));

% Extract indices of selected subjects:
subIndex    = find(ismember(subNames,iSub));

%% Create output file:

% Create output file name:
outputFile = fullfile(dirs.rejTrials,sprintf('MRIpav_%03d_rejTrials.mat',iSub));

if exist(outputFile,'file')
    fprintf('Subject %03d: rejectedtrials file for already exists\n',iSub);
    continue;
end

%% Create input file:

inputFile       = fullfile(dirs.postICA, fileList{subIndex}); %#ok<*FNDSB> % go to subject directory
fprintf('EEG input file is %s\n',inputFile);

%% Load data:

fprintf('Subject %03d: Loading data ... (can take a few seconds)\n',iSub);
load(inputFile) % 
fprintf('Subject %03d: Loading completed\n',iSub);

%% Initialize EEGLAB object:

EEG                 = [];
EEG.setname         = ''; % give name
EEG.filename        = '';
EEG.filepath        = '';
EEG.subject         = '';
EEG.group           = '';
EEG.condition       = '';
EEG.session         = [];
EEG.comments        = 'Original file: '; % 'Original file: ...'
EEG.nbchan          = []; % #
EEG.trials          = []; % #
EEG.pnts            = []; % # timepoints in trial
EEG.srate           = []; % hz
EEG.xmin            = [];
EEG.xmax            = [];
EEG.times           = []; % 1xpnts array with timings
EEG.data            = []; % chan x time x trial data
EEG.icaact          = [];
EEG.icawinv         = [];
EEG.icasphere       = [];
EEG.icaweights      = [];
EEG.icachansind     = [];
EEG.chanlocs        = []; 
EEG.urchanlocs      = [];
EEG.chaninfo        = []; 
EEG.ref             = ''; 
EEG.event           = [];
EEG.urevent         = [];
EEG.eventdescription= {};
EEG.epoch           = [];
EEG.epochdescription= {};
EEG.reject          = [];
EEG.stats           = [];
EEG.specdata        = [];
EEG.specicaact      = [];
EEG.splinefile      = '';
EEG.icasplinefile   = '';
EEG.dipfit          = [];
EEG.history         = '';
EEG.saved           = 'no';
EEG.etc             = [];

% Prepare data:
EEG.nbchan      = size(data.trial{1},1);
EEG.trials      = numel(data.trial); % #
EEG.pnts        = size(data.trial{1},2); % # timepoints in trial
EEG.srate       = data.fsample; %hz
EEG.xmin        = data.time{1}(1); 
EEG.xmax        = data.time{1}(end);
EEG.times       = data.time{1}; % 1xpnts array with timings

% Transfer data from Fieldtrip format to EEGlab format:
for iTrial = 1:EEG.trials
    EEG.data(:,:,iTrial) = data.trial{iTrial}(EEG.eloc_file,:); %
end

% Plot in EEGLAB for trial rejection:
eeglab redraw 
pop_eegplot(EEG, 1, 1, 0); % 'eloc_file', [1:31 33:64 67]);

%% Retrieve vector with rejected trials.

rejectedtrials = []; % delete
rejectedtrials = EEG.reject.rejmanual; % retrieve from EEGlab
fprintf('Subject %03d: rejected %d trials: %s\n',iSub,sum(rejectedtrials),strjoin(string(find(rejectedtrials)),', '));

%% Save:

save(outputFile,'rejectedtrials');
fprintf('Subject %03d: Saved\n',iSub)

fprintf('Done :-)\n');

end % end of function.
