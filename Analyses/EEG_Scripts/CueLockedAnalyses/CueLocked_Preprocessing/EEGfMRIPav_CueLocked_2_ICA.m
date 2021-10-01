function EEGfMRIPav_CueLocked_2_ICA(rootdir,subVec)

% EEGfMRIPav_CueLocked_2_ICA(rootdir,subVec)
% 
% Perform ICA on data, save components.
% 
% INPUTS:
% rootdir           = string, root directory of project
% subVec            = vector of integers, numbers of subjects to process
%
% OUTPUTS:
% saves components from ICA to disk under dirs.ICA.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% By J.C.Swart, 2016/ J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Preprocessing/

%% Complete input settings:

% clear all; close all; clc
% dbstop if error

if ~exist('rootdir','var')
    rootdir = preprocessing_set_rootdir(); % '/project/3017042.02';
    fprintf('rootdir unspecified, assume %s\n',rootdir)
end

if ~exist('subVec','var')
    subVec = 1:36; 
    fprintf('subVec unspecified, assume %s\n',strjoin(string(subVec),', '));
end

%% Set directories:

% Add helper files to path:
fprintf('Add helper files to path\n');
addpath(fullfile(rootdir,'Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Helpers'));

% Set directories and parameters:
dirs    = set_dirs(rootdir);
par     = set_par();

%% Load channels to be excluded:

load(fullfile(dirs.data,'channel2repair.mat'));

%% Detect input files:

% Detect all available input folders (all subjects):
tmp         = dir(fullfile(dirs.preICA,'*.mat'));
fileList    = {tmp.name};
nInput      = length(fileList);
fprintf('Found data files from %d subjects\n',nInput);

% Extract subject numbers from folders available:
subNames    = str2double(extractBetween(fileList,8,10));

% Extract indices of selected subjects:
selIndices  = find(ismember(subNames,subVec));
nSub        = length(selIndices);

%% Loop over subjects:

for iSub = selIndices % iSub = 1;

    %% Create output file name:
    
    outputFile      = fullfile(dirs.ICA, sprintf('MRIpav_%03d_ICs.mat',iSub));
    fprintf('EEG output file is %s\n',outputFile);
    
    % Check if output file already exists:
    if exist(outputFile,'file')
        warning('Subject %03d: File %s already exists, skip subject\n', iSub, outputFile);
        continue;
    end

    %% Create input file name:
    
    inputFile       = fullfile(dirs.preICA, fileList{iSub}); % go to subject directory
    fprintf('EEG input file is %s\n',inputFile);
    
    %% Load data:
    
    fprintf('Subject %03d: Loading data ... (can take a few seconds)\n',iSub);
    load(inputFile) % 
    fprintf('Subject %03d: Loading completed\n',iSub);

    %% Determine valid and invalid channels:
    
    invalidChan             = cell2mat(channel2repair(iSub,2)); % previously discarded from visual inspection
    validChan               = setdiff(1:par.chan.nChan,[par.chan.ECG par.chan.HR par.chan.RESP invalidChan]); % not ECG, not HR, not RESP, not discarded
    %
    
    %% Run ICA (runica is default):
    
    fprintf('Subject %03d: Run ICA (takes several minutes)...\n', iSub)
    
    cfg                     = [];
    cfg.numcomponent        = length(validChan) - 1; % numbers channels - 1
    cfg.channel             = validChan; % only non-rejected channels
    comp                    = ft_componentanalysis(cfg, data);
    
    %% Save data: 

    fprintf('Subject %03d: Start saving data\n',iSub);
    
    hist.ICA.dirs       = dirs;
    hist.ICA.par        = par;
    
    save(outputFile,'comp','hist','-v7.3');
    clear data comp hist

    fprintf('Subject %03d: Finished saving data\n',iSub);

end % end iSub-loop.

fprintf('Done :-) \n');

end % end of function.
