function EEGfMRIPav_CueLocked_6_finalPP(rootdir,subVec)

% EEGfMRIPav_CueLocked_6_finalPP(rootdir,subVec)
% 
% Reject trials, compute surface LaPlacian filter (and at the same time
% interpolate bad channels). Final steps of pre-processing.
% 
% INPUTS:
% rootdir           = string, root directory of project
% subVec            = vector of integers, numbers of subjects to process
%
% OUTPUTS:
% saves data to disk under dirs.finalPP.
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

%% Load channels to be interpolate:

load(fullfile(dirs.data,'channel2repair.mat'));

%% Load list of complete channels:

load(fullfile(dirs.polhemus,'posLabels.mat'));

%% Detect input files:

% Detect all available input folders (all subjects):
tmp         = dir(fullfile(dirs.postICA,'*.mat'));
fileList    = {tmp.name};
nInput      = length(fileList);
fprintf('Found data files from %d subjects\n',nInput);

% Extract subject numbers from folders available:
subNames    = str2double(extractBetween(fileList,8,10));

% Extract indices of selected subjects:
selIndices  = find(ismember(subNames,subVec));
nSub        = length(selIndices);

%% Loop over subjects:

for iSub = selIndices % %iSub = 1;
    
    %% Create output file name:
    
    outputFile      = fullfile(dirs.finalPP, sprintf('MRIpav_%03d_finalPP.mat',iSub));
    fprintf('EEG output file is %s\n',outputFile);
    
    % Check if output file already exists:
    if exist(outputFile,'file')
        warning('Subject %03d: File %s already exists, skip subject\n', iSub, outputFile);
        continue;
    end

    %% Create input file name:
    
    inputFile       = fullfile(dirs.postICA, fileList{iSub}); % go to subject directory
    fprintf('EEG input file is %s\n',inputFile);
    
    %% Load data:
    
    fprintf('Subject %03d: Loading data ... (can take a few seconds)\n',iSub);
    load(inputFile) % 
    fprintf('Subject %03d: Loading completed\n',iSub);
    
    %% Load rejected trials:
    
    fprintf('Subject %03d: Load rejected trials\n',iSub)
    
    filenameRejTrials = fullfile(dirs.rejTrials,sprintf('MRIpav_%03d_rejTrials.mat',iSub));
    load(filenameRejTrials);

    % Provide information on screen:
    fprintf('Subject %03d: Removing rejected trials\n',iSub)

    % Remove post-IC rejected trials:
    if ~isempty(rejectedtrials)
        fprintf('Subject %03d: Removing %d out of %d trials\n',iSub,sum(rejectedtrials),size(data.trial,2))
        cfg                     = [];
        cfg.trials              = find(~rejectedtrials);
        data                    = ft_selectdata(cfg,data);

    end
    
    %% Interpolate channels (replace with NaN):

    invalidChan         = cell2mat(channel2repair(iSub,2)); % previously discarded from visual inspection

    % Provide information on screen:
    fprintf('Subject %03d: Create missing channels as NaNs\n',iSub)

    % Recover deleted channels (dropped before ICA) as NaN.
    cfg                 = [];
    cfg.method          = 'nan';
    cfg.missingchannel  = ft_channelselection(invalidChan,posLabels);
    data                = ft_channelrepair(cfg,data);  
    
    %% Compute surface Laplacian:

    % Provide information on screen:
    fprintf('Subject %03d: Compute scd via Surface Laplacian\n',iSub)

    % Calculate current source density:
    cfg                     = [];
    cfg.method              = 'spline';
    
    if ismember(iSub, [4 32:34]) % Polhemus data corrupted for sub 4, missing fors subs 32-34
        
        fprintf('Subject %03d: Invalid Polhemus data, load grand average instead\n',iSub);
        load(fullfile(dirs.polhemus,'PolhemusAvg.mat'));
        cfg.elec            avgGrand;
        
    else
        
        fprintf('Subject %03d: Found Polhemus data, loading...\n',iSub);
        cfg.elec            = ft_read_sens(fullfile(dirs.polhemus,sprintf('MRIpav_001_%03d.pos',iSub)));
        
        if iSub == 5
            cfg.elec.label(1:66) = posLabels([1:65,67]); % RESP skipped in Polhemus recording for subject 5
        else
            cfg.elec.label(1:67) = posLabels; % overwrite to have channel labels (e.g. FP1) instead of numbers    
        end
        
    end
    
    cfg.badchannel          = ft_channelselection(invalidChan,cfg.elec.label);
    scd                     = ft_scalpcurrentdensity(cfg,data);
    
    %% Save data: 

    fprintf('Subject %03d: Start saving data\n',iSub);
    
    hist.finalPP.dirs       = dirs;
    hist.finalPP.par        = par;
    
    save(outputFile,'scd','hist','-v7.3');
    clear data scd hist avgGrand

    % Save data before IC rejection:
    clear data
    
    fprintf('Subject %03d: Finished saving data\n',iSub);
    
end % end iSub-loop.

%% Check number of rejected trials:

nRejTrials = nan(nSub,1);
for iSub = 1:nSub % iSub = 1;
    fprintf('Start subject %03d\n',iSub);
    load(fullfile(dirs.rejTrials,sprintf('MRIpav_%03d_rejTrials.mat',iSub)));
    nRejTrials(iSub) = sum(rejectedtrials);
    clear rejectedtrials
end

fprintf('Number of rejected trials: M = %.03f, SD = %.03f, range %d - %d \n',mean(nRejTrials),std(nRejTrials),min(nRejTrials),max(nRejTrials));
% Number of rejected trials: M = 29.583, SD = 27.704, range 2 - 112 

fprintf('Done :-)\n');

end % end of function
