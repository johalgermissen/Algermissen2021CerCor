function EEGfMRIPav_CueLocked_8_time_create_BOLD(rootdir,subVec,lockSettings,baselineSettings,baselineTimings,...
    ROI2use)

% EEGfMRIPav_CueLocked_8_TF_grouplevel_create(rootdir,subVec,lockSettings,baselineSettings,baselineTimings,...
%     ROI2use,nBin)
% 
% Load EEG, rejected trials, behavior; apply baseline correction; sort
% trials into highest/lowest tertile based on trial-by-trial BOLD signal; 
% aggregate within bins, save. 
% 
% INPUTS:
% rootdir           = string, root directory of project.
% subVec            = vector of integers, numbers of subjects to include
% (default: 1:36).
% lockSettings      = string, type of time-locking before TF decomposition,
% either 'stimlocked' (default) or 'resplocked'.
% baselineSettings  = string, type of baseline correction, 'trend'
% (regression-based), 'ft' (with Fieldtrip's ft_timelockbaseline).
% baselineTimings   = vector of 1 or 2, timing of baseline correction
% which to use.
% TFtype            = type of TF decomposition, 'hanning' or 'morlet'.
% nFreq             = numeric, number of frequency bins decomposed
% (default: 15).
% ROI2use           = string, name of ROI for which to load BOLD data, use
% 'GLM1vmPFCValenceMan' (default).
%
% OUTPUTS:
% saves data to disk under dirs.TFgroup.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% By J.C.Swart, 2016/ J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel/

%% Complete input settings:

% clear all; close all; clc
% dbstop if error

if ~exist('rootdir','var')
    rootdir = grouplevel_set_rootdir(); % '/project/3017042.02';
    fprintf('rootdir unspecified, assume %s\n',rootdir)
end

if ~exist('subVec','var')
    subVec = 1:36; 
    fprintf('subVec unspecified, assume %s\n',strjoin(string(subVec),', '));
end

if ~exist('lockSettings','var')
    lockSettings = 'stimlocked'; % stimlocked resplocked
    fprintf('lockSettings unspecified, assume %s\n',lockSettings);
end

if ~exist('baselineSettings','var')
    baselineSettings = 'trend'; % trend, all, condition, trial, ft_trial
    fprintf('baselineSettings unspecified, assume %s\n',baselineSettings);
end

if ~exist('baselineTimings','var')
    baselineTimings = [-0.250 -0.050]; % 
    fprintf('baselineTimings unspecified, assume %s\n',strjoin(string(baselineTimings),' - '));
end

if ~exist('accSettings','var')
    accSettings = 'correct'; % correct incorrect bothAcc noAcc reqAct
    fprintf('accSettings unspecified, assume %s\n',accSettings);
end

if ~exist('ROI2use','var')
    ROI2use = 'GLM1vmPFCValenceMan';
    fprintf('ROI2use unspecified, assume %s\n',ROI2use);
end

%% Set directories:

% Add helper files to path:
fprintf('Add helper files to path\n');
addpath(fullfile(rootdir,'Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Helpers'));

% Set directories and parameters:
dirs        = set_dirs(rootdir);
par         = set_par();

if strcmp(accSettings,'noACC')
    allAcc = [1 0];

elseif strcmp(accSettings,'correct')
    allAcc = 1;

elseif strcmp(accSettings,'incorrect')
    allAcc = 0;

else
    fprintf('Unknown accuracy setting')
end

% Add SPM for upsampling:
addpath('/home/common/matlab/spm12');

% Add TAfT for upsampling:
addpath(fullfile(rootdir,'Analyses/EEG_Scripts/CueLockedAnalyses/TAfT'));

%% Detect input files:

% Detect all available input folders (all subjects):
tmp         = dir(fullfile(dirs.finalPP,'*_finalPP.mat'));
fileList    = {tmp.name};
nInput      = length(fileList);
fprintf('Found data files from %d subjects\n',nInput);

% Extract subject numbers from folders available:
subNames    = str2double(extractBetween(fileList,8,10));

% Extract indices of selected subjects:
selIndices  = find(ismember(subNames,subVec));
nSub        = length(selIndices);
nCond       = 2; % number of bins, fixed to 2 because only high/low

%% Create final file name:

outputFile = sprintf('EEGfMRIPav_time_%s_baseCor_%s_baseTime_%03dms_%03dms_BOLD_%s', ...
    lockSettings,baselineSettings,baselineTimings(1)*-1000,baselineTimings(end)*-1000,...
    ROI2use);

% Add directory and .mat ending:
outputFile = fullfile(dirs.timegroup,sprintf('%s.mat',outputFile));

fprintf('Selected file is %s\n',outputFile)

%% Loop over subjects:

if ~exist(outputFile,'file')

    subCount = 1; % relative subject number, i.e. also works if certain numbers are missing

    ERPdata     = cell(nSub,1);
    nTrialCond  = nan(nSub,nCond);
    
    for iSub = selIndices % %iSub = 1;

        %% 1) Load EEG data:

        fprintf('Subject %03d: Load time-domain EEG data\n',iSub);
        load(fullfile(dirs.finalPP,fileList{iSub}))
        scd = rmfield(scd,'elec'); % remove elec field
        fprintf('Subject %03d: Finished loading\n',iSub);
               
        %% 2) Load rejected trials:

        fprintf('Subject %03d: Load rejected trials\n',iSub)
        filenameRejTrials = fullfile(dirs.rejTrials,sprintf('MRIpav_%03d_rejTrials.mat',iSub));
        load(filenameRejTrials);

        %% 3) Load and recode behavior:
        
        fprintf('Subject %03d: Load behavioral data\n',iSub)
        behav = load_recode_behav(rootdir,iSub,rejectedtrials); 
               
        %% 4) Baseline correction:
        
        timeIdx = dsearchn(scd.time{1}',baselineTimings'); % select indices for baseline time window
        if length(timeIdx) > 1; timeIdx = timeIdx(1):timeIdx(end); end % extend to interval if start and end given
        
        if strcmp(baselineSettings,'ft')
            fprintf('Subject %03d: Subtract baseline per trial using Fieldtrip\n',iSub)
            cfg = [];
            cfg.baseline = baselineTimings; % -250 - 50 ms
            scd = ft_timelockbaseline(cfg,scd);
        elseif strcmp(baselineSettings,'trend')
            fprintf('Subject %03d: Subtract baseline per trial by detrending\n',iSub)
            scd.trial   = regress_baseline(scd.trial,scd.trialinfo,timeIdx(1),timeIdx(end),'subtraction');          
        end
        
        %% 5) Response lock the data (if selected):
           
        if strcmp(lockSettings,'resplocked')
            
            fprintf('Response-lock data\n');
            scd = redefine_resplock(par,behav,scd);
            
            % Remove any numerical imprecision by subtracting time point
            % closest to 0:
            for k = 1:numel(scd.time)
                ix = nearest(scd.time{k},0);
                scd.time{k} = scd.time{k} - scd.time{k}(ix);
            end
            
        end
        
        %% 6) Load BOLD data:
        
        % Initialize job with taft_initialize_job:

        job             = taft_preprocess_initialize_job('TF',iSub,ROI2use);
        job.goodTrlIdx  = 1:640; % select all trials
        X               = taft_preprocess_load_fMRI(job); % load upsampled, epoched, deconvolved data (one HRF amplitude per trial):
        
        % Select valid trials:
        X = X(scd.trialinfo);       
                
        %% 7) Sort trials into quantiles based on BOLD:
        
        % Determine quantiles:
        lowQ = quantile(X,1/3); % lowest third cutoff
        highQ = quantile(X,2/3); % highest third cutoff
        Xlow = X <= lowQ; % lowest third indices
        Xhigh = X >= highQ; % highest third indices
        Xsel = Xlow + 2*Xhigh; % low becomes 1, high becomes 2
        
        %% 6) Loop over conditions:
        
        for iCondi = 1:nCond
            cfg                             = [];
            cfg.trials                      = find(Xsel==iCondi);
            nTrialCond(iSub,iCondi)         = length(cfg.trials);
            fprintf('Subject %03d condition %d of %s: Found %d trials \n',iSub,iCondi,char(ROI2use),length(cfg.trials));
            if ~isempty(cfg.trials) % if any trials in this condition
                ERPdata{iSub,iCondi}        = ft_timelockanalysis(cfg,scd);
            else % otherwise fill in NaNs
                fprintf('Subject %03d condition %d: Found no trials, store NaNs! \n',iSub, iCondi);
                ERPdata{iSub,iCondi}.time   = scd.time{1}; % timing of first trial
                ERPdata{iSub,iCondi}.label  = scd.label;
                ERPdata{iSub,iCondi}.avg    = zeros(length(scd.label),length(scd.time{1}));
                ERPdata{iSub,iCondi}.var    = zeros(length(scd.label),length(scd.time{1}));
                ERPdata{iSub,iCondi}.dof    = zeros(length(scd.label),length(scd.time{1}));
                ERPdata{iSub,iCondi}.dimord = 'chan_time';
                ERPdata{iSub,iCondi}.cfg    = scd.cfg;                         
            end % end of check for empty trials
        end
                    
        fprintf('Subject %03d: Finished\n',iSub);
        clear scd behav 
        
    end % end for-loop
    
    %% Save data for all subjects:
    
    fprintf('Save file under %s\n',outputFile)
    save(outputFile,'ERPdata','nTrialCond','-v7.3');
    fprintf('Finished :-)\n')

end % end if ~exist(outputFile)

fprintf('Done :-)\n');

end % end of function.
