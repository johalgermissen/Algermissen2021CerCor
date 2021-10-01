function [Y,job,betas] = taft_preprocess_load_EEG(job)

% [Y,job,betas] = taft_preprocess_load_EEG(job)
%
% Loads EEG data, reformat into matrix Y. 
% Also update job for which (good) trials are kept.
%
% INPUTS:
% job           = structure with settings to be specified via taft_preprocess_initialize_job
% .EEGdomain 	= string, either 'TF' (data in TF domain) or 'time' (data in time domain).
% .EEGfile      = string, full file path of subject-specific EEG data set to load.
% .selTrialsFile= string, name of file containing indices of selected trials to load. 
% .subID        = numeric scalar, subject ID.
%  
% OUTPUTS:
% Y             = a 3D (trials x channels x time) or 4D (trials x channels x frequency x time) matrix
% job           = structure with settings for creating TAfT object, new field added:
% .goodTrlIdx 	= numeric vector, indices of trials to include (after EEG trial rejection).
% betas         = trial-averaged Fieldtrip object, used as template to insert beta-map later while forwarding info on timing, channels, etc.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Adapted from Tobias Hauser (https://github.com/tuhauser/TAfT).
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/CueLockedAnalyses/TAfT/

fprintf('Load EEG data\n');

%% Load EEG file:

% Load file:
inputEEG   = load(job.EEGfile);

% Extract:
subEEG      = inputEEG.(char(job.EEGname));
    
% Remove electrode positions if existing field:
if isfield(subEEG,'elec')
    subEEG = rmfield(subEEG,'elec');
end

%% Retrieve good trials:

job.goodTrlIdx = subEEG.trialinfo; % indices of trials to-be-kept in fMRI

% Check if other trials provided for selection; if yes, intersect:
if isfield(job,'selTrialsFile')
    
    fprintf('Found file indicating selected %s: Load and adjust good trials\n',job.selTrials);
    selTrlIdx       = load(job.selTrialsFile); % load selected trials
    job.goodTrlIdx  = intersect(job.goodTrlIdx,selTrlIdx.selTrials); % intersection with EEG trials

    % Need to also re-select EEG data:
    cfg             = [];
    cfg.trials      = find(ismember(subEEG.trialinfo,job.goodTrlIdx));
    subEEG          = ft_selectdata(cfg,subEEG); % select time samples

end

%% Select EEG data:

% Select time:
if job.toi ~= [-inf inf] 
    cfg             = [];
    cfg.latency     = job.toi;
    subEEG          = ft_selectdata(cfg,subEEG); % select time samples
end

% Select frequencies:
if isfield(subEEG,'freq') && any(job.foi ~= [-inf inf])
    cfg             = [];
    cfg.frequency   = job.foi;
    subEEG          = ft_selectdata(cfg,subEEG); % select frequencies
end    

%% Create 3- or 4-D EEG object (channels, samples, (frequencies), trials): 

fprintf('Reshape EEG data\n');
Y                       = []; % initialize empty object

if strcmp(job.EEGdomain, 'time')

    cfg                 = [];
    betas               = ft_timelockanalysis(cfg,subEEG); % average over trials to obtain settings for final betascome object
    Y                   = zeros(size(betas.label,1),size(betas.time,2),size(subEEG.trialinfo,1)); % initialize
    for iTrial = 1:size(subEEG.trialinfo,1) 
        Y(:,:,iTrial)   = subEEG.trial{iTrial}; % trial becomes last dimension
    end
    
elseif strcmp(job.EEGdomain, 'TF')

    cfg                 = [];
    betas               = ft_freqdescriptives(cfg,subEEG); % average over trials to obtain settings for final betascome object
    Y                   = zeros(size(betas.label,1),size(betas.freq,2),size(betas.time,2),size(subEEG.trialinfo,1)); % initialize
    for iTrial = 1:size(subEEG.trialinfo,1) 
            Y(:,:,:,iTrial) = squeeze(subEEG.powspctrm(iTrial,:,:,:)); % trial becomes last dimension
    end
    
else
    error('Unknown EEGdomain data format')
end

end % end of function.
