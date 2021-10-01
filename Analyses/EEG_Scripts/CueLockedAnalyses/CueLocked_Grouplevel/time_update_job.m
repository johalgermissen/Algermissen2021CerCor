function [job] = time_update_job(job)

% [job] = TF_update_job(job)
%
% Complete settings for job on time-domain data to create.
% (data sets to create; plotting settings).
% 
% INPUTS:
% job                   = cell, created via time_update_job.m, needs at 
% least fields:
%   .nSub               = integer, number of subjects.
%   .sub2exclude        = numeric vector, numbers of subjects to exclude.
%   .lockSettings       = string, type of event-locking, 'stimlocked' or
%   'resplocked'.
%   baselineSettings  = string, type of baseline correction, 'trend'
%   (regression-based), 'ft' (with Fieldtrip's ft_timelockbaseline).
%   .baselineTimings    = vector of 1 or 2, timing of baseline correction
%   for which to load data.
%   .responseSettings   = string, type of response setting for which
%   conditions are split, 'Go' (Go/ NoGo) or 'Hand' (Left Go/ Right Go/
%   NoGo.)
%   .actionSettings     = string, type of action setting for which
%   conditions are split, 'exeAct' (executed action) or 'reqAct' (required
%   action).
%   .accSettings        = string, type of accuracy setting for which 
%   conditions are split, 'correct', 'incorrect', 'bothAcc' (both
%   accuracies separately), 'noAcc' (no distinction), 'reqAct' (required
%   action).
%   .chanArea           = string, area of channels to select, 
%   'frontopolar', 'frontal', 'midfrontal', 'leftfrontal', 'midcentral', 
%   'central', 'leftmotor', 'rightmotor', 'centroparietal', 'leftparietal', 
%   'rightparietal', 'posterior', 'midoccipital', 'rightoccipital', 
%   'occipital'.
%   .contrastType       = string, contrast to be used: 'Congruency', 'Go',
%   'Valence', 'GoValence', 'GoAccuracy', 'GoLeftRight', 'Accuracy',
%   'CongruentAccuracy', or 'IncongruentAccuracy'.
% 	.bin.Type           = string, split up by what variable
% 	('RT', 'GLM5FvmPFCValence'), optional.
%   .bin.Num            = numeric, number of bins (2 or 3), optional (only
%   for RTs, for BOLD always tertile split and outer bins retained).
% 
% OUTPUTS:
% job                   = cell, with the following fields added:
%   .inputFileName      = name of file to load: 
%   .layout             = string, cap to use for topoplots.
%   .timetiming         = numeric vector of 2 elements, timing for TF
%   plots.
%   .sigTime            = numeric vector of 2 elements, timing for
%   topoplots.
%   .channels           = vector of strings, selected channel names.
%   .chanName           = string, name of selected channel area
%   capitalized.
%   .invalidSubs        = numeric vector, subjects to exclude.
%   .validSubs          = numeric vector, subjects to include.
%   .nValidSubs         = numeric, number of subjects to include.
%   .lockName           = string, name of event-locking, fully written out.
%   .allACC             = numeric vector, values of behavioral accuracy to
%   loop over.
%   .nCond              = numeric, number of conditions.
%   .condOrder          = numeric vector, order of conditions to plot.
%   .colMat             = matrix, colors of lines to plot per condition.
%   .lineStyle          = vector of strings, line types to plot per
%   condition.
%   .condNames          = vector of strings, condition labels to use in
%   plots.
%   .plotName           = string, important settings to be included in plot
%   names.
%   .twoLineColor       = 2 x 3 matrix, line colors to use for twoLinePlot.
%   .twoLineLabels      = vector of 2 strings, labels for conditions.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel/

fprintf('Update and complete job settings\n')

%% File to load:

% Default input file:
job.inputFileName = sprintf('EEGfMRIPav_time_%s_baseCor_%s_baseTime_%03dms_%03dms', ...
    job.lockSettings,job.baselineSettings,job.baselineTimings(1)*-1000,job.baselineTimings(end)*-1000);

% If binning (RT or BOLD) involved:
if isfield(job,'bin')
    if strcmp(job.bin.Type,'RT')
        job.inputFileName   = [job.inputFileName sprintf('_RT%d', job.bin.Num)]; % add number of bins (response/valence/accuracy irrelevant)
        job.contrastType    = 'RT';
    else
        job.inputFileName   = [job.inputFileName sprintf('_BOLD_%s', job.bin.Type)]; % add ROI name (always 2 bins) (response/valence/accuracy irrelevant)
        job.contrastType    = 'BOLD';
    end
else
    % Standard input file:
    job.inputFileName = [job.inputFileName sprintf('_resp_%s_%s_%s', ...
    job.responseSettings,job.actionSettings,job.accSettings)];
end

% Add final ending:
job.inputFileName = sprintf('%s.mat',job.inputFileName); % add final ending

%% Plotting options:

% Color map:
set(0, 'DefaultFigureColormap', jet(64))    

% Cap layout:
job.layout   = 'easycapM11.mat'; %

% Timing of time plot:
job.timetiming = [-0.25 1.0];

%% Two line plots: color, labels, time of significant differences:

% Timing for topoplot:
if isfield(job,'sigTime') % save in case it is given as input; return back in again later
    sigTime = job.sigTime;
end
job.sigTime = []; % initialize/delete

% Check for contrast type:
if isfield(job,'bin') 

    job.twoLineColor    = [1 0 0; 0 0 1]; % red, blue 

    if strcmp(job.bin.Type,'RT')
        job.twoLineLabels   = {'Fast RTs','Slow RTs'};
        
    elseif strcmp(job.bin.Type,'GLM5FvmPFCValenceMan')
        job.twoLineLabels   = {'High vmPFC','Low vmPFC'};
    
    else
        job.twoLineLabels   = {'High BOLD','Low BOLD'};
    end
    
elseif strcmp(job.contrastType,'Congruency')

    job.twoLineColor    = [1 0.5 0; 0.37 0.15 0.54]; % orange, purple 
    job.twoLineLabels   = {'Incongruent','Congruent'};

    job.sigTime         = []; % data.line2 > data.line1

elseif strcmp(job.contrastType,'Go')

    job.twoLineColor    = [1 0 0; 0 0 1]; % red, blue 
    job.twoLineLabels   = {'Go','NoGo'};

    if strcmp(job.lockSettings,'stimlocked')
        job.sigTime = [];

    elseif strcmp(job.lockSettings,'resplocked')
        job.sigTime = [];

    else
        error('Invalid locking setting')
    end

elseif strcmp(job.contrastType,'Valence')

    job.twoLineColor    = [0 0.6 0.2;0.8 0 0]; % green, red
    job.twoLineLabels = {'Win','Avoid'};

    if strcmp(job.lockSettings,'stimlocked')
        job.sigTime = [];

    elseif strcmp(job.lockSettings,'resplocked')
        job.sigTime = [];

    else
        error('Invalid locking setting')
    end

elseif strcmp(job.contrastType,'GoValence')

job.twoLineColor    = [0 0.6 0.2;0.8 0 0]; % green, red
job.twoLineLabels   = {'Go2Win','Go2Avoid'};

    if strcmp(job.lockSettings,'stimlocked')
        job.sigTime = [];

    elseif strcmp(job.lockSettings,'resplocked')
        job.sigTime = [];

    else
        error('Invalid locking setting')
    end

elseif strcmp(job.contrastType,'GoLeftRight')

    job.twoLineColor    = [1 0 0; 0 0 1]; % red, blue 
    job.twoLineLabels   = {'Left','Right'};

    if strcmp(job.lockSettings,'stimlocked')

        job.sigTime = [];

    elseif strcmp(job.lockSettings,'resplocked')

        job.sigTime = [];

    else
        error('Invalid locking setting')
    end

elseif strcmp(job.contrastType,'Accuracy')

    job.twoLineColor    = [1 0 0; 0 0 1]; % red, blue 
    job.twoLineLabels   = {'Correct','Incorrect'};

elseif strcmp(job.contrastType,'GoAccuracy')

    job.twoLineColor    = [1 0 0; 0 0 1]; % red, blue 
    job.twoLineLabels   = {'Go Correct','Go Incorrect'};

elseif strcmp(job.contrastType,'CongruentAccuracy')

    job.twoLineColor    = [1 0 0; 0 0 1]; % red, blue 
    job.twoLineLabels = {'Congruent correct','Congruent incorrect'};

else
    warning('unknown sigTime')
end

if exist('sigTime','var') % if previously saved:
    job.sigTime = sigTime; % add back in again
end

%% Channels and labels:

if strcmp(job.chanArea,'frontopolar')
    job.channels = {'Fpz','Fp1','Fp2'};
    job.chanName = 'Frontopolar';
    
elseif strcmp(job.chanArea,'frontal')
    job.channels = {'Fpz','AFz','Fz'};
    job.chanName = 'Frontal';
    
elseif strcmp(job.chanArea,'midfrontal')
    job.channels = {'Fz','FCz','Cz'};
    job.chanName = 'Midfrontal';
    
elseif strcmp(job.chanArea,'leftfrontal')
    job.channels = {'F1','F3','FCz','FC1','FC3','Cz','C1','C3'};
    job.chanName = 'Leftfrontal';

elseif strcmp(job.chanArea,'midcentral')
    job.channels = {'FCz','Cz'};
    job.chanName = 'Midcentral';
    
elseif strcmp(job.chanArea,'central')
    job.channels = {'FCz','Cz','CPz','FC1','FC2','C1','C2','CP1','CP2'};
    job.chanName = 'Central';
    
elseif strcmp(job.chanArea,'leftmotor')
    job.channels = {'C3','C5','CP3','CP5'};
    job.chanName = 'Left motor';
    
elseif strcmp(job.chanArea,'rightmotor')
    job.channels = {'C4','C6','CP4','CP6'};
    job.chanName = 'Right motor';
    
elseif strcmp(job.chanArea,'centroparietal')
    job.channels = {'CPz','Pz'};
    job.chanName = 'Centroparietal';
    
elseif strcmp(job.chanArea,'leftparietal')
    job.channels = {'CP1','CP3','P1','P3'};
    job.chanName = 'Left parietal';
    
elseif strcmp(job.chanArea,'rightcentral')
    job.channels = {'Cz','C2','CPz','CP2'};
    job.chanName = 'Right central';
    
elseif strcmp(job.chanArea,'posterior')
    job.channels = {'Oz','O1','O2','O3','O4','POz','PO1','PO2'};
    job.chanName = 'Posterior';
    
elseif strcmp(job.chanArea,'midoccipital')
    job.channels = {'Pz','POz','Oz'};
    job.chanName = 'Midoccipital';
    
elseif strcmp(job.chanArea,'rightoccipital')
    job.channels = {'P4','PO4','O2'};
    job.chanName = 'Right occipital';
    
elseif strcmp(job.chanArea,'occipital')
    job.channels = {'Oz','O1','O2','POz','PO3','PO4'};
    job.chanName = 'Occipital';
    
else
    fprintf('No channels set\n')
end

%% Invalid subjects:

% Update which subjects to drop:
job.invalidSubs = []; 
% job.invalidSubs = [1 11 15 19 21 25 26]; % outliers in TAfT and fMRI

% If any subjects given manually, via job.sub2exclude
if ~isempty(job.sub2exclude)
    fprintf('Will exclude manually specified subjects %s\n',num2str(job.sub2exclude));
    job.invalidSubs = job.sub2exclude;
end

if ~(isfield(job,'bin')) % unless split into bins:
    if strcmp(job.responseSettings,'Go') && strcmp(job.actionSettings,'exeAct') && (strcmp(job.accSettings,'bothAcc') || strcmp(job.accSettings,'reqAct') || strcmp(job.accSettings,'incorrect'))
        job.invalidSubs = [job.invalidSubs 11 24];

    elseif strcmp(job.responseSettings,'Go')
        job.invalidSubs = job.invalidSubs;

    elseif strcmp(job.responseSettings,'Hand') && strcmp(job.actionSettings,'exeAct') && (strcmp(job.accSettings,'bothAcc') || strcmp(job.accSettings,'reqAct'))
        job.invalidSubs = [job.invalidSubs 6 11 18 24 31];

    elseif strcmp(job.responseSettings,'Hand') % && strcmp(job.accSettings,'correct')
        job.invalidSubs = [job.invalidSubs 6 18 31];

    else
        error('Invalid response setting')
    end
end

% Remove repetitions:
job.invalidSubs = unique(job.invalidSubs); % remove repetitions

% Infer valid subjects:
job.validSubs   = setdiff(1:job.nSub,job.invalidSubs);
job.nValidSubs  = length(job.validSubs);

%% Label for lock settings:

if strcmp(job.lockSettings,'stimlocked')
    job.lockName = 'stimulus-locked';

elseif strcmp(job.lockSettings,'resplocked')
    job.lockName = 'response-locked';
    
else
    error('Unknown locking setting')
end

%% Complete settings based on accuracy setting:

if ~(isfield(job,'bin')) % unless split into bins:
    if strcmp(job.accSettings,'bothAcc') || strcmp(job.accSettings,'reqAct') || strcmp(job.accSettings,'noAcc')
        job.allAcc = [1 0];

    elseif strcmp(job.accSettings,'correct')
        job.allAcc = 1;

    elseif strcmp(job.accSettings,'incorrect')
        job.allAcc = 0;

    else
        fprintf('Unknown accuracy setting')
    end
end

%% Complete settings based on response setting:

% If RT or BOLD bins:
if isfield(job,'bin')
    
    if strcmp(job.bin.Type,'RT')
        
        if job.bin.Num == 2

            job.nCond       = 4;
            job.condOrder   = 1:4;    
            job.colMat      = [0 0.6 0.2;0.8 0 0;0 0.6 0.2;0.8 0 0];
            job.lineStyle   = {'-','-','--','--'};
            job.condNames = {'Go2Win fast','Go2Win slow','Go2Avoid fast','Go2Avoid slow'};

        elseif job.bin.Num == 3

            job.nCond       = 6;
            job.condOrder   = 1:6;
            job.colMat      = [0 0.6 0.2;0.8 0 0;0 0.6 0.2;0.8 0 0;0 0.6 0.2;0.8 0 0];
            job.lineStyle   = {'-','-',':',':','--','--'};
            job.condNames = {'Go2Win fast','Go2Win medium','Go2Win slow','Go2Avoid fast','Go2Avoid medium','Go2Avoid slow'};
            
        else
            error('Bin number %d not implemented yet',job.bin.Num);
        end
        
    else % else BOLD
        
            job.nCond       = 2;
            job.condOrder   = 1:2;
            job.colMat      = [1 0 0; 0 0 1];
            job.lineStyle   = {'-','-'};
            job.condNames = {'High BOLD','Low BOLD'};
        
    end
    
% If Go/NoGo distinction only:
elseif strcmp(job.responseSettings,'Go')
    
    if strcmp(job.accSettings,'correct') || strcmp(job.accSettings,'incorrect') % only one accuracy level
        
        job.nCond       = 4;
        job.condOrder   = 1:4;
        job.condNames   = {'Go2Win','Go2Avoid','NoGo2Win','NoGo2Avoid'};
        job.colMat      = [0 0.6 0.2;0.8 0 0;0 0.6 0.2;0.8 0 0];
        job.lineStyle   = {'-','-','--','--'};
        
    elseif strcmp(job.accSettings,'bothAcc') % irrespective of accuracy
        
        job.nCond       = 8;
        job.condOrder   = 1:8;
        job.condNames   = {'Go2Win correct','Go2Avoid correct','NoGo2Win correct','NoGo2Avoid correct','Go2Win incorrect','Go2Avoid incorrect','NoGo2Win incorrect','NoGo2Avoid incorrect'};
        job.colMat      = [0 0.6 0.2;0.8 0 0;0 0.6 0.2;0.8 0 0;.93 .69 0.13;0 0 1;.93 .69 0.13;0 0 1];
        job.lineStyle   = {'-','-','--','--','-','-','--','--'};
        
    elseif strcmp(job.accSettings,'reqAct') % if required action
        
        job.nCond       = 4;
        job.condOrder   = 1:4;
        job.condNames   = {'NoGo2Win -- NoGo required','Go2Avoid -- Go required','NoGo2Win -- Go required','Go2Avoid -- NoGo required'};
        job.colMat      = [0 0.6 0.2;0.8 0 0;0 0.6 0.2;0.8 0 0];
        job.lineStyle   = {'-','-','--','--'};
        
    else
        error('Unknown accuracy setting')
    end
    
% If leftGo/ right Go/ NoGo distinction only:
elseif strcmp(job.responseSettings,'Hand')
    
    if strcmp(job.accSettings,'correct') || strcmp(job.accSettings,'incorrect') % only one accuracy level

        job.nCond       = 6;
        job.condOrder   = 1:6;
        job.condNames   = {'LeftGo2Win','LeftGo2Avoid','RightGo2Win','RightGo2Avoid','NoGo2Win','NoGo2Avoid'};
        job.colMat      = [0 0.6 0.2;0.8 0 0;0 0.6 0.2;0.8 0 0;0 0.6 0.2;0.8 0 0];
        job.lineStyle   = {'-','-',':',':','--','--'};
        
    elseif strcmp(job.accSettings,'bothAcc') % irrespective of accuracy

        job.nCond       = 12;
        job.condOrder   = 1:12;
        job.condNames   = {'LeftGo2Win correct','LeftGo2Avoid correct','RightGo2Win correct','RightGo2Avoid correct','NoGo2Win correct','NoGo2Avoid correct','LeftGo2Win incorrect','LeftGo2Avoid incorrect','RightGo2Win incorrect','RightGo2Avoid incorrect','NoGo2Win incorrect','NoGo2Avoid incorrect'};
        % Decide on colors...?
        % Decide on line types...?
        
    else
        error('Unknown accuracy setting')
    end
else
    fprintf('Unknown response setting')
end

%% Handle for saving:

% Standard:
job.plotName = sprintf('%s_%s_baseCor_%s_baseTime_%d_%d',...
    job.contrastType,job.lockSettings,job.baselineSettings,job.baselineTimings(1)*1000,job.baselineTimings(end)*1000); % start

% If binned: add accuracy, type and number of bins:
if isfield(job,'bin')
    job.plotName = sprintf('%s_%s_%s%d',...
        job.plotName,job.accSettings,job.bin.Type,job.bin.Num);
else % else add response, action, accuracy
    job.plotName = sprintf('%s_%s_%s_%s',...
        job.plotName,job.responseSettings,job.actionSettings,job.accSettings);
end

% Add channels:
job.plotName = sprintf('%s_%s',...
    job.plotName,job.chanArea);

end % end of function.
