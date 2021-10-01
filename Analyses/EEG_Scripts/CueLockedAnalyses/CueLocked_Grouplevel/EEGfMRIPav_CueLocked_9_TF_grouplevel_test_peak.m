% EEGfMRIPav_CueLocked_9_TF_grouplevel_test_peak

% This is an interactive script---execute it step-by-step.
% Set the appropriate job settings for loading TF data aggregated across
% subjects, then perform tests comparing peak height and latency of
% conditions.
% 
% EXPLANATION OF SETTINGS:
% rootdir               = string, root directory of project.
% job                   = cell, created via TF_update_job.m, needs at least
% fields:
%   .nSub               = integer, number of subjects.
%   .dirs               = cell, directories
%   .sub2exclude        = numeric vector, numbers of subjects to exclude.
%   .lockSettings       = string, type of event-locking, 'stimlocked' or
%   'resplocked'.
%   .baselineSettings   = string, type of baseline correction, 'trend'
%   (regression-based), 'all' (grand mean oer subject), 'condition'
%   (grand mean separately per condition), 'trial' (per trial, own
%   implementation), 'ft_trial' (with Fieldtrip's ft_baselinecorrect).
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
%   .TFtype             = type of TF decomposition, 'hanning' or 'morlet'.
%   .nFreq              = numeric, number of frequency bins decomposed
%   (default: 15).
%   .chanArea           = string, area of channels to select, 'midfrontal',
%   'frontal', 'leftmotor' or 'rightmotor'.
% 	.band               = string, frequency band to select, 'theta' or 
%   'alpha' or 'beta' or 'broad' 
%   .contrastType       = string, contrast to be used, either 
%   'Congruency', 'Go', 'GoLeftRight', 'GoValence', 'GoAccuracy', 
%   'Accuracy', 'CongruentAccuracy' 'IncongruentAccuracy'.
% 	.stimERPcor         = Boolean, read data corrected for stimulus-locked 
% 	ERP (true) or not (false), default false.
% 	.respERPcor         = Boolean, read data corrected for response-locked 
% 	ERP (true) or not (false), default false.
% 	.bin.Type           = string, split up by what variable ('RT',
%   'vmPFCValenceMan'), optional.
%   .bin.Num            = numeric, number of bins (2 or 3), optional.
%   .bin.Type           = string, variable after which bins are created,
%   'RT' or 'GLM1vmPFCValenceMan', optional.
%
% OUTPUTS:
% no outputs, just plots.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% By J.C.Swart, 2016/ J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel/

%% Set directories:

rootdir     = grouplevel_set_rootdir(); % '/project/3017042.02';

% Add helper files to path:
fprintf('Add helper files to path\n');
addpath(fullfile(rootdir,'Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Helpers'));
addpath(fullfile(rootdir,'Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel'));

% Set directories and parameters:
dirs        = set_dirs(rootdir);
par         = set_par();

%% Initialize job:

job = []; % initialize empty job

job.dirs = dirs; % add directories

% Data settings:
job.nSub                = 36; % necessary for validSubs before loading data
job.sub2exclude         = [];
% job.sub2exclude         = [11 12 15 23 25 26 30]; % exclusion in TAfT-outcome-locked

job.lockSettings        = 'stimlocked'; % stimlocked resplocked

job.baselineSettings    = 'trend'; % 'all' or 'condition' or 'trial' or 'ft_trial' or 'trend'    
job.baselineTimings     = [-0.250 -0.050];

job.responseSettings    = 'Go'; % 'Go' or 'Hand'  
job.actionSettings      = 'exeAct'; % 'reqAct' or 'exeAct'
job.accSettings         = 'correct'; % 'bothAcc' or 'correct' or 'incorrect' or 'noAcc' or 'reqAct'
% for bothAcc, always use reqAct

job.TFtype              = 'hanning'; % hanning morlet
job.nFreq               = 15; % 64 15 6

% Channel settings:
job.chanArea            = 'midfrontal'; % 'midfrontal' or 'frontal' or 'central' or 'leftmotor' or 'rightmotor'

% Band settings:
job.band                = 'theta'; % 'theta' or 'alpha' or 'beta' or 'broad' 

% Contrast of interest:
job.contrastType        = 'Go'; % 'Congruency' or 'Go' or 'GoLeftRight' or 'GoValence' or 'GoAccuracy' or 'Accuracy' or 'CongruentAccuracy' or 'IncongruentAccuracy'

% ----------------------------------------------------------------------- %
% ERP-corrected:
job.stimERPcor = false;
job.respERPcor = false;

% ----------------------------------------------------------------------- %
% Split up in bins:
job.bin.Type    = 'RT'; % RT GLM5FvmPFCValenceMan
job.bin.Num     = 3; % only relevant for RT

% Load and prepare data:

job         = TF_update_job(job); % initialize job
[job, data] = TF_load_data(job); % load data
[job, data] = TF_prepare_generic_data(job,data); % prepare generic objects

job         = TF_update_job(job); % update job
[job, data] = TF_prepare_contrast_data(job,data); % prepare objects specific to contrast

%% TEST FOR TIMING DIFFERENCES:

% CIs for Cohen's d: use Cousineau-Morey method to compute SE, 1.96*SE
% (one-sided: -Inf - 1.64*SE)
fprintf('Extract peak per condition per subject\n');

% Retrieve indices of maximum per subject per condition:
peakIdx     = nan(size(data.SubCondTime,1),size(data.SubCondTime,2));
peakTime    = nan(size(data.SubCondTime,1),size(data.SubCondTime,2));
peakTHeight = nan(size(data.SubCondTime,1),size(data.SubCondTime,2));

% Minimal and maximum time for extraction:
if strcmp(job.lockSettings,'stimlocked')
    iTimeStart = find(round(data.mu.time,3)==0.5); % 41
    iTimeStop = find(round(data.mu.time,3)==1.3); % 93
elseif strcmp(job.lockSettings,'resplocked')
    iTimeStart = find(round(data.mu.time,3)==-0.5); % 41
    iTimeStop = find(round(data.mu.time,3)==0.5); % 93
else
    error('Unknown lock setting')
end

% Extract peak height and latency per condition per subject:
for iSub = 1:size(data.SubCondTime,1)
    for iCond = 1:size(data.SubCondTime,2)
        maxValue                = max(data.SubCondTime(iSub,iCond,iTimeStart:iTimeStop));
        peak                    = iTimeStart + find(data.SubCondTime(iSub,iCond,iTimeStart:iTimeStop)==maxValue);
%         peak                    = find(data.SubCondTime(iSub,iCond,:)==maxValue);
        peakIdx(iSub,iCond)     = peak(1);
        peakTime(iSub,iCond)    = data.mu.time(peakIdx(iSub,iCond));
        peakHeight(iSub,iCond)  = data.SubCondTime(iSub,iCond,peakIdx(iSub,iCond));
    end
end

%% Go2Win vs. Go2Avoid:
% Load regular 4 condition data:

fprintf('Go2Win: M = %.03f; Go2Avoid: M = %.03f\n',mean(peakTime(:,1)),mean(peakTime(:,2)));
fprintf('Go2Win: STD = %.03f; Go2Avoid: STD = %.03f\n',std(peakTime(:,1)),std(peakTime(:,2)));
dif = peakTime(:,2)-peakTime(:,1); % Go2Avoid - Go2Win
[~,p,~,stats] = ttest(dif); % index
fprintf('t(%d) = %.03f, p = %.03f, d = %.02f\n',stats.df,stats.tstat,p,mean(dif)/std(dif));

[~] = permutation_test(peakIdx(:,2),peakIdx(:,1),1000);

%% Permutation tests:

% Permutation test: 2 bins:
% [~] = permutation_test(peakIdx(:,2),peakIdx(:,1),1000);
% [~] = permutation_test(peakIdx(:,4),peakIdx(:,3),1000);

% Permutation test: 3 bins:
[~] = permutation_test(peakIdx(:,2),peakIdx(:,1),1000);
[~] = permutation_test(peakIdx(:,3),peakIdx(:,2),1000);
[~] = permutation_test(peakIdx(:,5),peakIdx(:,4),1000);
[~] = permutation_test(peakIdx(:,6),peakIdx(:,5),1000);

%% t-test (both indices and timing give exactly the same result):

% 2 bins:
% dif = peakIdx(:,2)-peakIdx(:,1);
% dif = peakIdx(:,4)-peakIdx(:,3);
% [~,p,~,stats] = ttest(peakIdx(:,4)-peakIdx(:,3)); % index
% fprintf('t(%d) = %.03f, p = %.03f, d = %.02f\n',stats.df,stats.tstat,p,mean(dif)/std(dif)); % divide  by 2 because one-sided

% 3 bins:
dif = peakIdx(:,2)-peakIdx(:,1);
[~,p,~,stats] = ttest(dif); % index
fprintf('t(%d) = %.03f, p = %.03f, d = %.03f\n',stats.df,stats.tstat,p,mean(dif)/std(dif)); % divide  by 2 because one-sided

dif = peakIdx(:,3)-peakIdx(:,2);
[~,p,~,stats] = ttest(dif); % index
fprintf('t(%d) = %.03f, p = %.03f, d = %.03f\n',stats.df,stats.tstat,p,mean(dif)/std(dif)); % divide  by 2 because one-sided

% dif = peakIdx(:,3)-peakIdx(:,1);
% [~,p,~,stats] = ttest(dif); % index
% fprintf('t(%d) = %.02f, p = %.03f, d = %.03f\n',stats.df,stats.tstat,p,mean(dif)/std(dif)); % divide  by 2 because one-sided

dif = peakIdx(:,5)-peakIdx(:,4);
[~,p,~,stats] = ttest(dif); % index
fprintf('t(%d) = %.03f, p = %.03f, d = %.03f\n',stats.df,stats.tstat,p,mean(dif)/std(dif)); % divide  by 2 because one-sided

dif = peakIdx(:,6)-peakIdx(:,5);
[~,p,~,stats] = ttest(dif); % index
fprintf('t(%d) = %.03f, p = %.03f, d = %.03f\n',stats.df,stats.tstat,p,mean(dif)/std(dif)); % divide  by 2 because one-sided

% dif = peakIdx(:,6)-peakIdx(:,4);
% [~,p,~,stats] = ttest(dif); % index
% fprintf('t(%d) = %.03f, p = %.03f, d = %.03f\n',stats.df,stats.tstat,p,mean(dif)/std(dif)); % divide  by 2 because one-sided

%% Loop to systematically contrast relevant conditions with each other:

% Stimulus-locked:
startIdx    = [2 3 5 6];
endIdx      = startIdx -1;

% Response-locked:
startIdx    = [2 3 3 5 6 6];
endIdx      = [1 2 1 4 5 4];

% Accuracy:
% startIdx    = [1 2];
% endIdx      = [5 6];

for iTime = 1:length(startIdx)
    fprintf('For index %d against %d:\n',startIdx(iTime),endIdx(iTime));
    
    % Compute difference vector:
    testVec           = peakIdx(:,startIdx(iTime))-peakIdx(:,endIdx(iTime)); % test peak index
%     testVec         = peakTime(:,startIdx(iTime))-peakTime(:,endIdx(iTime)); % test peak time (same as index)
%     testVec         = peakHeight(:,startIdx(iTime))-peakHeight(:,endIdx(iTime)); % test peak height
    
    % Perform test:
    [~,p,~,stats]   = ttest(testVec); % perform t-test
    d               = mean(testVec)/std(testVec); % compute Cohen's d
    
    % correct if one-sided:
    if strcmp(job.lockSettings,'stimlocked')
        p = p/2; % one-tailed
    end
    
    % Print to console:
    fprintf('Result is t(%d) = %.03f, p = %.03f, d = %0.3f\n',stats.df,stats.tstat,p,d);
end

%% Permutation test:

for iTime = 1:length(startIdx)
    fprintf('For index %d against %d:\n',startIdx(iTime),endIdx(iTime)); % test peak index
    p = permutation_test(peakIdx(:,startIdx(iTime)),peakIdx(:,endIdx(iTime)),1000); % test peak time (same as index)
%     p = permutation_test(peakHeight(:,startIdx(iTime)),peakHeight(:,endIdx(iTime)),1000); % test peak height
end

% find(testVec==min(testVec))

% Peak latency stimulus-locked:
% For index 2 against 1:
% Result is t(35) = 1.745, p = 0.045, d = 0.291
% For index 3 against 2:
% Result is t(35) = 2.577, p = 0.007, d = 0.430
% For index 5 against 4:
% Result is t(35) = 1.662, p = 0.053, d = 0.277
% For index 6 against 5:
% Result is t(35) = 5.220, p = 0.000, d = 0.870

% Peak latency response-locked:
% For index 2 against 1:
% Result is t(35) = -0.784, p = 0.438, d = -0.131
% For index 3 against 2:
% Result is t(35) = 0.896, p = 0.376, d = 0.149
% For index 3 against 1:
% Result is t(35) = 0.135, p = 0.893, d = 0.023
% For index 5 against 4:
% Result is t(35) = 0.014, p = 0.989, d = 0.002
% For index 6 against 5:
% Result is t(35) = 0.347, p = 0.731, d = 0.058
% For index 6 against 4:
% Result is t(35) = 0.245, p = 0.808, d = 0.041

% Peak height response-locked:
% For index 2 against 1:
% Result is t(35) = -1.138, p = 0.263, d = -0.190
% For index 3 against 2:
% Result is t(35) = 1.003, p = 0.323, d = 0.167
% For index 3 against 1:
% Result is t(35) = 0.206, p = 0.838, d = 0.034
% For index 5 against 4:
% Result is t(35) = -0.017, p = 0.987, d = -0.003
% For index 6 against 5:
% Result is t(35) = 0.195, p = 0.846, d = 0.033
% For index 6 against 4:
% Result is t(35) = 0.195, p = 0.846, d = 0.033

% END.
