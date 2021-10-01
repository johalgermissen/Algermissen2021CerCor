% EEGfMRIPav_CueLocked_9_time_grouplevel_test

% This is an interactive script---execute it step-by-step.
% Set the appropriate job settings for loading TF data aggregated across
% subjects, then perform permutation test contrasting two conditions or
% just simple t-test (by hand or with Fieldtrip) contrasting two
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
%   .chanArea           = string, area of channels to select, 'midfrontal',
%   'frontal', 'leftmotor' or 'rightmotor'.
%   .contrastType       = string, contrast to be used, either 
%   'Congruency', 'Go', 'GoLeftRight', 'GoValence', 'GoAccuracy', 
%   'Accuracy', 'CongruentAccuracy' 'IncongruentAccuracy'.
% 	.stimERPcor         = Boolean, read data corrected for stimulus-locked 
% 	ERP (true) or not (false), default false.
% 	.respERPcor         = Boolean, read data corrected for response-locked 
% 	ERP (true) or not (false), default false.
% 	.bin.Type           = string, split up by what variable ('RT',
%   'GLM1vmPFCValenceMan'), optional.
%   .bin.Num            = numeric, number of bins (2 or 3), optional.
%   .bin.Type           = string, variable after which bins are created,
%   'RT' or 'BOLD', optional.
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

% Set root directory:
rootdir = grouplevel_set_rootdir(); % '/project/3017042.02';

% Add helper files to path:
fprintf('Add helper files to path\n');
addpath(fullfile(rootdir,'Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Helpers'));
addpath(fullfile(rootdir,'Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel'));

% Set directories and parameters:
dirs        = set_dirs(rootdir);
par         = set_par();

%% Initialize job:

job = []; % initialize empty job

job.nSub = 36; % necessary for validSubs
job.dirs = dirs; % add directories

% Data settings:
job.sub2exclude         = [];
% job.sub2exclude         = [11 12 15 23 25 26 30]; % exclusion in TAfT-outcome-locked

job.lockSettings        = 'stimlocked'; % stimlocked resplocked
job.baselineSettings    = 'trend'; % 'all' or 'condition' or 'trial' or 'ft_trial' or 'trend'    
job.baselineTimings     = [-0.250 -0.050];

job.responseSettings    = 'Go'; % 'Go' or 'Hand'  
job.actionSettings      = 'exeAct'; % 'reqAct' or 'exeAct'
job.accSettings         = 'correct'; % 'bothAcc' or 'correct' or 'incorrect' or 'noAcc' or 'reqAct'
% for bothAcc, always use reqAct

% Channel settings:
job.chanArea            = 'midfrontal'; % 'midfrontal' or 'frontal' or 'central' or 'Polania' or 'leftparietal' or 'leftmotor' or 'rightmotor'
% job.chanArea            = 'frontal'; % 'midfrontal' or 'frontal' or 'central' or 'Polania' or 'leftparietal' or 'leftmotor' or 'rightmotor'
% job.chanArea            = 'central'; % 'midfrontal' or 'frontal' or 'central' or 'Polania' or 'leftparietal' or 'leftmotor' or 'rightmotor'
% job.chanArea            = 'occipital'; % 'midfrontal' or 'frontal' or 'central' or 'Polania' or 'leftparietal' or 'leftmotor' or 'rightmotor'

% Contrast of interest:
job.contrastType        = 'Go'; % 'Congruency' or 'Go' or 'GoLeftRight' or 'GoValence' or 'GoAccuracy' or 'Accuracy' or 'IncongruentAccuracy' or 'CongruentAccuracy' or 'CorrectIncongruentIncorrectCongruent'

% ----------------------------------------------------------------------- %
% Split up in bins:
% job.bin.Type    = 'vmPFCValenceMan'; % RT vmPFCValenceMan
% job.bin.Num     = 3; % relevant only for RT

% ----------------------------------------------------------------------- %
% Load and prepare data:
job         = time_update_job(job); % initialize job

[job, data] = time_load_data(job); % load data
[job, data] = time_prepare_generic_data(job,data); % prepare generic objects

job         = time_update_job(job); % update job
[job, data] = time_prepare_contrast_data(job,data); % prepare objects specific to contrast

%% Permutation test:

% Initialize neighbors (needed only of not averaged over channels):
rng(70); % set seed for reproducibility of permutation test.

% Retrieve neighbours:
cfg_neighb          = [];
cfg_neighb.method   = 'distance';
cfg_neighb.elecfile = 'easycap-M1.txt';
cfg.neighbours      = ft_prepare_neighbours(cfg_neighb);
cfg.channel         = job.channels; 

% Average over channels, but not time bins:
% Channels:
cfg.avgoverchan     = 'no';
cfg.channel         = job.channels; % set by job.chanArea
% cfg.channel         = {'Fz','FCz','Cz'}; % midfrontal
% cfg.channel         = {'F1','F3','FCz','FC1','FC3','FC5','Cz','C1','C3'}; % vmPFC time domain
% Time:
cfg.latency         = [0.0 0.7];
cfg.avgovertime     = 'no';

% Other settings:
cfg.parameter       = 'avg';
cfg.method          = 'montecarlo';
cfg.statistic       = 'ft_statfun_depsamplesT';
cfg.alpha           = 0.05;
cfg.correctm        = 'cluster';
cfg.correcttail     = 'prob';
cfg.numrandomization = 1000;
cfg.design(1,1:2*job.nValidSubs)  = [ones(1,job.nValidSubs) 2*ones(1,job.nValidSubs)];
cfg.design(2,1:2*job.nValidSubs)  = [1:job.nValidSubs 1:job.nValidSubs];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number
fprintf('Perform permutation test over channels %s, time range %.03f - %.03f\n',...
    strjoin(string(cfg.channel),', '),cfg.latency(1),cfg.latency(end));
[stat]                  = ft_timelockstatistics(cfg,data.time1{job.validSubs},data.time2{job.validSubs});

%% Evaluate output:

evaluate_stat(stat,0.70)

%% T-test with Fieldtrip:

rng(70); % set seed for reproducibility of permutation test.

cfg             = [];
cfg.channel     = job.channels; % default 'all'
cfg.avgoverchan = 'yes';
cfg.latency     = [0.2 0.3];
cfg.avgovertime = 'yes';
cfg.parameter   = 'avg';
cfg.method      = 'analytic';
cfg.statistic   = 'ft_statfun_depsamplesT';
cfg.alpha       = 0.05;
cfg.correctm    = 'no';
cfg.design(1,1:2*job.nValidSubs)  = [ones(1,job.nValidSubs) 2*ones(1,job.nValidSubs)];
cfg.design(2,1:2*job.nValidSubs)  = [1:job.nValidSubs 1:job.nValidSubs];
cfg.ivar                = 1; % the 1st row in cfg.design contains the independent variable
cfg.uvar                = 2; % the 2nd row in cfg.design contains the subject number
stat = ft_timelockstatistics(cfg,data.time1{job.validSubs},data.time2{job.validSubs}); 
fprintf('T-test of %s vs. %s in range %.03f - %.03f seconds, channels %s:\nt(%d) = %.03f, p = %.03f\n',...
    job.twoLineLabels{1},job.twoLineLabels{2},cfg.latency(1),cfg.latency(end),strjoin(string(cfg.channel),'/'),...
    stat.df,stat.stat,stat.prob);

%% T-test in selected time range for selected channels:

% Select channels to average over:
selChans = {'Fz','FCz','Cz'};

% Select timerange to average over:
selTime = [0.2 0.3];

selTimeIdx = dsearchn(data.mu.time',selTime');
selTimeIdx = selTimeIdx(1):selTimeIdx(2);

% Extract mean signal in defined time/channel range for each subject:
difVec = nan(size(data.ERPdata,1),1); % Initialize empty vector
for iSub = 1:size(data.ERPdata,1)
    selChanIdx = find(ismember(data.ERPdata{iSub,1}.label, selChans)); % determine indices of channels for each subject
    % Define contrast per hand:
    % Average over time (2) and channels (1)
    % Go:
    difVec(iSub) = mean(mean(data.ERPdata{iSub,1}.avg(selChanIdx,selTimeIdx) + data.ERPdata{iSub,2}.avg(selChanIdx,selTimeIdx) - ...
        data.ERPdata{iSub,3}.avg(selChanIdx,selTimeIdx) - data.ERPdata{iSub,4}.avg(selChanIdx,selTimeIdx),2),1);
end

% Perform t-test:
[~,p,~,stats] = ttest(difVec); % only electrode [1 2] show the effect significantly.
fprintf('T-test of %s vs. %s in range %.03f - %.03f seconds, channels %s:\nt(%d) = %.03f, p = %.03f\n',...
    job.twoLineLabels{1},job.twoLineLabels{2},selTime(1),selTime(end),strjoin(string(selChans),'/'),...
    size(data.ERPdata,1)-1,stats.tstat,p);

% END
