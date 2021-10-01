% EEGfMRIPav_CueLocked_9_time_grouplevel_plot

% This is an interactive script---execute it step-by-step.
% Set the appropriate job settings for loading TF data aggregated across
% subjects, then create plots.
% - two-line plot contrasting two conditions of cotnrast
% - multi-line plot contrasting all conditions within data set
% - Topoplot contrasting two conditions of contrast
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
job.chanArea            = 'midfrontal'; % 'midfrontal' or 'frontopolar' or 'frontal' or 'leftfrontal' or 'midcentral' or 'central' or 'leftmotor' or 'rightmotor' or ...

% Contrast of interest:
job.contrastType        = 'Go'; % 'Congruency' or 'Go' or 'GoLeftRight' or 'GoValence' or 'GoAccuracy' or 'Accuracy' or 'IncongruentAccuracy' or 'CongruentAccuracy' or 'CorrectIncongruentIncorrectCongruent'

% ----------------------------------------------------------------------- %
% Split up in bins:
% job.bin.Type    = 'GLM5FvmPFCValenceMan'; % RT GLM5FvmPFCValenceMan
% job.bin.Num     = 3; % only relevant for RT

% ----------------------------------------------------------------------- %
% Load and prepare data:
job         = time_update_job(job); % initialize job

[job, data] = time_load_data(job); % load data
[job, data] = time_prepare_generic_data(job,data); % prepare generic objects

job         = time_update_job(job); % update job
[job, data] = time_prepare_contrast_data(job,data); % prepare objects specific to contrast

%% 1) TWOLINEPLOT ON MAIN CONTRAST:

% job = rmfield(job,'sigTime'); % remove sigTime if inappropriate
twoLinePlot(job,data)

% set(gca,'xtick',0:0.05:1,'xtickLabel',0:0.05:1,'fontsize',12) % finer time scale

saveas(gcf,fullfile(dirs.timeplot,sprintf('twoLinePlot_%s.png',job.plotName)));
pause(3)
close gcf

%% 2) MULTILINEPLOT FOR ALL CUE CONDITIONS:

multiLinePlot(job,data)

% set(gca,'ylim',[-0.025 0.035]); % for 'all'
% set(gca,'xtick',0:0.05:1,'xtickLabel',0:0.05:1,'fontsize',12) % finer time scale

saveas(gcf,fullfile(dirs.timeplot,sprintf('multiLineplot_%s.png',job.plotName)))
close gcf

%% 3A) SINGLE TOPO PLOTS.

% zlim:
zlim = 0.015; % task effects
% zlim = 0.010; % BOLD effects
% zlim = 0.005; % BOLD effects even weaker

% Timing: If no sigTime given in job, adjust manually
job.sigTime = [0.341 0.412]; % N2/P3

singleTopoPlot(job,data,zlim)

% Save:
saveas(gcf,fullfile(dirs.timeplot,sprintf('Topoplot_%s_TOI_%03d_%03dms_zlim%.03f.png',...
    job.plotName, job.sigTime(1)*1000, job.sigTime(end)*1000, zlim)))
pause(3)
close gcf

%% 3B) MULTIPLE TOPO PLOTS OVER TIME.

zlim = 0.015; % task effects
% zlim = 0.010; % BOLD effects
% zlim = 0.005; % BOLD effects (stronger)

% Settings for grid:
% Stim-locked:
nRows = 2; startTime = 0; endTime = 0.45; steps = 0.05; % 2 rows (50 ms)
% nRows = 2; startTime = 0; endTime = 0.90; steps = 0.10; % 2 rows (100 ms)
% Resplocked:
% nRows = 2; startTime = -0.20; endTime = 0.25; steps = 0.05; % 2 rows (50 ms)
% nRows = 2; startTime = -0.40; endTime = 0.55; steps = 0.10; % 2 rows (100 ms)

% Plot
multiTopoPlot(job,data,zlim,nRows,startTime,endTime,steps) % 

% Save:
saveas(gcf,fullfile(dirs.timeplot,sprintf('multiTopoPlot_%s_TOISequence_%d-%dms_zlim%.03f.png',...
    job.plotName,1000*startTime,1000*(endTime+steps),zlim)))
pause(3)
close gcf

% END.
