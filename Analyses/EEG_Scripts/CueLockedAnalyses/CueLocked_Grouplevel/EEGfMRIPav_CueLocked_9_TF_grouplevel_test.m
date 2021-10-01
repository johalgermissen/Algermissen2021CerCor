% EEGfMRIPav_CueLocked_9_TF_grouplevel_test

% This is an interactive script---execute it step-by-step.
% Set the appropriate job settings for loading TF data aggregated across
% subjects, then perform tests.
% - cluster-based permutation test with Fieldtrip
% - evaluate stat output object
% - plot stat output object as topoplot or TF plot
% - save clusters above threshold as mask for trial-by-trial extraction.
% - plot mask.
% - perform t-test per electrode and averaged over electrodes.
% - export mean of time/frequency/channel range per subject per condition.
% - bar-plot mean of time/frequency/channel range per subject per condition.
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
job.accSettings         = 'bothAcc'; % 'bothAcc' or 'correct' or 'incorrect' or 'noAcc' or 'reqAct'
% for bothAcc, always use reqAct

job.TFtype              = 'hanning'; % hanning morlet
job.nFreq               = 15; % 64 15 6

% Channel settings:
job.chanArea            = 'midfrontal'; % 'midfrontal' or 'frontal' or 'central' or 'leftmotor' or 'rightmotor'

% Band settings:
job.band                = 'alpha'; % 'theta' or 'alpha' or 'beta' or 'broad' 

% Contrast of interest:
job.contrastType        = 'Congruency'; % 'Congruency' or 'Go' or 'GoLeftRight' or 'GoValence' or 'GoAccuracy' or 'Accuracy' or 'CongruentAccuracy' or 'IncongruentAccuracy'

% ----------------------------------------------------------------------- %
% ERP-corrected:
job.stimERPcor = false;
job.respERPcor = false;

% ----------------------------------------------------------------------- %
% Split up in bins:
% job.bin.Type    = 'RT'; % RT GLM5FvmPFCValenceMan
% job.bin.Num     = 3; % only relevant for RT

% Load and prepare data:

job         = TF_update_job(job); % initialize job
[job, data] = TF_load_data(job); % load data
[job, data] = TF_prepare_generic_data(job,data); % prepare generic objects

job         = TF_update_job(job); % update job
[job, data] = TF_prepare_contrast_data(job,data); % prepare objects specific to contrast

%% PERMUTATION TEST:

%% 1) Averaging over channels and/ or frequencies:

% Option A: Average over nothing, retrieve both separate channels and 
% frequencies that are part of cluster above threshold
cfg                 = [];
cfg.avgoverchan = 'no'; % 
cfg.avgoverfreq = 'no'; % 
fprintf('Average neither over channels nor frequencies\n');

% Option B: average over channels, only retrieve separate frequencies that 
% are part of cluster above threshold (--> topoplot)
cfg                 = [];
cfg.avgoverchan     = 'no'; % 
cfg.avgoverfreq     = 'yes'; %
% fprintf('Average over frequencies, but not channels\n');

% Option C: average over frequencies, only retrieve separate channels that 
% are part of cluster above threshold (--> TF plot)
cfg                 = [];
cfg.avgoverchan     = 'yes'; % 
cfg.avgoverfreq     = 'no'; %
% fprintf('Average over channels, but not frequencies\n');

% Option D: average over both channels and frequencies, only obtain p-value
cfg                 = [];
cfg.avgoverchan = 'yes'; % 
cfg.avgoverfreq = 'yes'; % 
% fprintf('Average over both channels and frequencies\n');

%% 2) Frequencies band to select:

Band = 'theta'; % delta theta alpha beta broad
fprintf('Selected frequency band is %s\n',Band);

%% 3) Run permutation test:

rng(70); % set seed for reproducibility of permutation test.

% Initialize neighbours:
cfg_neighb          = [];
cfg_neighb.method   = 'distance';
cfg_neighb.elecfile = 'easycap-M1.txt';
cfg.neighbours      = ft_prepare_neighbours(cfg_neighb);

% Select channels:
cfg.channel         = job.channels; %  {'Fz','FCz','Cz'};
if length(cfg.channel) == 1 % if single channel, must average:
    cfg.avgoverchan = 'yes';
end

% Select frequency band:
if strcmp(Band,'broad')
    cfg.frequency       = [0.5 15];
elseif strcmp(Band,'delta')
    cfg.frequency       = [1.5 4];
elseif strcmp(Band,'theta')
    cfg.frequency       = [4 8];
elseif strcmp(Band,'alpha')
    cfg.frequency       = [8 13];
elseif strcmp(Band,'beta')
    cfg.frequency       = [13 15];
end

% Time:
if strcmp(job.lockSettings,'stimlocked')
    cfg.latency         = [0 1.3]; % 0 1.3
elseif strcmp(job.lockSettings,'resplocked')
    cfg.latency         = [-1 0.5];
end

% Permutation test settings:
cfg.method          = 'montecarlo';
cfg.statistic       = 'ft_statfun_depsamplesT';
cfg.correctm        = 'cluster';
cfg.clusteralpha    = 0.05; % threshold at t > 2
% cfg.clusteralpha    = 0.001; % threshold at t > 3.1
cfg.clusterstatistic= 'maxsum'; % maxsum, maxsize
cfg.minnbchan       = 1;
cfg.tail            = 0; % 0 for two-sided, 1 for one-sided positive, -1 for one-sided negative; see http://www.fieldtriptoolbox.org/reference/ft_statistics_montecarlo/
cfg.clustertail     = cfg.tail; % must correspond to cfg.tail
% cfg.correcttail = 'alpha';
cfg.alpha           = 0.05; % cut-off for thresholding t-values
cfg.numrandomization = 500; % 500 by default, can be set higher
subj                = length(job.validSubs); % number of subjects
design              = zeros(2,2*subj); % initialize design matrix
design(1,:)         = repmat(1:subj,[1,2]); % subject numbers twice (congruent/ incongruent)
design(2,:)         = [ones(1,subj) 2*ones(1,subj)]; % first 1 for each subject (congruent), then 2 for each subject (incongruent)
cfg.design          = design;
cfg.uvar            = 1; % unit variables: (1 for) within 1 subject, not (2 for) across subjects
cfg.ivar            = 2; % independent variables: 2 conditions
fprintf('>>> Perform permutation test for channels %s, %.1f - %.1f Hz, %.1f - %.1f sec\n',strjoin(job.channels,'/'),cfg.frequency(1),cfg.frequency(end),cfg.latency(1),cfg.latency(end))
% General purpose statement:
[stat]          = ft_freqstatistics(cfg,data.TF1{job.validSubs},data.TF2{job.validSubs});

%% 4) Evaluate stat:

evaluate_stat(stat,0.05)

%% 5) Plot STATS output as TF plot with Fieldtrip (averaged across channels):
% Needs stat being averaged over channels, but not (!) over frequencies
% Will ignore limited frequency range; always broadband.

zlim = 2;

% Settings for topoplot:
cfg = [];
if strcmp(job.lockSettings,'stimlocked')
    cfg.xlim = [-0 1.3]; % -0.25 1.3
elseif strcmp(job.lockSettings,'resplocked')
    cfg.xlim = [-1 1]; % -0.25 1.3
else
    error('Invalid response settings')
end

% Start plot:
cfg.ylim = [1 15];
cfg.zlim = [zlim*-1 zlim*1];
cfg.maskparameter = 'mask';
cfg.maskalpha = 0.05; % 0.025;
cfg.parameter = 'stat';
ft_singleplotTFR(cfg, stat);
% Save:
saveas(gcf,fullfile(dirs.TFplot,sprintf('TFplot_PermutationCluster_%s_%s_zlim%d.png',...
    job.plotName,strjoin(sort(job.channels),''),zlim)));
pause(3)
close gcf

%% 6) Plot STATS output as topoplot with Fieldtrip (not averaged across channels):

% Topoplot with clusterplot (needs non-averaged channels):
% needs time or frequency be a singleton dimension (but not channels!).

zlim = 6;

cfg             = [];
cfg.alpha       = 0.05;
cfg.parameter   = 'stat';
cfg.zlim        = [-1*zlim 1*zlim];
cfg.layout      = job.layout; 
ft_clusterplot(cfg, stat);
saveas(gcf,fullfile(dirs.TFplot,sprintf('Topoplot_permutation_%s_Band_%s_zlim%d.png',...
    job.plotName,Band,zlim)));
pause(3)
close gcf

%% SAVE MASK:

% In preparation for saving mask:
% Run permutation test with cfg.avgoverchan = 'no'; cfg.avgoverfreq = 'no'; % 
% one-sided (only get all-negative or all-positive mask);
% band-limited (only focus on one cluster)

% Outputs from permutation test used here:
% stat.mask % 1 if within significant cluster, 0 outside
% stat.posclusterslabelmat % > 0 if part of respective cluster number, 0 outside
% stat.negclusterslabelmat % > 0 if part of respective cluster number, 0 outside
% stat.stat % t-values unthresholded

% Threshold t-values based on significance:
statMask    = stat.mask .* stat.stat; % mask [significance: 0/1] * t-values (continuous)

% Detect channel, time and frequency indices of permutation test settings (not result!) performed:
maskChanIdx = find(ismember(data.mu.label, cfg.channel)); % retrieve indices of channels used
maskTimeIdx = dsearchn(data.mu.time',stat.time'); % retrieve indices of timing of stat
maskFreqIdx = dsearchn(data.mu.freq',stat.freq'); % retrieve indices of frequencies of stat

% Create empty mask, insert t-values:
freqMask = zeros(size(data.mu.powspctrm)); % create zeros as default
freqMask(maskChanIdx,maskFreqIdx,maskTimeIdx) = statMask; % add statMask to it

% Impose additional restriction based on time:
% selTimeIdx = dsearchn(data.mu.time',0.5);
% freqMask(:,:,1:selTimeIdx) = 0; % till selTimeIdx 
% freqMask(:,:,selTimeIdx:end) = 0; % from selTimeIdx onwards

% If negative: flip
if sum(freqMask(:)) < 0; freqMask = freqMask*-1; end

% Determine channel name:
if numel(stat.label) > 1
    chanName = strjoin(sort(cfg.channel),'');
else
    chanName = cfg.channel;
end
fprintf('Create mask for contrast %s, channels %s, frequency band %s\n', ...
    job.contrastType, chanName, Band);

% ----------------------------------------------------------------------- %
% Plot TF power (averaged over channels):
figure('Position',[100 100 1000 600]); hold on
contourf(data.mu.time,data.mu.freq,squeeze(nanmean(freqMask,1)),40,'linestyle','none')
% Settings:
fontSize    = 16;
lineWidth   = 3;
set(gca,'xlim',[data.mu.time(1) data.mu.time(end)],'ylim',[data.mu.freq(1) data.mu.freq(end)],...
    'xtick',-3:0.5:3,... 
    'ytick',2:2:14,'yscale','lin',...
    'fontsize',fontSize,'Linewidth',lineWidth) %
% Labels:
ylabel('Frequency (Hz)','fontsize',fontSize);
xlabel('Time (s)','fontsize',fontSize);
% Title:
title(sprintf('Mask: %s contrast, channels %s, %s band, %s',...
    job.contrastType, chanName, Band, job.lockSettings),'fontsize',fontSize);

% ----------------------------------------------------------------------- %
% Save mask:
fileName = fullfile(dirs.mask,sprintf('Mask_%s_%s_%s_%s.mat',...
    job.contrastType, chanName, Band, job.lockSettings));
save(fileName,'freqMask'); % use 'Band' as used for permutation test

%% Load and plot mask, save plot:

% Specify settings:
contrast    = 'Go'; % Go Congruency
channels    = 'CzFCzFz'; % CzFCzFz Cz
band        = 'theta'; % theta alpha beta broad
lockSetting = 'resplocked'; % stimlocked resplocked

% Load:
load(fullfile(dirs.mask,sprintf('Mask_%s_%s_%s_%s.mat', contrast, channels, band, lockSetting)));

% Figure:
figure('Position',[100 100 800 800]); 
contourf(data.mu.time,data.mu.freq,squeeze(nanmean(freqMask,1)),40,'linestyle','none');
% Settings:
fontSize    = 16;
lineWidth   = 3;
set(gca,'xlim',[data.mu.time(1) data.mu.time(end)],'ylim',[data.mu.freq(1) data.mu.freq(end)],...
    'xtick',-3:0.5:3,... 
    'ytick',2:2:14,'yscale','lin',...
    'fontsize',fontSize,'Linewidth',lineWidth) %
% Labels:
ylabel('Frequency (Hz)','fontsize',fontSize);
xlabel('Time (s)','fontsize',fontSize);
% Title:
title(sprintf('Mask: %s contrast, channels %s, %s band, %s',...
    contrast, channels, band, lockSetting),'fontsize',fontSize);
% Save:
saveas(gcf,fullfile(dirs.mask,sprintf('Mask_%s_%s_%s_%s.png',...
    contrast, channels, band, lockSetting)));
pause(3)
close gcf

%% T-TEST separately per electrode & averaged over electrodes:

% Which electrodes contribute to this effect? Post-hoc alpha = .05/3 = .017

% Select time:
selTime = job.sigTime;
% selTime = [0.525 1.300];
selTimeIdx = dsearchn(data.mu.time',selTime'); % indices of start and end time
selTimeIdx = selTimeIdx(1):selTimeIdx(end); % all indices in between

% Select frequencies:
selFreq = job.freq';
% selFreq = [13 15]';
selFreqIdx = dsearchn(data.mu.freq',selFreq); % indices of lowest/highest frequency in that range
selFreqIdx = selFreqIdx(1):selFreqIdx(end); % all indices in between

% Compute differences per subject per channel in selected time and frequency range:
fprintf('Extract mean TF power for %.03f - %0.3f sec., %d-%d Hz:\n',selTime(1),selTime(end),selFreq(1),selFreq(end))
difData = nan(job.nSub,length(job.channels));
for iSub = 1:job.nSub % iSub = 1;
    selChanIdx = find(ismember(data.TF1{iSub}.label, job.channels)); % determine indices of channels per subject
    difData(iSub,:) = mean(mean(data.TF1{iSub}.powspctrm(selChanIdx,selFreqIdx,selTimeIdx),3),2) - ...
        mean(mean(data.TF2{iSub}.powspctrm(selChanIdx,selFreqIdx,selTimeIdx),3),2);
end

% Loop over channels, compute t-test:
fprintf('t-test per channel:\n')
for iChan = 1:length(job.channels)
    difVec         = difData(:,iChan);
    [~,P,~,STATS]  = ttest(difVec);
    fprintf('Channel %s: t(%d) = %.02f, p = %.03f, d = %.02f\n',job.channels{iChan},STATS.df,STATS.tstat,P/2,mean(difVec)/std(difVec)); % divide  by 2 because one-sided
end

% ----------------------------------------------------- %
% Overall t-test:
difVec          = mean(difData,2); % average over electrodes
[H,P,CI,STATS]  = ttest(difVec);
fprintf('Overall: t(%d) = %.02f, p = %.03f, d = %.02f\n',STATS.df,STATS.tstat,P/2,mean(difVec)/std(difVec)); % divide  by 2 because one-sided

%% EXPORT DATA in selected time/frequency/channel range per subject per condition:

% Select time:
selTime = job.sigTime;
selTimeIdx = dsearchn(data.mu.time',selTime'); % indices of start and end time
selTimeIdx = selTimeIdx(1):selTimeIdx(end); % all indices in between

% Select frequency:
selFreq = job.freq;
selFreqIdx = dsearchn(data.mu.freq',selFreq'); % indices of lowest/highest frequency in that range
selFreqIdx = selFreqIdx(1):selFreqIdx(end); % all indices in between

% Select channels:
selChan = job.channels;

% Extract mean data in selected time/frequency/channel range per subject:
% Compute differences per subject per channel in selected time and frequency range:
fprintf('Extract mean TF power for %.03f - %0.3f sec., %d-%d Hz, %s:\n',...
    selTime(1),selTime(end),selFreq(1),selFreq(end),strjoin(selChan,'/'));
selData = nan(job.nSub,job.nCond);
for iSub = 1:job.nSub % iSub = 1;
    selChanIdx = find(ismember(data.TFall{iSub}.ValxAct{iCond}.label, selChan)); % determine indices of channels per subject
    for iCond = 1:job.nCond % iCond = 1;
        selData(iSub,iCond) = nanmean(nanmean(nanmean(data.TFall{iSub}.ValxAct{iCond}.powspctrm(selChanIdx,selFreqIdx,selTimeIdx))));
    end
end

% Save as csv:
fileName = fullfile(dirs.TFgroup,sprintf('TFSubCond_%s.csv',job.plotName));
fprintf('Save data under %s\n',fileName);
csvwrite(fileName,selData);

%% BAR PLOT data in selected time/frequency/channel range per subject per condition:

% x-labels need to be adjusted based on particular contrast and
% accSettings

points          = true; % plot individual data points

% Determine timing:
selTime         = job.sigTime; % retrieving timing of ROI
% selTime = [0.1750 0.3250];
selTimeIdx      = dsearchn(data.mu.time',selTime'); % retrieve edge indices in matrix
selTimeIdx      = selTimeIdx(1):selTimeIdx(end); % insert intermediate indices

% Baseline:
baseTimeIdx     = dsearchn(data.mu.time',0); % retrieve index at t = 0
baseline        = -1*nanmean(nanmean(nanmean(data.SubCondTime(:,:,baseTimeIdx)))); % retrieve baseline data

% Baseline correction:
subCondMean     = baseline + squeeze(nanmean(data.SubCondTime(:,:,selTimeIdx),3));
job.nValidSubs       = size(subCondMean,1);

% Cousineau-Morey SE: average over time, subtract subject overall mean, add sample overall mean
nCond = job.nCond;
subMean = nanmean(subCondMean,2);
grandMean = nanmean(nanmean(nanmean(data.SubCondTime)));
condMean = nan(job.nCond,1); condSE = nan(job.nCond,1);
for iCond = 1:nCond
    condMean(iCond) = nanmean(subCondMean(:,iCond));
    condSE(iCond) = nCond/(nCond-1)*nanstd(subCondMean(:,iCond)-subMean+repmat(grandMean,job.nValidSubs,1))./sqrt(job.nValidSubs);
end

% Settings:
colMat      = repmat([0 .6 .2; .8 0 0],job.nCond/2,1);
xLoc        = [1 2 3 4 5.5 6.5 7.5 8.5];
yLim        = [-0.6 0.6]; yTick = -1:.2:1;

lineWidth   = 3;
fontSize    = 24;

figure('Position',[100 100 1200 800]); hold on
p   = cell(job.nCond,1);
for iCond = 1:job.nCond
    p{iCond} = bar(xLoc(iCond),condMean(iCond),.75,'FaceColor',colMat(iCond,:));
    errorbar(xLoc(iCond),condMean(iCond),condSE(iCond),'k','linestyle','none','Linewidth',lineWidth);
    if points
        s = scatter(repmat(xLoc(iCond),1,job.nValidSubs),subCondMean(:,iCond)',[],'k', 'jitter','on', 'jitterAmount',0.15); hold on % was 0.05
        set(s,'MarkerEdgeColor',[0.4 0.4 0.4],'linewidth',3); % was 1 
        yLim = [round(min(subCondMean(:)),2)-.01 round(max(subCondMean(:)),2)+.01];
        yTick = -10:1:10;
    end
end
ylabel(sprintf('%s Power',job.bandName),'FontSize',fontSize);
xlabel('Performed action x Accuracy','FontSize',fontSize); % adjust based on settings
set(gca,'xlim',[0 10],'xtick',[1.5 2.5 3.5 6 7 8],...
    'xticklabel',{'Go','\newline Correct','NoGo','Go','\newline Incorrect','NoGo'},...
    'ylim',yLim,'ytick',yTick,...
    'FontSize',fontSize) % adjust based on settings
legend([p{1},p{2}],{'Win','Avoid'}); legend boxoff % adjust based on settings
% title(sprintf('Valence x Executed Action x Accuracy'),'FontSize',fontSize);

% Save:
saveas(gcf,fullfile(dirs.TFplot,'Barplot_%s_%s.png',job.plotName));
pause(3)
close gcf

% END.
