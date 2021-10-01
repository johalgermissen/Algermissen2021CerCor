% FigureS13.m

% Plots for Figure S13 in supplementary material.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/FiguresCueLocked

% Set root directory:
rootdir     = figures_set_rootdir(); % '/project/3017042.02';

addpath(fullfile(rootdir,'Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Helpers'));
addpath(fullfile(rootdir,'Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel'));

dirs        = set_dirs(rootdir);
par         = set_par();

%% Load data:

% Figure S13 A: job.lockSettings = 'stimlocked';
% Figure S13 B: job.lockSettings = 'resplocked';

job = []; % initialize empty job

job.dirs = dirs; % add directories

% Data settings:
job.nSub                = 36; % necessary for validSubs before loading data
job.sub2exclude         = []; % include all subjects
% job.sub2exclude         = [11 12 15 23 25 26 30]; % exclusion in TAfT-outcome-locked
job.lockSettings        = 'stimlocked'; % stimlocked resplocked

job.baselineSettings    = 'trend'; % 'all' or 'condition' or 'trial' or 'ft_trial' or 'trend'    
job.baselineTimings     = [-0.250 -0.050];

job.responseSettings    = 'Hand'; % 'Go' or 'Hand'  
job.actionSettings      = 'exeAct'; % 'reqAct' or 'exeAct'
job.accSettings         = 'correct'; % 'bothAcc' or 'correct' or 'incorrect' or 'noAcc' or 'reqAct'
% for bothAcc, always use reqAct
job.TFtype              = 'hanning'; % hanning morlet
job.nFreq               = 15; % 15
% Channel settings:
job.chanArea            = 'midfrontal'; % 'midfrontal' or 'frontal' or 'central' or 'leftmotor' or 'rightmotor'
% Band settings:
job.band                = 'theta'; % 'theta' or 'alpha' or 'beta' or 'broad' 
% Contrast of interest:
job.contrastType        = 'Go'; % 'Congruency' or 'Go' or 'GoLeftRight' or 'GoValence' or 'GoAccuracy' or 'Accuracy' or 'CongruentAccuracy' or 'IncongruentAccuracy'
% ERP-corrected:
job.stimERPcor = false;
job.respERPcor = false;

% Load and prepare data:
job         = TF_update_job(job); % initialize job
[job, data] = TF_load_data(job); % load data
[job, data] = TF_prepare_generic_data(job,data); % prepare generic objects
job         = TF_update_job(job); % update job
[job, data] = TF_prepare_contrast_data(job,data); % prepare objects specific to contrast

%% TF plot relative to baseline:

fontSize = 24;

if strcmp(job.lockSettings,'stimlocked')
    baseTime    = 0;
    clim        = [-2 2];
elseif strcmp(job.lockSettings,'resplocked')
    baseTime    = -0.5;
    clim        = [-2 2];
else
    error('Unknown lock setting')
end

% Start plot:
% figure('units','normalized','outerposition',[0 0 1 1]); hold on % fullscreen
figure('Position',[100 100 1200 800]); 

for iPlot = 1:job.nCond
    subplot(2,job.nCond/2,job.condOrder(iPlot)); hold on
    job.chanIdx = find(ismember(data.Cond{iPlot}.label, job.channels)); % determine indices of channels
    dat2plot    = squeeze(mean(data.Cond{iPlot}.powspctrm(job.chanIdx,:,:))); % retrieve relevant data
    baseIdx     = find(round(data.Cond{iPlot}.time,3)==baseTime); % retrieve index of baseline
    baseline    = squeeze(mean(data.Cond{iPlot}.powspctrm(job.chanIdx,:,baseIdx))); % retrieve baseline vector
    dat2plot    = dat2plot - repmat(baseline',1,length(data.Cond{iPlot}.time)); % subtract baseline
    contourf(data.mu.time,data.mu.freq,dat2plot,40,'linestyle','none') % plot data
    % Axes settings:
    set(gca,'ytick',2:2:14,'yscale','lin','xtick',[-1 -0.5 0 .5 1],'TickLength',[0.01 0.1],... % yscale log
        'ylim',[1 14.9],'fontsize',28,'clim',[clim(1) clim(end)],'xlim',[job.TFtiming],'fontsize',fontSize) % -.25 1.3
    box off
    % Extra lines:
    plot([0 0], get(gca,'ylim'),'k')
    if strcmp(job.lockSettings,'stimlocked')
        plot([.815 .815], get(gca,'ylim'),':k')
    end
    % Title:
    title(job.condNames{iPlot});
    % Labels:
    if mod(job.condOrder(iPlot),job.nCond/2)==1
        xlabel('Time','fontsize',fontSize);
        ylabel('Frequency (Hz)','fontsize',fontSize);
%       Add legend:
%     elseif iPlot == job.nCond
%         colorbar
    end
end

% Save:
saveas(gcf,fullfile(rootdir,sprintf('Log/CueLockedPaperPlots/TFPlot_allConditions_%s.png',job.lockSettings)));
pause(3)
close all

% END.