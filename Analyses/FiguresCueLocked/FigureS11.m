% FigureS11.m

% Plots for Figure S11 in supplementary material.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/FiguresCueLocked

% Set root directory:
rootdir = figures_set_rootdir(); % '/project/3017042.02';

addpath(fullfile(rootdir,'Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Helpers'));
addpath(fullfile(rootdir,'Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel'));

dirs        = set_dirs(rootdir);
par         = set_par();

%% Load data:

% Figure S11A A: accSettings = 'correct'; lockSettings = 'stimlocked', band = 'theta', contrastType = 'Go';
% Figure S11A B: accSettings = 'correct'; lockSettings = 'resplocked', band = 'theta', contrastType = 'Go';
% Figure S11A C: lockSettings = 'resplocked', band = 'theta', contrastType = 'Go';
% load both accSettings = 'correct' and accSettings = 'incorrect'; concatenate

% Figure S11B A: accSettings = 'correct'; lockSettings = 'stimlocked', band = 'theta', contrastType = 'Go';
% Figure S11B B: accSettings = 'correct'; lockSettings = 'stimlocked', band = 'alpha', contrastType = 'Go';
% Figure S11B C: accSettings = 'correct'; lockSettings = 'stimlocked', band = 'beta', contrastType = 'Go';

job = []; % initialize empty job

job.dirs = dirs; % add directories

% Data settings:
job.nSub                = 36; % necessary for validSubs before loading data
job.sub2exclude         = []; % include all subjects
% job.sub2exclude         = [11 12 15 23 25 26 30]; % exclusion in TAfT-outcome-locked

job.lockSettings        = 'stimlocked'; % stimlocked resplocked

job.baselineSettings    = 'trend'; % 'all' or 'condition' or 'trial' or 'ft_trial' or 'trend'    
job.baselineTimings     = [-0.250 -0.050];

job.responseSettings    = 'Go'; % 'Go' or 'Hand'  
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
% Bin type:
job.bin.Type    = 'RT'; % RT
job.bin.Num     = 3;

% Load and prepare data:
job         = TF_update_job(job); % initialize job
[job, data] = TF_load_data(job); % load data
[job, data] = TF_prepare_generic_data(job,data); % prepare generic objects
job         = TF_update_job(job); % update job
[job, data] = TF_prepare_contrast_data(job,data); % prepare objects specific to contrast

%% Figure S11. EEG MULTILINEPLOT

% Baselines:
if strcmp(job.lockSettings,'stimlocked')
    iBaseline = find(data.mu.time == 0);  % iBaseline = 41;
    baseline = min(squeeze(nanmean(data.SubCondTime(:,:,iBaseline))));
    
elseif strcmp(job.lockSettings,'resplocked')
    iBaseline = 1;
    baseline = min(squeeze(nanmean(data.SubCondTime(:,:,iBaseline))));
end

% Start figure:
figure('Position',[100 100 1200 800]); hold on
clear p

% A) With error bars:
lineWidth = 5; until = 100; transp = 0.10;
for iCond = 1:job.nCond
    p{iCond} = boundedline(data.mu.time(1:until),...
        squeeze(nanmean(data.SubCondTime(:,iCond,1:until)))-baseline,...
        job.nCond/(job.nCond-1)*squeeze(nanstd(squeeze(data.SubCondTime(job.validSubs,iCond,1:until))-data.SubTime(job.validSubs,1:until)+repmat(data.GrandTime(1:until),job.nValidSubs,1)))'./sqrt(job.nValidSubs),...
        'cmap',job.colMat(iCond,:),'alpha','transparency',transp); % ,'linewidth',2); hold on
    set(p{iCond},'linestyle',job.lineStyle{iCond})
    set(p{iCond},'Linewidth',lineWidth)
end

% General axis settings:
ylabel('Power (dB)','fontweight','bold','fontsize',32); 
xlabel('Time (s)','fontweight','bold','fontsize',32)
yLim = get(gca,'ylim'); % yLim = [0 2.5];

if strcmp(job.contrastType,'Congruency')
    yMinLim = -0.75;
    yMaxLim = 0.95;

elseif strcmp(job.contrastType,'Accuracy')
    yMinLim = -0.25;
    yMaxLim = 2.3;

elseif strcmp(job.contrastType,'Go') && strcmp(job.lockSettings,'stimlocked')
    yMinLim = -0.25; % -0.25 for only theta; -1 to compare all
    yMaxLim = 2.1;

elseif strcmp(job.contrastType,'Go') && strcmp(job.lockSettings,'resplocked')
    yMinLim = -0.25;
    yMaxLim = 2.5; % or 2.3

elseif strcmp(job.contrastType,'GoLeftRight') && strcmp(job.lockSettings,'stimlocked') && strcmp(job.band,'theta')
    yMinLim = -0.4;
    yMaxLim = 2.3;

elseif strcmp(job.contrastType,'GoLeftRight') && strcmp(job.lockSettings,'stimlocked') && strcmp(job.band,'beta')
    yMinLim = -0.8;
    yMaxLim = 1.7;

elseif strcmp(job.contrastType,'GoLeftRight') && strcmp(job.lockSettings,'resplocked') && strcmp(job.band,'theta')
    yMinLim = -0.4;
    yMaxLim = 2.75;

elseif strcmp(job.contrastType,'GoLeftRight') && strcmp(job.lockSettings,'resplocked') && strcmp(job.band,'beta')
    yMinLim = -0.8;
    yMaxLim = 1.7;

else
    error('Unknown contrast')
end

% Add x-axis labels, vertical lines:
if strcmp(job.lockSettings,'stimlocked')
    set(gca,'xlim',[-0.25 1.3],'xtick',[0 0.5 0.815 1 1.3],'xtickLabel',{'Cue','0.5','AvgRT','1','+'},...
        'ylim',[yMinLim yMaxLim],'ytick',-0.5:0.5:3,'fontsize',32,'Linewidth',3) %,'ytick',-5:0.1:4)
    plot([0 0],get(gca,'ylim'),':k','LineWidth',3); % cue onset
    plot([.815 .815],get(gca,'ylim'),':k','LineWidth',3); % average RT 

elseif strcmp(job.lockSettings,'resplocked')
    set(gca,'xlim',[-1 1],'xtick',-1:0.5:1,'xtickLabel',{'-1','-0.5','Resp','0.5','1'},...
        'ylim',[yMinLim yMaxLim],'ytick',-0.5:0.5:3,'fontsize',32,'Linewidth',3) %,'ytick',-5:0.1:4)
    plot([0 0],get(gca,'ylim'),':k','LineWidth',3);

else
    error('invalid locking setting')
end

% Save:
if job.stimERPcor
    saveas(gcf,sprintf('/project/3017042.02/Log/CueLockedPaperPlots/Multilineplot_%s_%s_induced.png',job.contrastType,job.lockSettings))
elseif ~isempty(job.invalidSubs)
    saveas(gcf,sprintf('/project/3017042.02/Log/CueLockedPaperPlots/multiLineplot_%s_without6.png',job.multiPlotName));
else
    saveas(gcf,sprintf('/project/3017042.02/Log/CueLockedPaperPlots/multiLineplot_%s.png',job.multiPlotName));
end
close gcf

% END.
