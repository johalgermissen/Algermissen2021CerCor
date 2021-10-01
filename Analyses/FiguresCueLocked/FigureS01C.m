% FigureS01C.m

% Plots for Figure S01C in supplementary material.
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

% Figure 01C A, B: lockSettings = 'resplocked', band = 'theta', contrastType = 'Go';
% Figure 01C C, D: lockSettings = 'stimlocked', band = 'theta', contrastType = 'Go';
% Figure 01C E, F: lockSettings = 'stimlocked', band = 'alpha', contrastType = 'Conflict';

job = []; % initialize empty job

job.dirs = dirs; % add directories

% Data settings:
job.nSub                = 36; % necessary for validSubs before loading data
% job.sub2exclude         = []; % include all subjects
job.sub2exclude         = [11 12 15 23 25 26 30]; % exclusion in TAfT-outcome-locked

job.lockSettings        = 'resplocked'; % stimlocked resplocked

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

% Load and prepare data:
job         = TF_update_job(job); % initialize job
[job, data] = TF_load_data(job); % load data
[job, data] = TF_prepare_generic_data(job,data); % prepare generic objects
job         = TF_update_job(job); % update job
[job, data] = TF_prepare_contrast_data(job,data); % prepare objects specific to contrast

%% Figure S01C. EEG MULTILINEPLOT

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
        job.nCond/(job.nCond-1)*squeeze(nanstd(squeeze(data.SubCondTime(job.validSubs,iCond,1:until))-data.SubTime(job.validSubs,1:until)+repmat(data.GrandTime(1:until),length(job.validSubs),1)))'./sqrt(length(job.validSubs)),...
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

%% Figure S01C: EEG TF-Plots

% Settings:
addpath(fullfile(rootdir,'/Analyses/colorMaps'));
set(0, 'DefaultFigureColormap', jet(64)); cmap = 'jet';  
% set(0, 'DefaultFigureColormap', parula(64)); cmap = 'parula';  
% set(0, 'DefaultFigureColormap', redblue(64)); cmap = 'redblue'; 
% set(0, 'DefaultFigureColormap', turbo); cmap = 'turbo'; 

if strcmp(job.contrastType,'Go')
    zlim = 2;
else
    zlim = 1;
end

figure('Position',[100 100 1200 800]); hold on
contourf(data.mu.time,data.mu.freq,data.TF2plot,40,'linestyle','none'); hold on
if strcmp(job.lockSettings,'stimlocked')
    set(gca,'xlim',[job.TFtiming],'ylim',[1.2 15],'clim',[-1*zlim 1*zlim],'yscale','lin',... % yscale log
    'xtick',[0 0.5 0.815 1 1.3 1.5 2],'xtickLabel',{'Cue','0.5','AvgRT','1','+','1.5','2'},'ytick',[2 4 8 12 15],...
        'fontsize',32,'clim',[-1*zlim 1*zlim],'Linewidth',3) % -.25 1.3
    plot([0 0],get(gca,'ylim'),':k','LineWidth',3);
    plot([0.815 0.815],get(gca,'ylim'),':k','LineWidth',3);

elseif strcmp(job.lockSettings,'resplocked')
    set(gca,'xlim',[job.TFtiming],'ylim',[1.2 15],'clim',[-1*zlim 1*zlim],'yscale','lin',...
    'xtick',[-1:0.5:2],'xtickLabel',{'-1','-0.5','Resp','0.5','1','1.5','2'},'ytick',[2 4 8 12 15],... % yscale log
        'fontsize',32,'Linewidth',3) % -.25 1.3
    plot([0 0],get(gca,'ylim'),':k','LineWidth',3);

else
    error('Unknown lock setting')
end

colorbar('Ticks',[(-1*zlim):(zlim/2):zlim],'Fontsize',32)
ylabel('Frequency (Hz)','fontsize',32,'fontweight','bold');
xlabel('Time (s)','fontsize',32,'fontweight','bold');

% Save:
if job.stimERPcor
    saveas(gcf,sprintf('/project/3017042.02/Log/CueLockedPaperPlots/EEGTF_%s_%s_%s_zlim%d_induced.png',job.contrastType,job.lockSettings,cmap,zlim))
elseif ~isempty(job.invalidSubs)
    saveas(gcf,sprintf('/project/3017042.02/Log/CueLockedPaperPlots/EEGTF_%s_%s_%s_zlim%d_without6.png',job.contrastType,job.lockSettings,cmap,zlim))
else
    saveas(gcf,sprintf('/project/3017042.02/Log/CueLockedPaperPlots/EEGTF_%s_%s_%s_zlim%d.png',job.contrastType,job.lockSettings,cmap,zlim))
end
close gcf

%% Figure S01C: EEG Topoplots

% Settings:
addpath(fullfile(rootdir,'/Analyses/colorMaps'));
set(0, 'DefaultFigureColormap', jet(64)); cmap = 'jet';  
% set(0, 'DefaultFigureColormap', parula(64)); cmap = 'parula';  
% set(0, 'DefaultFigureColormap', redblue(64)); cmap = 'redblue'; 
% set(0, 'DefaultFigureColormap', turbo); cmap = 'turbo'; 

lineWidth = 5;

% Time:
sigTime = job.sigTime;

% z-axis:
if strcmp(job.contrastType,'Go')
    zlim = 2;
else
    zlim = 1;
end

% Start figure:
figure('Position',[100 100 1200 800]); hold on
cfg = []; cfg.figure = gcf; cfg.ylim = job.freq; cfg.zlim = [-1*zlim 1*zlim]; cfg.marker='on'; cfg.style='straight';
cfg.layout = 'easycapM11.mat'; cfg.comment = 'no'; cfg.xlim = sigTime;
cfg.highlight = 'on'; cfg.highlightchannel = job.channels; cfg.highlightsymbol = 'o'; cfg.highlightsize = 12; cfg.highlightfontsize = 12;
cfg.colorbar = 'yes'; % want 'no', i.e. do it yourself
ft_topoplotTFR(cfg,data.topo2plot);

% Linewidth:
% 1) ring, 2) left ear, 3) right ear, 4) nose, 5) line around, 6) small crosses, 
% 7) ring, 8) left ear, 9) right ear, 10) nose, 11) line around, 12) highlighted crosses,
% 13) left ear, 14) right ear, 15) nose, 16) around

% Increase channel highlights:
p = gca;
for i = 12:16
    t = p.Children(i); 
    t.LineWidth = lineWidth;
end

% Save:
if job.stimERPcor
    saveas(gcf,sprintf('/project/3017042.02/Log/CueLockedPaperPlots/EEGTopo_%s_%s_%d-%dms_%s_zlim%d_induced.png',job.contrastType,job.lockSettings,sigTime(1)*1000,sigTime(end)*1000,cmap,zlim))
elseif ~isempty(job.invalidSubs)
    saveas(gcf,sprintf('/project/3017042.02/Log/CueLockedPaperPlots/EEGTopo_%s_%s_%s_%d-%dms_%s_zlim%d_without6.png',job.contrastType,job.lockSettings,job.band,sigTime(1)*1000,sigTime(end)*1000,cmap,zlim))
else
    saveas(gcf,sprintf('/project/3017042.02/Log/CueLockedPaperPlots/EEGTopo_%s_%s_%s_%d-%dms_%s_zlim%d.png',job.contrastType,job.lockSettings,job.band,sigTime(1)*1000,sigTime(end)*1000,cmap,zlim))
end
close gcf

% END.
