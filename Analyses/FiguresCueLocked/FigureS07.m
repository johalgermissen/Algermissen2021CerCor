% FigureS07.m

% Plots for Figure S07 in supplementary material.
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

% Figure 3A, D: contrastType = 'Congruency';
% Figure 3B, E: contrastType = 'Go';
% Figure 3C, F: contrastType = 'Valence';

% Initialize job:

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

% Contrast of interest:
job.contrastType        = 'Congruency'; % 'Congruency' or 'Go' or 'GoLeftRight' or 'GoValence' or 'GoAccuracy' or 'Accuracy' or 'IncongruentAccuracy' or 'CongruentAccuracy' or 'CorrectIncongruentIncorrectCongruent'

% ----------------------------------------------------------------------- %
% Load and prepare data:
job         = time_update_job(job); % initialize job

[job, data] = time_load_data(job); % load data
[job, data] = time_prepare_generic_data(job,data); % prepare generic objects

job         = time_update_job(job); % update job
[job, data] = time_prepare_contrast_data(job,data); % prepare objects specific to contrast

%% Figure S07A. EEG TWOLINEPLOT

% Settings:
fontSize        = 32; % 24; % 12 
lineWidth       = 5; % actual lines in plot

iBaseline   = find(round(data.mu.time,3) == 0); % find baseline via time in sec
baseline    = min(nanmean(data.mat1(:,iBaseline)),nanmean(data.mat2(:,iBaseline)));
until       = length(data.mu.time);

% Plot:
figure('Position',[100 100 1000 800]);  % short

p{1} = boundedline(data.mu.time(1:until),...
    nanmean(data.mat1(:,1:until))-baseline,...
    2/(2-1)*nanstd(data.mat1(:,1:until)-data.matMean(:,1:until)+repmat(data.matGrandMean(1:until),job.nValidSubs,1))./sqrt(job.nValidSubs),...
    'cmap',job.twoLineColor(1,:),'alpha','transparency',0.2); %
set(p{1},'linewidth',lineWidth) % adjust linewidth (not possible within boundedline)
p{2} = boundedline(data.mu.time(1:until),...
    nanmean(data.mat2(:,1:until))-baseline,...
    2/(2-1)*nanstd(data.mat2(:,1:until)-data.matMean(:,1:until)+repmat(data.matGrandMean(1:until),job.nValidSubs,1))./sqrt(job.nValidSubs),...
    'cmap',job.twoLineColor(2,:),'alpha','transparency',0.2); %
set(p{2},'linewidth',lineWidth) % adjust linewidth (not possible within boundedline)

% Settings:
set(gca,'xlim',[-0.25 1.3],'xtick',[0 0.5 0.815 1 1.3],'xtickLabel',{'Cue','0.5','AvgRT  ','1','+'},...
    'ylim', [-0.025 0.025], ...% [-0.017,0.02],
    'fontsize',28,'Linewidth',3) %,'ytick',-5:0.1:4)
set(gca,'ylim',[-0.02 0.032]); % good for ERPs
ylabel('Amplitude (A/cm²)','fontweight','bold','fontsize',32); 
xlabel('Time (s)','fontweight','bold','fontsize',32)
plot([0 0],get(gca,'ylim'),':k','LineWidth',3); % cue onset
plot([.815 .815],get(gca,'ylim'),':k','LineWidth',3); % average RT 

% Legend:
legend([p{:}],job.twoLineLabels,'fontsize',fontSize); legend boxoff

% Save:
saveas(gcf,sprintf('/project/3017042.02/Log/CueLockedPaperPlots/ERPs_%s_%s.png',job.contrastType,strjoin(sort(job.channels),'')))
pause(1)
close gcf

%% Figure S07A. EEG TOPOPLOT

startTimeVec = 0.0:0.10:0.50; % 12 panels
endTimeVec  = startTimeVec+(startTimeVec(2)-startTimeVec(1));
nCol        = 3;
nRows       = ceil(length(startTimeVec)/nCol);

% Create topoplot:
figure('Position',[100 100 1000 800]);  % short
for iPlot = 1:length(endTimeVec)
    subplot(nRows,nCol,iPlot)
    cfg = []; cfg.figure = gcf; cfg.marker='on'; cfg.style='straight'; % cfg.zlim = [-1*zlim 1*zlim]; 
    cfg.layout = 'easycapM11.mat'; cfg.comment = 'no'; cfg.xlim = [startTimeVec(iPlot) endTimeVec(iPlot)];
    if strcmp(job.contrastType,'BOLD')
        cfg.zlim = [-0.010 0.010]; % for BOLD effects
    else        
        cfg.zlim = [-0.015 0.015]; % for task effects
    end
    cfg.colorbar = 'no'; % want 'no', i.e. do it yourself % --> add externally
    ft_topoplotER(cfg,data.topo2plot);
%     title(sprintf('%.2f to %.2f sec',startTimeVec(iPlot),endTimeVec(iPlot)),'fontsize',32) % 18 for 2 digit
    title(sprintf('%.1f to %.1f sec', startTimeVec(iPlot),endTimeVec(iPlot)),'fontsize',24) % 28 for 1 digit
end
% Save:
saveas(gcf,sprintf('/project/3017042.02/Log/CueLockedPaperPlots/ERPTopoplot_%s.png',job.contrastType));
pause(1)
close gcf

%% Figure S07B. EEG MULTILINEPLOT

% Settings:
fontSize        = 32; % 24; % 12 
fontSizeTitle   = 15; % keep overall title at 15
lineWidth       = 3; % actual lines in plot

% Baseline:
iBaseline   = round(data.mu.time,3) == 0;
baseline    = min(squeeze(nanmean(data.SubCondTime(:,:,iBaseline))));

% Plot:
figure('Position',[100 100 1600 800]);  % long

% 4 lines:
for iCond = 1:job.nCond
        p{iCond} = boundedline(data.mu.time(1:until),...
            squeeze(nanmean(data.SubCondTime(job.validSubs,iCond,1:until)))-baseline,...
            job.nCond/(job.nCond-1)*squeeze(nanstd(squeeze(data.SubCondTime(job.validSubs,iCond,1:until))-data.SubTime(job.validSubs,1:until)+repmat(data.GrandTime(1:until),job.nValidSubs,1)))'./sqrt(job.nValidSubs),...
            'cmap',job.colMat(iCond,:),'alpha','transparency',0.2); %
        set(p{iCond},'linewidth',lineWidth)
end
set(p{3},'linestyle','--');set(p{4},'linestyle','--'); % NoGo conditions

% Settings:
set(gca,'xlim',[-0.25 1.3],'xtick',[0 0.5 0.815 1 1.3],'xtickLabel',{'Cue','0.5','AvgRT  ','1','+'},...
    'ylim', [-0.025 0.025], ...
    'fontsize',28,'Linewidth',3) 
set(gca,'ylim',[-0.02 0.035]); 
ylabel('Amplitude (A/cm²)','fontweight','bold','fontsize',32); 
xlabel('Time (s)','fontweight','bold','fontsize',32)
plot([0 0],get(gca,'ylim'),':k','LineWidth',3); % cue onset
plot([.815 .815],get(gca,'ylim'),':k','LineWidth',3); % average RT 
ylim = get(gca,'ylim');
% Title:
title('ERPs for the different action x valence conditions','fontsize',fontSize);
% Legend:
legend([p{:}],job.condNames,'fontsize',fontSize); legend boxoff

% Save:
saveas(gcf,sprintf('/project/3017042.02/Log/CueLockedPaperPlots/ERPs_Conditions_%s.png',strjoin(sort(job.channels),'')))
pause(1)
close gcf

% END.
