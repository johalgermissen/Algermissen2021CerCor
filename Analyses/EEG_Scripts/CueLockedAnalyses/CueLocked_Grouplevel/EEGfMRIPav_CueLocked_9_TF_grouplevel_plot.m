% EEGfMRIPav_CueLocked_9_TF_grouplevel_plot

% This is an interactive script---execute it step-by-step.
% Set the appropriate job settings for loading TF data aggregated across
% subjects, then create plots.
% - two-line plot contrasting two conditions of cotnrast
% - multi-line plot contrasting all conditions within data set
% - Topoplot contrasting two conditions of contrast
% - TF plot contrasting two conditions of contrast
% - these plots on the condition averaged signal
% - these plots for each subject separately
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
% job.bin.Type    = 'RT'; % RT GLM5FvmPFCValenceMan
% job.bin.Num     = 3; % only relevant for RT

% Load and prepare data:

job         = TF_update_job(job); % initialize job
[job, data] = TF_load_data(job); % load data
[job, data] = TF_prepare_generic_data(job,data); % prepare generic objects

job         = TF_update_job(job); % update job
[job, data] = TF_prepare_contrast_data(job,data); % prepare objects specific to contrast

%% 01A) LINEPLOT ON MAIN CONTRAST:

% Plot:
twoLinePlot(job,data)

% Save:
saveas(gcf,fullfile(dirs.TFplot,sprintf('twoLinePlot_%s.png',...
    job.plotName)));
pause(3)
close gcf

%% 01B) MULTILINEPLOT FOR ALL CUE CONDITIONS:

% Plot:
multiLinePlot(job,data)

% Manually adjust y-axis:
% set(gca,'ylim',[-0.2 1.6]); % stimulus-locked
% set(gca,'ylim',[-0.2 1.4]); % response-locked
% set(gca,'ylim',[-1 1]); % response-locked

saveas(gcf,fullfile(dirs.TFplot,sprintf('multiLinePlot_%s.png',...
    job.plotName)));
pause(3)
close gcf

%% 01C) SINGLE TOPO PLOTS:

% Settings:
zlim = 2;

% Plot:
singleTopoPlot(job,data,zlim)

% Save:
saveas(gcf,fullfile(dirs.TFplot,sprintf('Topoplot_%s_TOI_%03d_%03dms_zlim%d.png',...
    job.plotName, job.sigTime(1)*1000, job.sigTime(end)*1000, zlim)))
pause(3)
close gcf

%% 01D) MULTIPLE TOPOPLOTS OVER TIME:

% Settings:
zlim = 1; % (zlim=2 for Go)
% zlim = 2; % (zlim=2 for Go)

nRows = 2; startTime = 0; endTime = 0.9; steps = 0.1; % 2 rows stimlocked
% nRows = 3; startTime = 0; endTime = 1.4; steps = 0.1; % 3 rows stimlocked
% nRows = 2; startTime = -0.5; endTime = 0.4; steps = 0.1; % 2 rows resplocked
% nRows = 3; startTime = -0.9; endTime = 0.4; steps = 0.1; % 3 rows resplocked

% Plot:
multiTopoPlot(job,data,zlim,nRows,startTime,endTime,steps) % stim-locked 

% Save:
saveas(gcf,fullfile(dirs.TFplot,sprintf('multiTopoPlot_%s_TOI_%03d_%03dms_zlim%d.png',...
    job.plotName,1000*startTime,1000*(endTime+steps),zlim)));
close gcf

% Manual:
% multiTopoPlot(job,data,1,2) % stim-locked (zlim=2 for Go)
% multiTopoPlot(job,data,2,2) % stim-locked (zlim=2 for Go)
% multiTopoPlot(job,data,1,3,-1,0.4,0.1) % response-locked
% multiTopoPlot(job,data,2,3,-1,0.4,0.1) % response-locked (zlim=2 for Go)

%% 01E) TF CONTRAST PLOT:

% Settings:
zlim = 1;
% zlim = 2;

% Plot:
TFplot(job,data,zlim)

% Save:
saveas(gcf,fullfile(dirs.TFplot,sprintf('TFplot_%s_freq_1-15Hz_zlim%d.png', ...
    job.plotName,zlim)));
close gcf

%% 02) CONDITION-AVERAGED PLOTS.
%% 02A) CONDITION-AVERAGED TOPOPLOT for certain frequency band relative to baseline:
% use data.mu

fontSize    = 24;
lineWidth   = 5;

figure('name',sprintf('%s %s power %d-%d Hz',job.chanName,job.bandName,job.freq(1),job.freq(2)),'units','normalized','outerposition',[0 0 1 1]); hold on % fullscreen
% figure('name',sprintf('%s %s power %d-%d Hz',job.chanName,job.bandName,job.freq(1),job.freq(2)),'Position',[50 100 200 200])

cfg = []; 
% Settings:
cfg.figure = gcf; cfg.ylim = job.freq; cfg.zlim = [-5 -3];
cfg.layout = job.layout; cfg.comment = 'no';
cfg.xlim = job.sigTime; cfg.style = 'straight'; cfg.marker='on';
cfg.highlight = 'on'; cfg.highlightchannel = job.channels; 
cfg.highlightsymbol = 'o';
cfg.highlightsize = fontSize/2; cfg.highlightfontsize = fontSize/2;
% Plot:
ft_topoplotTFR(cfg,data.mu);
% Title:
title(sprintf('Grand average across conditions:\n %s (%d-%d Hz), %.3f to %.3f sec',...
    job.bandName, job.freq(1), job.freq(end), job.sigTime(1),job.sigTime(end)),'fontsize',fontSize)
% Color bar:
colorbar('fontsize',fontSize)

% Change marker line width:
p = gca;
for i = 12:16
    t = p.Children(i); 
    t.LineWidth = lineWidth;
end

% Save:
saveas(gcf,fullfile(dirs.TFplot,sprintf('singleTopoPlot_TOI_%s_%d-%dms_zlim_%d_%d.png',...
    job.plotName,1000*job.sigTime(1),1000*job.sigTime(end),cfg.zlim(1),cfg.zlim(end))));
pause(3)
close gcf

%% 02B) CONDITION-AVERAGED TOPOPLOT for certain frequency band relative to baseline LOOPING through time:

fontSize    = 24;

startTime = 0.0; endTime = 1.4; steps = 0.1; nRows = 3;
startTimeVec = startTime:steps:endTime;
endTimeVec = (startTimeVec(1)+steps):steps:(startTimeVec(end)+steps);

% Baseline:
baseIdx = dsearchn(data.mu.time',0);
data.mu.powspctrm = data.mu.powspctrm - repmat(data.mu.powspctrm(:,:,baseIdx),1,1,size(data.mu.powspctrm,3));

% Plot:
figure('units','normalized','outerposition',[0 0 1 1]); hold on % fullscreen
for iPlot = 1:length(endTimeVec)
    subplot(nRows,ceil(length(endTimeVec)/nRows),iPlot)
    cfg = []; 
    cfg.figure = gcf; cfg.ylim = job.freq; cfg.zlim = [-2 2]; cfg.marker='on'; cfg.style='straight';
    cfg.layout = job.layout; cfg.comment = 'no'; cfg.xlim = [startTimeVec(iPlot) endTimeVec(iPlot)];
    cfg.colorbar = 'no'; % want 'no', i.e. do it yourself % --> add externally
    % Plot:
    ft_topoplotTFR(cfg,data.mu);
    % Title:
    title(sprintf('%.2f to %.2f sec', startTimeVec(iPlot),endTimeVec(iPlot)),'fontsize',fontSize/2);
%     if iPlot == length(endTimeVec)
%         colorbar('fontsize',fontSize/2);
%     end
end
sgtitle(sprintf('Grand average across conditions: %s (%d-%d Hz)',...
    job.bandName, job.freq(1), job.freq(end)),'fontsize',fontSize)

saveas(gcf,fullfile(dirs.TFplot,sprintf('multiTopoPlot_%s_TOISequence_%d-%dms_zlim_%d_%d.png',...
    job.plotName,1000*startTime,1000*(endTime+steps),cfg.zlim(1),cfg.zlim(end))));
pause(3)
close gcf

%% 03) CONDITION-SPECIFIC PLOTS.
%% 03A) CONDITION-SPECIFIC TF-PLOT relative to pre-trial baseline.

fontSize    = 12;

if strcmp(job.lockSettings,'stimlocked')
    clim = [-1 2];
elseif strcmp(job.lockSettings,'resplocked')
    clim = [-2 1];
else
    error('Unknown lock setting')
end

% Start plot:
figure('name',sprintf('%s TF power all cue conditions',job.chanName),'units','normalized','outerposition',[0 0 1 1]); hold on % fullscreen
for iPlot = 1:job.nCond
    subplot(2,job.nCond/2,job.condOrder(iPlot)); hold on
    job.chanIdx = find(ismember(data.Cond{iPlot}.label, job.channels)); % determine indices of channels
    dat2plot = squeeze(mean(data.Cond{iPlot}.powspctrm(job.chanIdx,:,:))); % extract data
    % Baseline correction:
    baseTime = find(data.Cond{iPlot}.time==0);
    baseline = squeeze(mean(data.Cond{iPlot}.powspctrm(job.chanIdx,:,baseTime)));
    dat2plot = dat2plot - repmat(baseline',1,length(data.Cond{iPlot}.time)); % bring to 1
    % TF plot:
    contourf(data.mu.time,data.mu.freq,dat2plot,40,'linestyle','none')
    % Settings:
    set(gca,'xlim',[job.TFtiming],'xtick',[-1 -0.75 -0.5 -0.25 0 0.25 .5 0.75 1],...
        'ylim',[1 15],'ytick',2:2:14,'yscale','lin',... % yscale log
        'clim',[clim(1) clim(end)],'fontsize',fontSize,'clim',[clim(1) clim(end)])
    box off
    % Vertical lines:
    plot([0 0], get(gca,'ylim'),'k')
    plot(get(gca,'xlim'),[0 0],':k')
    plot([1.3 1.3], get(gca,'ylim'),'k')
    plot([2 2], get(gca,'ylim'),'k')
    plot([.815 .815], get(gca,'ylim'),':k')
    % Subplot title:
    title(job.condNames{iPlot},'fontsize',fontSize)
    % Axis labels:
    if iPlot == job.nCond/2+1 % bottom left
        ylabel('Frequency (Hz)','fontsize',fontSize);
        xlabel('Time','fontsize',fontSize);
%       Add legend:
%     elseif iPlot == job.nCond
%         colorbar
    end
end

% Save:
saveas(gcf,fullfile(dirs.TFplot,sprintf('TFplot_AllConditions_%s_freq_1-15Hz_zlim%d-%d.png', ...
    job.plotName,clim(1),clim(end))))
pause(3)
close gcf

%% 03B) CONDITION-SPECIFIC TOPO-PLOT relative to pre-trial baseline.

selTime = job.sigTime;
% selTime = [.5 .6]; % manually

% Start plot:
figure('Position',[600 100 600 500]); hold on
for iPlot = 1:job.nCond
    % Subplot:
    subplot(2,job.nCond/2,job.condOrder(iPlot))
    % Extract data;
    dat2plot = data.Cond{iPlot};
    % Prepare plot:
    cfg = []; cfg.figure = gcf; cfg.ylim = job.freq; cfg.zlim = [-5 -3]; % -2 2 ; -6 -2
    cfg.layout = job.layout; cfg.comment = 'no'; cfg.xlim = selTime;
    % Make plot:
    ft_topoplotTFR(cfg,dat2plot);
    % Title:
    title(job.condNames{iPlot})
    % Color bar:
%     if iPlot == job.nCond
%         colorbar
%     end
end

% Save:
saveas(gcf,fullfile(dirs.TFplot,sprintf('Topoplot_AllConditions_%s_TOI_%d-%dms_zlim%d-%d.png', ...
    job.plotName,cfg.xlim(1),cfg.xlim(end),clim(1),clim(end))));
pause(3)
close gcf

%% 04) SUBJECT-SPECIFIC PLOTS.
%% 04A) TWOLINEPLOT PER SUBJECT IN RASTER:

lineWidth   = 2;
fontSize    = 8;

% Start plot:
figure('name',sprintf('Line plot per subject: %s channels, %s band',job.chanName, job.band),'units','normalized','outerposition',[0 0 1 1]); hold on % fullscreen
for iSub = 1:job.nSub
    subplot(6,6,iSub)
    % Line plots:
    plot(data.mu.time,data.mat1(iSub,:),'color',job.twoLineColor(1,:),'linewidth',lineWidth); hold on
    plot(data.mu.time,data.mat2(iSub,:),'color',job.twoLineColor(2,:),'linewidth',lineWidth); 
    % Settings:
    set(gca,'xtick',[0 0.815 1.3 2 3],'xtickLabel',{'Cue','AvgRT','Fix','FB','ITI'},...
        'fontsize',fontSize,'xlim',[-.25 2])
%     set(gca,'ylim',[-4 4]) %  Consider adding y-lim to bring all subjects to same axes
    % Title:
    title(sprintf('Subject %03d',iSub),'fontsize',fontSize)
    box off; hold on
    % Vertical/horizontal lines:
    plot([0 0], get(gca,'ylim'),'k')
    plot(get(gca,'xlim'),[0 0],':k')
    plot([1.3 1.3], get(gca,'ylim'),'k')
    plot([.815 .815], get(gca,'ylim'),':k') % average response time
end
legend(job.twoLineLabels); legend boxoff

% Save:
saveas(gcf,fullfile(dirs.TFplot,sprintf('twoLinePlot_allSubjects_%s_%s.png',...
    job.plotName)));
pause(3)
close gcf

%% 04B) MULTILINEPLOT CUE CONDITIONS PER SUBJECT in subplot raster:

lineWidth   = 2;
fontSize    = 8;

% X-axis limits:
if strcmp(job.lockSettings,'stimlocked')
    xLim        = [-.25 1.3];
    xTick       = [-0.25 0 0.5 0.815 1.3];
    xTickLabel  = {'-0.25','0','0.50','AvgRT','1.3'};
elseif strcmp(job.lockSettings,'resplocked') 
    xLim        = [-1 1];
    xTick       = [-1 -0.5 0 .5 1];
    xTickLabel  = {'-1', '-0.5', 'Resp', '0.5', '1'};
else
    error('Unknown lock settings')
end

% Start plot:
figure('Position',[0 0 1500 1000]); hold on
for iSub = 1:job.nSub
    subplot(6,6,iSub); 
    p = cell(job.nCond,1);
    for iCond = 1:job.nCond
        p{iCond} = plot(data.mu.time,squeeze(data.SubCondTime(iSub,iCond,:))',...
            '-','Color',job.colMat(iCond,:),'linewidth',lineWidth); hold on
        set(p{iCond},'linestyle',job.lineStyle{iCond});
    end
    % Settings:
      set(gca,'xlim',xLim,'xtick',xTick,'xticklabel',xTickLabel,'ytick',-1*10:1:10,'fontsize',fontSize) % joint y-axis
%     set(gca,'ylim',[-4 4]) %  Consider adding y-lim to bring all subjects to same axes
    % Title:
    title(sprintf('Subject %03d',iSub),'fontsize',fontSize)
    % Vertical/horizontal lines:
    plot([0 0], get(gca,'ylim'),'k')
    plot([1.3 1.3], get(gca,'ylim'),'k')
    plot([.815 .815], get(gca,'ylim'),':k') % average response time
%     plot(get(gca,'xlim'),[0 0],'k-')
end
xlabel('Time (s)');
ylabel('Power (dB)'); 
legend([p{:}],job.condNames); legend boxoff

% Save:
saveas(gcf,fullfile(dirs.TFplot,sprintf('multiLinePlot_allSubjects_%s_%s.png',...
    job.plotName)));
pause(3)
close gcf

%% 04C) MULTILINEPLOT CUE CONDITIONS PER SUBJECT, loop with waitforbuttonpress:

lineWidth = 5;
fontSize = 12;

% X-axis limits:
if strcmp(job.lockSettings,'stimlocked')
    xLim        = [-.25 1.3];
    xTick       = [-0.25 0 0.5 0.815 1.3];
    xTickLabel  = {'-0.25','0','0.50','AvgRT','1.3'};
elseif strcmp(job.lockSettings,'resplocked') 
    xLim        = [-1 1];
    xTick       = [-1 -0.5 0 .5 1];
    xTickLabel  = {'-1', '-0.5', 'Resp', '0.5', '1'};
else
    error('Unknown lock settings')
end

baseIdx = dsearchn(data.mu.time',0);

for iSub = 1:job.nSub % iSub = 1
%     figure('units','normalized','outerposition',[0 0 1 1]); hold on % fullscreen
    figure('Position',[100 100 1200 600]); hold on
    p = cell(job.nCond,1);
    for iCond = 1:job.nCond
        p{iCond} = plot(data.mu.time,squeeze(data.SubCondTime(iSub,iCond,:))-data.SubCondTime(iSub,iCond,baseIdx),...
            '-','Color',job.colMat(iCond,:),'linewidth',lineWidth); hold on
        set(p{iCond},'linestyle',job.lineStyle{iCond});
    end
    % Settings:
    set(gca,'xlim',xLim,...
        'xtick',xTick,'xticklabel',xTickLabel,'ytick',-5:1:5,...
        'fontsize',fontSize)
%     set(gca,'ylim',[-3 7]) % consider keeping y-axis constant across people
    % Axis labels:
    xlabel('Time (s)','fontweight','bold');
    ylabel('Power (dB)','fontweight','bold');
    % Vertical lines:
    plot([.815 .815],get(gca,'ylim'),':k')
    % Title:
    title(sprintf('Subject %02d, %s %s (%d - %d Hz) power', ...
        iSub, job.chanName,job.band,job.freq(1),job.freq(2)), ...
        'fontsize',fontSize)
    legend([p{:}],job.condNames); legend boxoff
   w = waitforbuttonpress;
   % Save:
%    saveas(gcf,fullfile(dirs.TFplot,sprintf('multiLineplot_%s_sub%02d.png',job.plotName,iSub)))
   close gcf
end

%% 04D) MULTILINEPLOT CUE CONDITIONS without certain subject (LEAVE-ONE OUT):

lineWidth1  = 5;
lineWidth2  = 3;
fontSize    = 12;
transp      = 0.1;

% X-axis limits:
if strcmp(job.lockSettings,'stimlocked')
    xLim        = [-.25 1.3];
    xTick       = [-0.25 0 0.5 0.815 1.3];
    xTickLabel  = {'-0.25','0','0.50','AvgRT','1.3'};
elseif strcmp(job.lockSettings,'resplocked') 
    xLim        = [-1 1];
    xTick       = [-1 -0.5 0 .5 1];
    xTickLabel  = {'-1', '-0.5', 'Resp', '0.5', '1'};
else
    error('Unknown lock settings')
end

until       = find(round(data.mu.time,3) == 1.3); % end for cue-locked
iBaseline   = round(data.mu.time,3) == 0;
baseline    = min(squeeze(nanmean(data.SubCondTime(:,:,iBaseline))));

% Go effect robust
% Congruency effect in alpha: Influential cases: 6, 17, 25, 28, 34

for iSub = 1:job.nSub % iSub = 1
    
    % Option a: Drop only one subject:
%     dropSubs    = iSub;
%     selSubs     = setdiff(job.validSubs,iSub);
    % Option b: Drop k subjects    
    k           = 10; % how many to leave out
    dropSubs    = datasample(job.validSubs,k,'Replace',false); % Randomly select subjects to drop
    % Subject selection:
    selSubs     = setdiff(job.validSubs,dropSubs);
    fprintf('Drop subjects %s\n',strjoin(string(sort(dropSubs)),', '));

    % Plot:
    %     figure('units','normalized','outerposition',[0 0 1 1]); hold on % fullscreen
    figure('Position',[100 100 1200 600]); hold on
    p = cell(job.nCond,1);
    for iCond = 1:job.nCond
        p{iCond} = boundedline(data.mu.time(1:until),...
            squeeze(nanmean(data.SubCondTime(selSubs,iCond,1:until)))-baseline,...
            job.nCond/(job.nCond-1)*squeeze(nanstd(squeeze(data.SubCondTime(selSubs,iCond,1:until))-data.SubTime(selSubs,1:until)+repmat(data.GrandTime(1:until),length(selSubs),1)))'./sqrt(length(selSubs)),...
            'cmap',job.colMat(iCond,:),'alpha','transparency',transp); % hold on
        set(p{iCond},'linestyle',job.lineStyle{iCond}) % adjust style (solid, dashed)
        set(p{iCond},'Linewidth',lineWidth1) % adjust line width
    end
    % Settings:
    set(gca,'xlim',xLim,...
        'xtick',xTick,'xticklabel',xTickLabel,'ytick',-5:1:5,...
        'linewidth',lineWidth2,'fontsize',fontSize)
%     set(gca,'ylim',[-3 7]) % consider keeping y-axis constant across people
    % Axis labels:
    xlabel('Time (s)','fontweight','bold');
    ylabel('Power (dB)','fontweight','bold');
    % Vertical lines:
    plot([.815 .815],get(gca,'ylim'),':k')
    % Title:
    title(sprintf('Without subjects %s,\n %s %s (%d - %d Hz) power', ...
        strjoin(string(sort(dropSubs)),', '), job.chanName,job.band,job.freq(1),job.freq(2)), ...
        'fontsize',fontSize)
    legend([p{:}],job.condNames); legend boxoff
   w = waitforbuttonpress;
   % Save:
   saveas(gcf,fullfile(dirs.TFplot,sprintf('multiLineplot_%s_without_sub%s.png',...
       job.plotName,strjoin(string(sort(dropSubs)),'_'))))
   close gcf
end 

%% 04E) TOPOPLOT PER SUBJECT in subplot raster:
% Needs manual coding, thus change for respective contrast:

fontSize    = 10;

selTime     = [.5 .7]; % Go > NoGo; Go2Win > Go2Avoid
zLim        = 1;

% Start plot:
figure('units','normalized','outerposition',[0 0 1 1]); hold on % fullscreen
for iSub = 1:job.nSub % iSub = 1
    
    % Create topoplot data object per subject:
    % Needs manual coding, thus change for respective contrast:
    topo2plot = data.TFall{iSub}.ValxAct{1}; % initialize
    topo2plot.powspctrm = (data.TFall{iSub}.ValxAct{1}.powspctrm + data.TFall{iSub}.ValxAct{2}.powspctrm)./2 ...
        - (data.TFall{iSub}.ValxAct{3}.powspctrm + data.TFall{iSub}.ValxAct{4}.powspctrm)/2;

    % Plot:
    subplot(6,6,iSub)
    cfg = []; cfg.figure = gcf; cfg.ylim = job.freq; cfg.marker='on'; cfg.style='straight';
    cfg.layout = job.layout; cfg.comment = 'no'; cfg.xlim = selTime;
    cfg.colorbar = 'yes'; % want 'no', i.e. do it yourself
    cfg.zlim = [-1*zLim 1*zLim]; % set same z-axis for all subjects or not
    ft_topoplotTFR(cfg,topo2plot);
    title(sprintf('Subject %02d',iSub),'fontsize',fontSize)
end
sgtitle(sprintf('%s power %.1f to %.1f sec',...
    job.bandName, selTime(1), selTime(2)),'fontsize',fontSize*2)

saveas(gcf,fullfile(dirs.TFplot,sprintf('Topoplot_allSubjects_raster_%s_TOI_%03d_%03dms_zlim%d.png',...
    job.plotName,selTime(1),selTime(end),zLim)));
% pause(3)
% close gcf

%% 04F) TOPOPLOT PER SUBJECT, loop with waitforbuttonpress:
% Needs manual coding, thus change for respective contrast:

fontSize    = 12;

selTime     = [.5 .7]; % Go > NoGo; Go2Win > Go2Avoid

% Start plot:
for iSub = 1:job.nSub % iSub = 1
    
    figure('Position',[600 100 600 500]); hold on
    
    % Create topoplot data object per subject:
    % Needs manual coding, thus change for respective contrast:
    topo2plot = data.TFall{iSub}.ValxAct{1}; % initialize
    topo2plot.powspctrm = (data.TFall{iSub}.ValxAct{1}.powspctrm + data.TFall{iSub}.ValxAct{2}.powspctrm)./2 ...
        - (data.TFall{iSub}.ValxAct{3}.powspctrm + data.TFall{iSub}.ValxAct{4}.powspctrm)/2;

    % Plot:
    cfg = []; cfg.figure = gcf; cfg.ylim = job.freq; cfg.marker='on'; cfg.style='straight';
    cfg.layout = job.layout; cfg.comment = 'no'; cfg.xlim = selTime;
    cfg.highlight = 'on'; cfg.highlightjob.channel = job.channels; cfg.highlightsymbol = 'o';
    cfg.colorbar = 'yes'; % want 'no', i.e. do it yourself
    set(gca,'clim',[-1*zLim 1*zlim]);
    ft_topoplotTFR(cfg,topo2plot);
    title(sprintf('Subject %02d, %s power %.1f to %.1f sec',...
        iSub,job.band, selTime(1),selTime(2)),'fontsize',fontSize);
    w = waitforbuttonpress;
    % Save:
    saveas(gcf,fullfile(dirs.TFplot,sprintf('Topoplot_%s_TOI_%03d_%03dms_sub%02d.png',...
    job.plotName,selTime(1),selTime(end),iSub)));
    % pause(3)
    close gcf
end

%% 04G) TF PLOT PER SUBJECT in subplot raster:
% Needs manual coding, thus change for respective contrast:

lineWidth   = 2;
fontSize    = 8;

zlim        = 1;

% X-axis limits:
if strcmp(job.lockSettings,'stimlocked')
    xLim        = [-.25 1.3];
    xTick       = [-0.25 0 0.5 0.815 1.3];
    xTickLabel  = {'-0.25','0','0.50','AvgRT','1.3'};
elseif strcmp(job.lockSettings,'resplocked') 
    xLim        = [-1 1];
    xTick       = [-1 -0.5 0 .5 1];
    xTickLabel  = {'-1', '-0.5', 'Resp', '0.5', '1'};
else
    error('Unknown lock settings')
end

% Start plot:
figure('units','normalized','outerposition',[0 0 1 1]); hold on % fullscreen
for iSub = 1:job.nSub % iSub = 1
    
    % Create TF data object per subject:
    % Needs manual coding, thus change for respective contrast:
    TF2plot = squeeze(nanmean(data.TFall{iSub}.ValxAct{1}.powspctrm(job.chanIdx,:,:) + data.TFall{iSub}.ValxAct{2}.powspctrm(job.chanIdx,:,:) ...
            - data.TFall{iSub}.ValxAct{3}.powspctrm(job.chanIdx,:,:) - data.TFall{iSub}.ValxAct{4}.powspctrm(job.chanIdx,:,:)))/2;
            
    % Plot:
    subplot(6,6,iSub)
    contourf(data.mu.time,data.mu.freq,TF2plot,40,'linestyle','none'); hold on
    set(gca,'xlim',xLim,'xtick',xTick,'xticklabel',xTickLabel,'ylim',[1 15],'ytick',2:2:14,'yscale','lin',... % yscale log
        'fontsize',fontSize) % -.25 1.3
%     set(gca,'clim',[-1*zlim 1*zlim]); % same color axis for each subject
    colormap('jet'); colorbar
    % Axis labels:
    xlabel('Time (s)','fontsize',fontSize);
    ylabel('Frequency (Hz)','fontsize',fontSize);
    % Title:
    title(sprintf('Subject %02d',iSub),'fontsize',fontSize);
    
end
saveas(gcf,fullfile(dirs.TFplot,sprintf('TFplot_allSubjects_raster_%s_zlim%d.png',...
    job.plotName,zlim)));
close gcf

%% 04H) TF PLOT PER SUBJECT, loop with waitforbuttonpress:

fontSize    = 12;

zlim        = 1;

% X-axis limits:
if strcmp(job.lockSettings,'stimlocked')
    xLim        = [-.25 1.3];
    xTick       = [-0.25 0 0.5 0.815 1.3];
    xTickLabel  = {'-0.25','0','0.50','AvgRT','1.3'};
elseif strcmp(job.lockSettings,'resplocked') 
    xLim        = [-1 1];
    xTick       = [-1 -0.5 0 .5 1];
    xTickLabel  = {'-1', '-0.5', 'Resp', '0.5', '1'};
else
    error('Unknown lock settings')
end

for iSub = 1:job.nSub % iSub = 1
    
    % Create TF data object per subject:
    % Needs manual coding, thus change for respective contrast:
    TF2plot = squeeze(nanmean(data.TFall{iSub}.ValxAct{1}.powspctrm(job.chanIdx,:,:) + data.TFall{iSub}.ValxAct{2}.powspctrm(job.chanIdx,:,:) ...
            - data.TFall{iSub}.ValxAct{3}.powspctrm(job.chanIdx,:,:) - data.TFall{iSub}.ValxAct{4}.powspctrm(job.chanIdx,:,:)))/2;

    % Plot:
    figure('units','normalized','outerposition',[0 0 1 1]); hold on % fullscreen
    contourf(data.mu.time,data.mu.freq,TF2plot,40,'linestyle','none'); hold on
    set(gca,'xlim',xLim,'xtick',xTick,'xticklabel',xTickLabel,'ylim',[1 15],'ytick',2:2:14,'yscale','lin',... % yscale log
        'fontsize',fontSize) % -.25 1.3
%     set(gca,'clim',[-1*zlim 1*zlim]); % same color axis for each subject
    colormap('jet'); colorbar
    % Axis labels:
    xlabel('Time (s)','fontsize',fontSize);
    ylabel('Frequency (Hz)','fontsize',fontSize);
    % Vertical lines:
    plot([0 0],get(gca,'ylim'),'k-')
    plot([.815 .815],get(gca,'ylim'),':k')
    plot([1.3 1.3],get(gca,'ylim'),'k-')
    % Title:
    title(sprintf('Subject %02d, %s %s (%d - %d Hz) power',...
        iSub,job.chanName,job.band,job.freq(1),job.freq(2)),'fontsize',fontSize*2)
    w = waitforbuttonpress;
    % Save:
    saveas(gcf,fullfile(dirs.TFplot,sprintf('TFplot_%s_sub%03d.png',...
    job.plotName,iSub)));
    close gcf
end

%% 04I) TF PLOT LEAVE K SUBJECTS OUT:

for iSub = 1:job.nSub

    % Option a: Drop only one subject:
%     dropSubs    = iSub;
%     selSubs     = setdiff(job.validSubs,iSub);
    % Option b: Drop k subjects    
    k           = 10; % how many to leave out
    dropSubs    = datasample(job.validSubs,k,'Replace',false); % Randomly select subjects to drop
    % Subject selection:
    selSubs     = setdiff(job.validSubs,dropSubs);
    fprintf('Drop subjects %s\n',strjoin(string(sort(dropSubs)),', '));
    
    % Compute TF again by averaging over selected subjects:
    TFtmp   = cell(job.nCond,1);
    for iCondi = 1:job.nCond
        cfg              = [];
        cfg.parameter    = 'powspctrm';
        TFtmp{iCondi} = ft_freqgrandaverage(cfg,data.tmpSubCond{selSubs,iCondi}); % average over all conditions
    end    
    % Go contrast:   
    TF2plot = squeeze(nanmean(TFtmp{1}.powspctrm(job.chanIdx,:,:) + TFtmp{2}.powspctrm(job.chanIdx,:,:) ...
            - TFtmp{3}.powspctrm(job.chanIdx,:,:) - TFtmp{4}.powspctrm(job.chanIdx,:,:)))/2;    
    
    % Plot:
    figure('units','normalized','outerposition',[0 0 1 1]); hold on % fullscreen
    contourf(data.mu.time,data.mu.freq,TF2plot,40,'linestyle','none'); hold on
    set(gca,'xlim',xLim,'xtick',xTick,'xticklabel',xTickLabel,...
        'ylim',[1 15],'ytick',2:2:14,'yscale','lin',... % yscale log
        'fontsize',fontSize) % -.25 1.3
%     set(gca,'clim',[-1*zlim 1*zlim]); % same color axis for each subject
    colormap('jet'); colorbar
    % Axis labels:
    xlabel('Time (s)','fontsize',fontSize);
    ylabel('Frequency (Hz)','fontsize',fontSize);
    % Vertical lines:
    plot([0 0],get(gca,'ylim'),'k-')
    plot([.815 .815],get(gca,'ylim'),':k')
    plot([1.3 1.3],get(gca,'ylim'),'k-')
    % Title:
    title(sprintf('Without subjects %s,\n %s %s (%d - %d Hz) power',...
        strjoin(string(sort(dropSubs)),', '),job.chanName,job.band,job.freq(1),job.freq(2)),...
        'fontsize',fontSize)
    w = waitforbuttonpress;
    % Save:
    saveas(gcf,fullfile(dirs.TFplot,sprintf('TFplot_%s_without_sub%s.png',...
    job.plotName,strjoin(string(sort(dropSubs)),'_'))));
    close gcf
end

% END
