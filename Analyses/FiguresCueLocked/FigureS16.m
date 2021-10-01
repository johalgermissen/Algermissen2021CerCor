% FigureS16.m

% Plots for Figure S16 in supplementary material.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/FiguresCueLocked

% Set root directory:
rootdir = figures_set_rootdir(); % '/project/3017042.02';
addpath(fullfile(rootdir,'Analyses/EEG_Scripts/CueLockedAnalyses/TAfT'));

%% Load data:

EEGdomain   = 'time';

ROIs2use    = {'GLM1StriatumAction','GLM1CingulateAnteriorAction','GLM1LeftMotorHand','GLM1RightMotorHand','GLM1vmPFCValenceMan'}; % final paper

% 2) Behavioral regressors:
behav2use   = {'isgo','iswin','isW2G'}; % main effects and interaction --> final paper

% 3) Perform only on selected trials:
selTrials   = 'all';

[job, ~, betas] = taft_postprocess_load_job(EEGdomain,ROIs2use,behav2use,selTrials);

%% Figure S16 A; B: Line plot 

iROI        = 5; % vmPFC
selChans = {'FCz','Cz'}; % Figure S15 A
% selChans = {'F1','F3','FCz','FC1','FC3','Cz','C1','C3'}; % Figure S15 B: identified as P2 modulation

% Settings:
thresh      = 2.0;
nP          = 1000;
ylim        = 5;
lineWidth   = 4;
fontSize    = 32;
colMat      = [0 0 1];

% Compute:
[sortBetas,~]   = taft_postprocess_time_selectData(job,betas,iROI);
sortBetas       = taft_postprocess_FisherZ(sortBetas); % Perform Fisher z-transform:

c               = [];
cfg             = [];
cfg.channel     = selChans;
cfg.avgoverchan = 'yes';
cfg.latency  = [0 0.7]; % time where ERPs occur
cfg.avgovertime = 'no';
fprintf('Regressor %s: Average over channels %s, reshape into 3D\n',char(job.regNames(iROI)),strjoin(selChans,'/'));
for iSub = 1:length(sortBetas) % iSub = 1;
    chanBetas{iSub} = ft_selectdata(cfg,sortBetas{iSub}); % average over selected channels
    c(iSub,:,:)     = chanBetas{iSub}.avg; % bring into 4-D, with subject as first dimension
end
[corrp,tg]  = clustertf(c,thresh,nP);

% Plot:
figure('Position',[100 100 1200 800]); 
endTime     = cfg.latency(end);
timeVec     = sortBetas{1}.time(1:find(sortBetas{1}.time==endTime));
Tvec = squeeze(tg);
pvec        = squeeze(corrp);
plot(chanBetas{1}.time,Tvec,'color',colMat,'linewidth',lineWidth); hold on

if any(pvec < .05 & Tvec < -2)
    sigTimeIdx  = (Tvec < -2 & pvec < 0.05); % negative
    sigTime     = find(sigTimeIdx==1);
    fprintf('Significant negative cluster from %.03f - %.03f sec.\n',timeVec(sigTime(1)),timeVec(sigTime(end)));
    patch([timeVec(sigTimeIdx) fliplr(timeVec(sigTimeIdx))],[repmat(-2,1,sum(sigTimeIdx)) fliplr(Tvec(sigTimeIdx)')],colMat,'EdgeColor','none'); hold on
end

if any(pvec < .05 & Tvec > 2)
    sigTimeIdx  = (Tvec > 2) & (pvec < 0.05); % positive
    sigTime     = find(sigTimeIdx==1);
    fprintf('Significant positive cluster from %.03f - %.03f sec.\n',timeVec(sigTime(1)),timeVec(sigTime(end)));
    patch([timeVec(sigTimeIdx) fliplr(timeVec(sigTimeIdx))],[repmat(2,1,sum(sigTimeIdx)) fliplr(Tvec(sigTimeIdx)')],colMat,'EdgeColor','none'); hold on
end

% Horizontal lines:
nTime   = length(chanBetas{iSub}.time);
plot(chanBetas{iSub}.time,repelem(2,nTime),'color',[0.5 0.5 0.5],'linestyle','--','linewidth', 2);
plot(chanBetas{iSub}.time,repelem(-2,nTime),'color',[0.5 0.5 0.5],'linestyle','--','linewidth', 2);

% Final settings:
set(gca,'xlim',cfg.latency,'xtick',0:0.1:1,'fontsize',fontSize); % hard-coded so far
set(gca,'ylim',[-ylim ylim],'ytick',-10:2:10,'fontsize',fontSize); % hard-coded so far
xlabel('Time (in s)','fontweight','bold','fontsize',fontSize); 
ylabel('T-values','fontweight','bold','fontsize',fontSize); 
% title(sprintf('T-values over %s',strjoin(selChans,'/')));

% Save:
saveas(gcf,sprintf('/project/3017042.02/Log/CueLockedPaperPlots/TAfT_ERplot_%s_%s_%dsec_%s_ylim%d.png',...
    strjoin(job.regNames(iROI),''),job.HRFtype,job.trialdur,strjoin(selChans,''),ylim))
pause(2)
close gcf
fprintf('Done :-)\n')

%% Figure S16 C: TOPOPLOT BASED ON TAfT:

% Get t-values:
[~,Tvalues]     = taft_postprocess_time_selectData(job,betas,iROI);

% Settings:
zlim            = 3;
nCol            = 4; % number of columns in subplot

if strcmp(job.lockSettings,'stimlocked')
    job.lockName= 'stim';
elseif strcmp(job.lockSettings,'resplocked')
    job.lockName= 'resp';
else
    error('Unknown lock settings\n')
end    

% Timing settings:
startTime   = 0.0;
endTime     = 0.5;
steps       = 0.1;
startTimeVec    = startTime:steps:endTime; % vector of start time for each subplot
endTimeVec      = (startTimeVec(1)+steps):steps:(startTimeVec(end)+steps);  % vector of end time for each subplot
nRow            = ceil(length(startTimeVec)/nCol); % number of rows following from bins

% Create topoplot:
figure('Position',[100 100 1000 800]);  % short
for iPlot = 1:length(endTimeVec) % loop over given time windows

    subplot(nRow,ceil(length(endTimeVec)/nRow),iPlot)
    
    cfg = []; cfg.figure = gcf; cfg.zlim = [-1*zlim 1*zlim]; cfg.marker='on'; cfg.style='straight';
    cfg.layout = 'easycapM11.mat'; cfg.comment = 'no'; cfg.xlim = [startTimeVec(iPlot) endTimeVec(iPlot)];
    cfg.colorbar = 'no'; % want 'no', i.e. do it yourself % --> add externally
    
    ft_topoplotER(cfg,Tvalues);
    
    title(sprintf('%.1f to %.1f sec', startTimeVec(iPlot),endTimeVec(iPlot)),'fontsize',24) % 1 digit: 28
%     title(sprintf('%.2f to %.2f sec', startTimeVec(iPlot),endTimeVec(iPlot)),'fontsize',24) % 2 digits: 18

end % end iPlot

% Save:
saveas(gcf,fullfile(rootdir,sprintf('Log/CueLockedPaperPlots/Topoplot_time_%s_%s_%s_%dsec_%s-%sms.png',...
    char(job.regNames(iROI)),job.lockSettings,job.HRFtype,job.trialdur,...
    num2str(1000*startTimeVec(1)),num2str(1000*endTimeVec(end)))));
pause(1)
close gcf

%% Figure S16 D, E: LINEPLOT BASED ON EEG BINNED BY fMRI:

% Figure S16 D: job.chanArea = 'midfrontal'; 
% Figure S16 E: job.chanArea = 'leftfrontal'; 

rootdir = '/project/3017042.02';

addpath(fullfile(rootdir,'Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Helpers'));
addpath(fullfile(rootdir,'Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel'));

dirs        = set_dirs(rootdir);
par         = set_par();

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
% job.baselineTimings     = 0;

job.responseSettings    = 'Go'; % 'Go' or 'Hand'  
job.actionSettings      = 'exeAct'; % 'reqAct' or 'exeAct'
job.accSettings         = 'correct'; % 'bothAcc' or 'correct' or 'incorrect' or 'noAcc' or 'reqAct'
% for bothAcc, always use reqAct

% Channel settings:
job.chanArea            = 'midfrontal'; % Figure S16 D
% job.chanArea            = 'leftfrontal'; % Figure S16 E

% Contrast of interest:
job.contrastType        = 'Valence'; % 'Congruency' or 'Go' or 'GoLeftRight' or 'GoValence' or 'GoAccuracy' or 'Accuracy' or 'IncongruentAccuracy' or 'CongruentAccuracy' or 'CorrectIncongruentIncorrectCongruent'

job.bin.Type            = 'GLM1vmPFCValenceMan'; % RT GLM1vmPFCValenceMan
job.bin.Num             = 2; % always 2 for BOLD

% ----------------------------------------------------------------------- %
% Load and prepare data:
job         = time_update_job(job); % initialize job

[job, data] = time_load_data(job); % load data
[job, data] = time_prepare_generic_data(job,data); % prepare generic objects

job         = time_update_job(job); % update job
[job, data] = time_prepare_contrast_data(job,data); % prepare objects specific to contrast

%% Line plot:

% Settings:
fontSize        = 32; % 24; % 12 
lineWidth       = 5; % actual lines in plot

iBaseline       = find(round(data.mu.time,3) == 0); % find baseline via time in sec
baseline        = min(nanmean(data.mat1(:,iBaseline)),nanmean(data.mat2(:,iBaseline)));
until           = length(data.mu.time);

% Plot:
figure('Position',[100 100 1000 800]);  % short
p       = {};
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
saveas(gcf,sprintf('/project/3017042.02/Log/CueLockedPaperPlots/ERPs_%s_%s.png',job.bin.Type,strjoin(sort(job.channels),'')))
pause(1)
close gcf

%% Figure S16 F: TOPOPLOT BASED ON EEG BINNED BY fMRI:

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
saveas(gcf,sprintf('/project/3017042.02/Log/CueLockedPaperPlots/ERPTopoplot_%s.png',job.bin.Type));
pause(1)
close gcf

% END
