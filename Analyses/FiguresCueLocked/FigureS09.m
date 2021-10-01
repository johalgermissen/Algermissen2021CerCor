% FigureS09.m

% Plots for Figure S09 in manuscript.
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
job.accSettings         = 'bothAcc'; % 'bothAcc' or 'correct' or 'incorrect' or 'noAcc' or 'reqAct'
% for bothAcc, always use reqAct
job.TFtype              = 'hanning'; % hanning morlet
job.nFreq               = 15; % 15
% Channel settings:
job.chanArea            = 'midfrontal'; % 'midfrontal' or 'frontal' or 'central' or 'leftmotor' or 'rightmotor'
% Band settings:
job.band                = 'alpha'; % 'theta' or 'alpha' or 'beta' or 'broad' 
% Contrast of interest:
job.contrastType        = 'Congruency'; % 'Congruency' or 'Go' or 'GoLeftRight' or 'GoValence' or 'GoAccuracy' or 'Accuracy' or 'CongruentAccuracy' or 'IncongruentAccuracy'
% ERP-corrected:
job.stimERPcor = false;
job.respERPcor = false;

% Load and prepare data:
job         = TF_update_job(job); % initialize job
[job, data] = TF_load_data(job); % load data
[job, data] = TF_prepare_generic_data(job,data); % prepare generic objects
job         = TF_update_job(job); % update job
[job, data] = TF_prepare_contrast_data(job,data); % prepare objects specific to contrast

%% Bar plot:

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
for iCond = 1:nCond
    condMean(iCond) = nanmean(subCondMean(:,iCond));
    condSE(iCond) = nCond/(nCond-1)*nanstd(subCondMean(:,iCond)-subMean+repmat(grandMean,job.nValidSubs,1))./sqrt(job.nValidSubs);
end

% Settings:
colMat      = repmat([0 .6 .2; .8 0 0],job.nCond/2,1);
xLoc        = [1 2 3 4 5.5 6.5 7.5 8.5];
lineWidth   = 3;
yLim        = [-0.6 0.6]; yTick = [-1:.2:1];

% Figure:
close all
p   = {};
figure('Position',[100 100 1200 800]); hold on
for iCond = 1:job.nCond
    p{iCond} = bar(xLoc(iCond),condMean(iCond),.75,'FaceColor',colMat(iCond,:));
    errorbar(xLoc(iCond),condMean(iCond),condSE(iCond),'k','linestyle','none','Linewidth',lineWidth);
    if points
        s = scatter(repmat(xLoc(iCond),1,job.nValidSubs),subCondMean(:,iCond)',[],'k', 'jitter','on', 'jitterAmount',0.15); hold on % was 0.05
        set(s,'MarkerEdgeColor',[0.4 0.4 0.4],'linewidth',3); % was 1 
        yLim = [round(min(subCondMean(:)),2)-.01 round(max(subCondMean(:)),2)+.01];
        yTick = [-10:1:10];
    end
end
ylabel('Alpha Power','FontSize',32,'FontName','Arial','Color',[0 0 0]);
xlabel('Required action x Accuracy','FontSize',32,'FontName','Arial','Color',[0 0 0]);
set(gca,'FontName','Arial','FontSize',32,'xlim',[0 10],'TickLength',[0 0],...
    'xtick',[1.5 2.5 3.5 6 7 8],'xticklabel',{'Go','\newline Correct','NoGo','Go','\newline Incorrect','NoGo'})
set(gca,'ylim',yLim,'ytick',yTick,'YColor',[0 0 0]); %
legend([p{1},p{2}],{'Win','Avoid'}); legend boxoff
% title(sprintf('Valence x Executed Action x Accuracy'),'FontName','Arial','FontSize',15,'FontWeight','normal','Color',[0 0 0])

% Save:
if points
    saveas(gcf,'/project/3017042.02/Log/CueLockedPaperPlots/Barplot_Alpha_ValenceActionAccuracy_points.png')
else
    saveas(gcf,'/project/3017042.02/Log/CueLockedPaperPlots/Barplot_Alpha_ValenceActionAccuracy.png')
%     saveas(gcf,fullfile(dirs.plot,sprintf('Barplot_%s_baseCor_%s_baseTime_%d-%d_selTime_%d-%dms_selFreq_%d-%dHz_ValenceExecutedActionAccuracy.jpg',job.lockSettings,job.baselineSettings,job.baselineTimings(1)*1000,job.baselineTimings(end)*1000,selTime(1)*1000,selTime(end)*1000,job.freq(1),job.freq(end))))
end
close gcf

% END.
