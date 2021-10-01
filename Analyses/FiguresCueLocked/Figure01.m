% Figure01.m

% Plots for Figure 1 in manuscript.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/FiguresCueLocked

% Set root directory:
rootdir = figures_set_rootdir(); % '/project/3017042.02';

%% Figures 1A-C. NO DATA, JUST TASK DEPICTION. 

%% Figure 1. PREPROCESS BEHAVIOR

% Add additional paths:
addpath(fullfile(rootdir,'Analyses/Behavior_Scripts/Matlab_Plots'))

% Retrieve recoded and pre-processed behaviora data:
out = EEGfMRIPav_aggr_data(); % execute to get out

% Extract relevant data objects:
pGoCond         = out.pGoCond;
pCorrectCond    = out.pCorrectCond;
RTCond          = out.RTCond;
RTAcc           = out.RTAcc;
nSub            = size(pGoCond,1);
nCond           = size(pGoCond,2);
nCondRT         = size(RTAcc,2);
nRep            = size(pGoCond,3);
pSimp           = 0.05;
pBonf           = pSimp/nRep;

% P-values:
sigValence = nan(nSub,1); sigAction = nan(nSub,1); sigAccuracy = nan(nSub,1);
for iRep = 1:nRep
    difVec = pGoCond(:,1,iRep) - pGoCond(:,2,iRep) + pGoCond(:,3,iRep) - pGoCond(:,4,iRep);
    [~,P,~,~]  = ttest(difVec);
    sigValence(iRep) = P;
    difVec = pGoCond(:,1,iRep) + pGoCond(:,2,iRep) - pGoCond(:,3,iRep) - pGoCond(:,4,iRep);
    [~,P,~,~]  = ttest(difVec);
    sigAction(iRep) = P;
    difVec = mean(pCorrectCond(:,:,iRep),2)-0.33;
    [~,P,~,~]  = ttest(difVec);
    sigAccuracy(iRep) = P/2; % one-sided test
end

%% Subject selection:

invalidSubs     = []; % include all subjects
% invalidSubs     = [1 11 15 19 21 25]; % outliers in TAfT and fMRI
fprintf('Exclude subjects %s\n',strjoin(string(invalidSubs),', '));
validSubs       = setdiff(1:nSub,invalidSubs);
nSubValid       = length(validSubs);

%% Figure 1D. LINE PLOT OVER TRIALS.

% Select subjects, grand means:
subMean     = squeeze(nanmean(pGoCond(validSubs,:,:),2)); % average across conditions
grandMean   = squeeze(nanmean(subMean,1)); % average across subjects

% Average per condition:
condMean    = nan(nCond,nRep);
condSE      = nan(nCond,nRep);
for iCond = 1:nCond
    condMean(iCond,:) = squeeze(nanmean(pGoCond(validSubs,iCond,:),1));
    condSE(iCond,:)   = nCond/(nCond-1)*nanstd(squeeze(pGoCond(validSubs,iCond,:))-subMean+repmat(grandMean,nSubValid,1))./sqrt(nSubValid);
end

% Plot settings:
colMat      = [0 .6 .2; .8 0 0;0 .6 .2; .8 0 0];
fontSize    = 38;
xWidth      = 1200; % 1200

% Start plot:
close all
figure('Position',[0 0 xWidth 800]); 
for iCond = 1:nCond
    p{iCond} = boundedline(1:nRep,condMean(iCond,:),condSE(iCond,:),'cmap',colMat(iCond,:),'alpha');
    set(p{iCond},'Linewidth',4)
    if iCond > nCond/2
        set(p{iCond},'Linestyle','--');
    end
end

lineWidth = 2.5; markerSize = 15; yValence = 0.08; yAction = 0.08; yAccuracy = 0.05;
xLoc = find(sigValence < pSimp); nLoc = length(xLoc);
plot(xLoc,yValence*ones(nLoc,1),'*','Color',[1 0.8 0.6],'Linewidth',lineWidth,'MarkerSize',markerSize);
xLoc = find(sigValence < pBonf); nLoc = length(xLoc);
plot(xLoc,yValence*ones(nLoc,1),'*','Color',[1 0.5 0],'Linewidth',lineWidth,'MarkerSize',markerSize);
xLoc = find(sigAccuracy< pSimp); nLoc = length(xLoc);
plot(xLoc,yAccuracy*ones(nLoc,1),'*','Color',[0.6 0.8 1],'Linewidth',lineWidth,'MarkerSize',markerSize);
xLoc = find(sigAccuracy < pBonf); nLoc = length(xLoc);
plot(xLoc,yAccuracy*ones(nLoc,1),'*','Color',[0 0 1],'Linewidth',lineWidth,'MarkerSize',markerSize);

% Add plot features:
set(gca,'xlim',[0 nRep],'ylim',[0 1],'Linewidth',3)
set(gca,'xtick',0:10:nRep,'ytick',0:.2:1,'FontSize',fontSize,'Linewidth',4)
ylabel('p(Go)','FontSize',fontSize,'FontName','Arial');
xlabel('Trial','FontSize',fontSize,'FontName','Arial');
box off

% Save:
if isempty(invalidSubs)
    saveas(gcf,fullfile(rootdir,sprintf('Log/CueLockedPaperPlots/Behavior_CueConditions_%d.png',xWidth)))
else
    saveas(gcf,fullfile(rootdir,sprintf('Log/CueLockedPaperPlots/Behavior_CueConditions_withoutInvalid_%d.png',xWidth)))
end
% pause(2)
close gcf

%% Figure 1E. P(Go) BAR PLOTS.

% Select subjects, average across trials, grand means:
pGoCondMean         = squeeze(nanmean(pGoCond(validSubs,:,:),3)); % Go responses
pGoCorrectCondMean  = squeeze(nanmean(pCorrectCond(validSubs,:,:),3)); % correct responses
subMean             = squeeze(nanmean(pGoCondMean,2)); % average across conditions
grandMean           = squeeze(nanmean(subMean,1)); % average across subjects

% Mean per condition:
condMean            = nan(nCond,1);
condMeanCorrect     = nan(nCond,1);
condSE              = nan(nCond,1);
for iCond = 1:nCond
    condMean(iCond)         = squeeze(nanmean(pGoCondMean(:,iCond),1));
    condSE(iCond)           = nCond/(nCond-1)*nanstd(squeeze(pGoCondMean(:,iCond))-subMean+repmat(grandMean,nSubValid,1))./sqrt(nSubValid);
    condMeanCorrect(iCond)  = squeeze(nanmean(pGoCorrectCondMean(:,iCond),1));
end

% General plot settings:
colMat    = [0 .6 .2; .8 0 0; 0 .6 .2; .8 0 0];
colMatCor = [0.43 .75 .39; .95 .44 .17]; % 110 192 101; 243 112 43
xLoc      = [1 2 3.5 4.5];
lineWidth = 4; capSize = 12; fontSize = 36;

% Start plot:
close all
figure('Position',[100 100 800 800]); hold on

% a) Plot bars with errorbars:
for iCond = 1:nCond
    bar(xLoc(iCond),condMean(iCond),.75,'FaceColor',colMat(iCond,:)); % bar plot
    errorbar(xLoc(iCond),condMean(iCond),condSE(iCond),'k','linestyle','none','linewidth',lineWidth,'Capsize',capSize); % error bars
end

% b) Plot bars with correct Gos:
for iCond = 1:2
    bar(xLoc(iCond),condMeanCorrect(iCond),.75,'FaceColor',colMatCor(iCond,:)); % set(p1,'FaceAlpha',0.5);
end

% c) Points:
for iCond = 1:nCond
    s = scatter(repmat(xLoc(iCond),1,nSubValid),pGoCondMean(:,iCond)',[],'k', 'jitter','on', 'jitterAmount',0.15); hold on % was 0.05
    set(s,'MarkerEdgeColor',[0.4 0.4 0.4],'linewidth',3); % was 1 
end

% Add plot features:
set(gca,'xlim',[.5 5.5],'ylim',[0 1],'xtick',0:10:nRep,...
    'xtick',[1.5 4],'xticklabel',{'Go','NoGo'},'ytick',0:.2:1,...
    'FontSize',fontSize,'FontName','Arial','FontWeight','normal','Linewidth',4)
ylabel('p(Go)','FontSize',fontSize,'FontName','Arial');
xlabel('Required Action','FontSize',fontSize,'FontName','Arial');
box off

% Save:
if isempty(invalidSubs)
    saveas(gcf,fullfile(rootdir,'Log/CueLockedPaperPlots/Behavior_pGoBars.png'))
else
    saveas(gcf,fullfile(rootdir,'Log/CueLockedPaperPlots/Behavior_pGoBars_withoutInvalid.png'))
end
close gcf

%% Figure 1F. BAR PLOT REACTION TIMES.

% Select subjects, grand means:
RTCondMean          = RTAcc(validSubs,:); % select subjects
subMean             = squeeze(nanmean(RTCondMean,2)); % average across conditions
grandMean           = squeeze(nanmean(subMean,1)); % average across subjects

% Mean per condition per subject:
condMean            = nan(nCondRT,1);
condMeanCorrect     = nan(nCondRT,1);
condSE              = nan(nCondRT,1);
for iCond = 1:nCondRT
    condMean(iCond)         = squeeze(nanmean(RTCondMean(:,iCond),1));
    condSE(iCond)           = nCond/(nCond-1)*nanstd(squeeze(RTCondMean(:,iCond))-subMean+repmat(grandMean,nSubValid,1))./sqrt(nSubValid);
end

% General plot settings:
yLimMax         = round(max(RTCondMean(:)),2)+.01;
colMat          = [0.43 .75 .39; .95 .44 .17; 0 .6 .2; .8 0 0; 0 .6 .2; .8 0 0];
xLoc            = [1 2 3.5 4.5 6 7];
fontSize        = 38;

% Start plot:
close all
figure('Position',[100 100 1000 800]); hold on
for iCond = 1:nCondRT % loop over conditions to create bar and scatter plots
    bar(xLoc(iCond),condMean(iCond),.75,'FaceColor',colMat(iCond,:)); % bar plot
    errorbar(xLoc(iCond),condMean(iCond),condSE(iCond),'k','linestyle','none','linewidth',lineWidth,'Capsize',capSize); % error bars
    s = scatter(repmat(xLoc(iCond),1,nSubValid),RTCondMean(:,iCond)',[],'k', 'jitter','on', 'jitterAmount',0.15); hold on % was 0.05
    set(s,'MarkerEdgeColor',[0.4 0.4 0.4],'linewidth',3); % was 1 
end

% Add plot features:
set(gca,'xlim',[xLoc(1)-0.5 xLoc(end)+0.5],'ylim',[0 yLimMax],'xtick',0:10:nRep,...
    'xtick',[1.5 4 5.25 6.5],'xticklabel',{'correct Go','\newline(other Go)','incorrect Go','\newline(NoGo)'},'ytick',0:.2:2,...
    'FontSize',fontSize,'FontName','Arial','FontWeight','normal','Linewidth',4)
ylabel('Reaction time','FontSize',fontSize,'FontName','Arial','Color',[0 0 0]);
% xlabel('Performed (required) action','FontSize',fontSize,'FontName','Arial','Color',[0 0 0]);
box off

% Save:
if isempty(invalidSubs)
    saveas(gcf,fullfile(rootdir,'Log/CueLockedPaperPlots/Behavior_RTBars.png'))
else
    saveas(gcf,fullfile(rootdir,'Log/CueLockedPaperPlots/Behavior_RTBars_withoutInvalid.png'))
end
close gcf

% END.
