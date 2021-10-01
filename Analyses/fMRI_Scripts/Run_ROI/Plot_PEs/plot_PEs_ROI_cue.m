function plot_PEs_ROI_cue(ROIname,ROItitle,isPoint,yLim)

% plot_PEs_ROI_cue(ROIname,ROItitle,isPoint,yLim)
%
% Creates and saves bar plot for 4 action x valence conditions for a given 
% ROI, with or without individual per-subject-per-condition data points.
% Expected order of conditions: Go2Win, Go2Avoid, NoGo2Win, NoGo2Avoid
%
% INPUTS:
% ROIname       = string, name of ROI (will search for file of type
% 'ROIname.txt')
% ROItitle      = string, title of plot (default: same as ROIname)
% isPoint       = Boolean, add individual subject data points or not
% (default: false)
% yLim          = vector of two floats, adjust y axis limits manually
% (optional)
%
% OUTPUTS:
% Saves plot under 'Log/fMRI/fMRI_ROIs/GLM%s_ROIPlots'
%
% Mind to change dirs.root to your own directory structure
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

if nargin < 2
    ROItitle = ROIname; % default: file name as title
end

if nargin < 3
    isPoint = false; % default: no individual data points
end

% yLim gets checked below

%% Set GLM:

GLMID = '1';
fprintf('GLM is GLM%s\n',GLMID);

%% Directories:

fprintf('Initialize directories\n');

dirs.root   = '/project/3017042.02'; 
dirs.ROI    = fullfile(dirs.root,'Log/fMRI/fMRI_ROIs');
dirs.data   = fullfile(dirs.ROI,sprintf('GLM%s_ROIData',GLMID));
dirs.plot   = fullfile(dirs.ROI,sprintf('GLM%s_ROIPlots',GLMID));
dirs.plot   = fullfile(dirs.root,'Log/CueLockedPaperPlots');
if ~exist(dirs.plot,'dir'); mkdir(dirs.plot); end

%% Load data;

dataFile    = sprintf('%s.txt',ROIname);

fprintf('Load file %s\n', dataFile);
data        = load(fullfile(dirs.data,dataFile));

%% Construct mean and SE per condition:

fprintf('Correct SE per condition for between-subjects variance\n');

nCond       = size(data,1);
subMean     = nanmean(data,1); % average over conditions --> one value per subject
grandMean   = nanmean(subMean,1); % average over subjects --> grand mean (1 value)

% Initialize:
condMean    = nan(nCond,1);
condSE      = nan(nCond,1);

% Correct SE by subtracting overall mean per subject and adding back grand average; correct for # conditions
for iCond = 1:nCond
    condMean(iCond) = nanmean(data(iCond,:));
    condSE(iCond) = nCond/(nCond-1)*nanstd(data(iCond,:)-subMean+repmat(grandMean,size(subMean,1)))./sqrt(length(subMean));
end
% condMean = [0.05 -0.025 0.025 -0.050]; ROItitle = 'Hypothesis'; % initial hypothesis

%% Fixed settings:

xLoc        = [0.5 1.5 3.0 4.0]; % locations of bars on x-axis
ColMat      = [0 .6 .2; .8 0 0;0 .6 .2; .8 0 0]; % colors
lwd         = 3; 
capSize     = 12; 
sz          = 24; 
fontSize    = 35; 

yTick       = -1:.05:1; 
yLimMargin  = 0.01; % if ylim has to be adjust automatically

%% Start plot:

figure('units','normalized','outerposition',[0 0 1 1]); hold on % fullscreen
% figure('Position',[0 0 1200 800]); hold on

% Bar plots:
for iCond = 1:nCond
    bar(xLoc(iCond),condMean(iCond),.75,'FaceColor',ColMat(iCond,:)) % bar plot
end

% Error bars:
for iCond = 1:nCond
    errorbar(xLoc(iCond),condMean(iCond),condSE(iCond),'k','linestyle','none','linewidth',lwd,'Capsize',capSize) % error bars
end

% Scatter plots:
if isPoint
    for iCond = 1:nCond
        scatter(repmat(xLoc(iCond),1,size(data,2)),data(iCond,:),sz,[0 0 0],'filled'); % scatter plot
    end
    yTick = [-1:.1:1];
end

% Settings:
set(gca, 'xlim',[0 max(xLoc)+0.6], 'xtick',xLoc,'xticklabel',{'Go2Win','Go2Avoid','NoGo2Win','NoGo2Avoid'},...
    'ytick',yTick,'XColor',[0 0 0],'YColor',[0 0 0],...
    'FontName','Arial','FontSize',fontSize-10,'Linewidth',lwd)
   
if exist('yLim','var') % adjust y axis if provided as input argument
    set(gca,'ylim',yLim)
else
    fprintf('Adjust ylim based on extreme values\n');
    if isPoint
        set(gca,'ylim',[min(data(:)) - yLimMargin, max(data(:)) + yLimMargin]) % min and max of individual data points
    else        
        set(gca,'ylim',[min((condMean(:)) - condSE(:)) - yLimMargin, max((condMean(:) + condSE(:))) + yLimMargin]) % means
    end
end

xlabel('Executed action x Valence','FontSize',fontSize,'FontName','Arial','Color',[0 0 0]);
ylabel('Parameter estimates (a.u.)','FontSize',fontSize,'FontName','Arial','Color',[0 0 0]);
title(sprintf('%s',ROItitle),'FontSize',fontSize,'FontName','Arial','Color',[0 0 0]);
legend({'Win','Avoid'},'fontsize',fontSize); legend boxoff;

%% Save:

figName = sprintf('Matlab_Barplot_GLM%s_%s',GLMID,ROIname);
if isPoint
    figName = [figName '_withPoints'];
end
% saveas(gcf,fullfile(dirs.plot,'Matlab_Barplot_Hypothesis.jpg'))

fprintf('Save figure %s\n',figName);
saveas(gcf,fullfile(dirs.plot,[figName '.png']));
pause(3)
close gcf

end % end of function.
