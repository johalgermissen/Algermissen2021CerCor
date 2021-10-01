function [tg,corrp] = taft_postprocess_time_ERplot(job,dirs,sortBetas,iROI,selChans,thresh,nP,isSave)

% [tg,corrp] = taft_postprocess_time_ERplot(job,dirs,sortBetas,iROI,selChans,thresh,nP,isSave)
%
% For fMRI-EEG regression weights for selected ROI, 
% select data for given channels/ time range, perform cluster-based 
% permutation test for mean signal across selected channels, output p-value,
% create lineplot of T-values with with "significant clusters" highlighted.
% 
% INPUTS:
% job           = cell with necessary settings for timing settings and name
% of image file to be saved:
% .HRFtype      = string, whether to perform estimation of HRF amplitude for each trial separately ('trial') or in a GLM featuring all trials of a given block ('block').
% .lockSettings	= string, type of event-locked data to include, either 'stimlocked' or 'resplocked'.
% .trialdur 	= numeric scalar, trial duration to use when epoching upsampled BOLD data in seconds (recommended by Hauser et al. 2015: 8 seconds).
% .regNames 	= cell, one or more string(s) of all regressors (fMRI and behavior) to include in design matrix plus selected trials.
% dirs          = cell, directories where to save image file:
% .TFplot       = string, directory where to save TF plots.
% sortBetas     = cell, regression weights for each ROI/ channel/ 
% time bin (with channels sorted) for each subject, for selected
% subjects.
% iROI          = numeric scalar, index of selected ROI (to retrieve correct name).
% selChans      = cell with strings, selected channel labels.
% thresh        = numeric scalar, T-value threshold for cluster-based permutation
%               test, default is 2.
% nP            = numeric scalar, number of permutations in cluster-based 
%               permutation test, default is 10000.
% zlim          = numeric scalar, limit for color axis (-zlim zlim) in units of T,
%               default is thresh.
% isSave        = Boolean, save plot (true) or not (false).
%
% OUTPUTS:
% tg            = vector with T-value for each time bin (averaged over
%               channels).
% corrp         = vector with p-value for cluster to which time bin belongs.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Adapted from Tobias Hauser (https://github.com/tuhauser/TAfT).
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/CueLockedAnalyses/TAfT/

%% Complete settings:

if nargin < 4
    iROI = 1;
    fprintf('No ROI specified -- use first ROI\n')    
end

if nargin < 5
    selChans = {'FCz','Cz'}; % Go/NoGo contrast
    fprintf('No channels selected -- use %s\n',strjoin(selChans,'/'));
end

if nargin < 6
    thresh = 2.0;
    fprintf('No cluster-forming threshold specified -- use 2.0\n')
end

if nargin < 7
    nP = 10000;
    fprintf('No number of permutations specified -- use 10,000\n')
end

if nargin < 8
    isSave = false;
    fprintf('Save by default\n')
end

%% Fixed settings:

lineWidth   = 3;
fontSize    = 24;
tCrit       = thresh;
pCrit       = 0.05;

%% Select channels, collapse over channels, reshape into 3-D:

chanBetas = cell(length(sortBetas),1); % Note: chanBetas might have different time axes

c = []; % initialize empty 4-D object

cfg = []; % empty cfg file for later data selection

% a) Select channels:
cfg.channel     = selChans;
cfg.avgoverchan = 'yes';

% b) Select latency:
if strcmp(job.lockSettings,'stimlocked')
    job.lockName    = 'stim';
   cfg.latency      = [0 0.7]; % time where ERPs occur
%    cfg.latency  = [0 1.3]; % whole trial
elseif strcmp(job.lockSettings,'resplocked')
    job.lockName    = 'resp';
    cfg.latency     = [-0.5 0.5];
%    cfg.latency  = [-1 0.5]; % maximum possible
else
    error('Unknown lock settings')
end
cfg.avgovertime = 'no';

fprintf('ROI %s: Average over channels %s, reshape\n',char(job.regNames(iROI)),strjoin(selChans,'/'));

%% Reshape to 4D object:

for iSub = 1:length(sortBetas) % iSub = 1;

    chanBetas{iSub} = ft_selectdata(cfg,sortBetas{iSub}); % average over selected channels
    c(iSub,:,:)     = chanBetas{iSub}.avg; % bring into 4-D, with subject as first dimension

end

%% Cluster-based permutation test:

[corrp,tg] = clustertf(c,thresh,nP); % default test statistic
% [corrp,tg] = clustertf(c,thresh,nP,'nExtent');
% [corrp,tg] = clustertf(c,thresh,nP,'maxT');
% [corrp,tg] = clustertf(c,thresh,nP,'sumT');

%% Joint plots:

figure('units','normalized','outerposition',[0 0 1 1]); hold on;

% First subplot: heat map of t-values:

timeVec = chanBetas{1}.time;
Tvec    = squeeze(tg);
pvec    = squeeze(corrp);
% ylim    = 4;
ylim    = max(4,ceil(max(abs(Tvec)))); % 4 or rounded-up maximal t-value

subplot(1,2,1);
plot(chanBetas{iSub}.time,Tvec,'k-','linewidth',lineWidth); hold on
set(gca,'xlim',cfg.latency,'ylim',[-ylim ylim],'fontsize',fontSize);
xlabel('Time (in s)','fontweight','bold','fontsize',fontSize); 
ylabel('T-values','fontweight','bold','fontsize',fontSize)
title('T-values')

% Highlight time ranges of "significant clusters" with patches:
if any(Tvec < -1*tCrit & pvec < pCrit)
    sigTimeIdx  = (Tvec < -1*tCrit & pvec < pCrit); % negative
    sigTime     = find(sigTimeIdx==1);
    fprintf('Significant negative cluster from %.03f - %.03f sec.\n',...
        timeVec(sigTime(1)),timeVec(sigTime(end)));
    patch([timeVec(sigTimeIdx) fliplr(timeVec(sigTimeIdx))],[repmat(-2,1,sum(sigTimeIdx)) fliplr(Tvec(sigTimeIdx)')],'b'); hold on
end
if any(Tvec > tCrit & pvec < pCrit)
    sigTimeIdx  = (Tvec > tCrit) & (pvec < pCrit); % positive
    sigTime     = find(sigTimeIdx==1);
    fprintf('Significant positive cluster from %.03f - %.03f sec.\n',timeVec(sigTime(1)),timeVec(sigTime(end)));
    patch([timeVec(sigTimeIdx) fliplr(timeVec(sigTimeIdx))],[repmat(2,1,sum(sigTimeIdx)) fliplr(Tvec(sigTimeIdx)')],'b'); hold on
end

% Horizontal lines:
nTime   = length(chanBetas{iSub}.time);
plot(chanBetas{iSub}.time,repelem(2,nTime),'color',[0.5 0.5 0.5],'linestyle','--','linewidth', lineWidth-1); hold on
plot(chanBetas{iSub}.time,repelem(-2,nTime),'color',[0.5 0.5 0.5],'linestyle','--','linewidth', lineWidth-1);

% P-values:
subplot(1,2,2);
plot(chanBetas{iSub}.time,pvec,'k-','linewidth',lineWidth)
set(gca,'ylim',[0 1],'fontsize',fontSize); % hard-coded so far
xlabel('Time (in s)','fontweight','bold','fontsize',fontSize); 
ylabel('p-values','fontweight','bold','fontsize',fontSize)
title('p-values')

% Overall title:
sgtitle(sprintf('ER-plot for regressor %s, %s-locked, GLM-type %s, channels %s',...
    char(job.regNames(iROI)),job.lockName,job.HRFtype,strjoin(selChans,'/')),'fontsize',fontSize)

fprintf('Finished ER plot\n');

%% Save figure:

if isSave
    saveas(gcf,fullfile(dirs.TFplot,sprintf('ERplot_%s_%s_%s_%dsec_%s.png',...
        char(job.regNames(iROI)),job.lockName,job.HRFtype,job.trialdur,strjoin(selChans,''))))
end

end % end of function.
