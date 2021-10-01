function meanAbsValue = taft_postprocess_TF_allSubjects(job,dirs,sortBetas,iROI,selChans,zlim,nCol,isSave)

% meanAbsValue = taft_postprocess_TF_allSubjects(job,dirs,sortBetas,iROI,selChans,zlim,nCol,isSave)
%
% Plot 2D (frequency x time) t-values given selected channels for each
% subject in one grid of subplots (similar to taft_postprocess_TF_TFplot.m,
% but no test and no p-value highlighting). 
% Handy to detect subjects that are outliers.
% Works after loading a sortBetas object.
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
% frequency/ time bin (with channels sorted) for each subject, for selected
% subjects.
% iROI          = numeric scalar, index of selected ROI (to retrieve correct name).
% selChans      = cell with strings, selected channel labels.
% thresh        = numeric scalar, T-value threshold for cluster-based permutation
%               test, default is 2.
% zlim          = numeric scalar, limit for color axis (-zlim zlim) in units of T,
%               default is thresh.
% nCol          = integer scalar, number of columns in grid.
% isSave        = Boolean, save plot (true) or not (false).
% 
% OUTPUTS:
% No output, just plotting.
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
    selChans = {'FCz','Cz'}; % 
    fprintf('No channels selected -- use %s\n',strjoin(selChans,'/'));
end
if nargin < 6
    zlim = .0001;
    fprintf('No zlim specified -- use %.4f\n',zlim)
end
if nargin < 7
    nCol = 6;
    fprintf('No number of columns specified -- use %d\n',nCol)
end
if nargin < 8
    isSave = false;
    fprintf('Do not save plot by default\n')
end

%% Select channels, collapse over channels, reshape into 4-D:

chanBetas       = cell(length(sortBetas),1); % initialize objects for averaged channels (not: can have different time/ frequency axis!)
cfg             = [];

% a) Select channels:
fprintf('Select channels %s\n',strjoin(string(selChans),', '));
cfg.channel     = selChans;
cfg.avgoverchan = 'yes';

% b) Select latency:
if strcmp(job.lockSettings,'stimlocked')
    cfg.latency = [0 1.3];
elseif strcmp(job.lockSettings,'resplocked')
    cfg.latency = [-1.0 0.5]; % maximum
else
    error('Unknown lock settings')
end
cfg.avgovertime = 'no';
fprintf('Select timing %.3f - %.3f\n',cfg.latency(1), cfg.latency(end));

% c) Select frequency:
% cfg.frequency = [1 15];
% cfg.frequency = [4 8];
% cfg.avgoverfreq = 'no';
% fprintf('Select freqencies %d - %d\n',cfg.frequency(1), cfg.frequency(end));

fprintf('ROI %s: Average over channels %s, reshape\n',char(job.regNames(iROI)),strjoin(selChans,'/'));

%% Intialize objects:

nSub            = length(sortBetas); % number subjects found
nRows           = ceil(nSub/nCol); % number rows following from number subjects and number columns
meanAbsValue    = nan(nSub,1); % mean weight across entire map to detect outliers

%% Plot each subject in one subplot:

figure('units','normalized','outerposition',[0 0 1 1]); hold on;

% Loop over subjects:
for iSub = 1:nSub
    
    fprintf('Plot subject %03d\n',iSub)
    
    % Select subject data:
    chanBetas{iSub} = ft_selectdata(cfg,sortBetas{iSub}); % average over selected channels
    
    % Plot:
    subplot(nRows,nCol,iSub)
    contourf(chanBetas{iSub}.time,chanBetas{iSub}.freq,squeeze(chanBetas{iSub}.powspctrm),30,'linestyle','none');
    % colorbar
    set(gca,'clim',[-1*zlim 1*zlim]);
    xlabel('Time (in s)','fontweight','bold','fontsize',6); 
    ylabel('Frequency (in Hz)','fontweight','bold','fontsize',6)
    title(sprintf('Subject %03d',iSub))
    
    % Save mean absolute b-value:
    tmp = squeeze(chanBetas{iSub}.powspctrm);
    meanAbsValue(iSub) = mean(abs(tmp(:)))*1000;
    
end

% Overall title:
sgtitle(sprintf('TF-plot per subject for ROI %s, %s-locked, GLM-type %s, channels %s',char(job.regNames(iROI)),job.lockSettings,job.HRFtype,strjoin(selChans,'/')),'fontsize',12)

%% Save as:

if isSave 
    fileName = sprintf('TFplot_%s_%s_%s_%dsec_%s.png',...
        char(job.regNames(iROI)),job.lockSettings,job.HRFtype,job.trialdur,strjoin(selChans,''));
    saveas(gcf,fullfile(dirs.TFplot,fileName));
end

end % end of function.