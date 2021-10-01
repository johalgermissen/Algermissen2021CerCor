function taft_postprocess_time_topoplot(job,dirs,Tvalues,iROI,startTime,endTime,steps,zlim,isSave) 

% taft_postprocess_TF_topoplot(job,dirs,Tvalues,iROI,startFreq,endFreq,startTime,endTime,steps,zlim,isSave)
% 
% For fMRI-EEG regression weights for selected ROI, 
% plot a topoplot for given time window ranges (in subplots).
%
% INPUTS:
% job           = cell with necessary settings for timing settings and name
% of image file to be saved:
% .HRFtype      = string, whether to perform estimation of HRF amplitude for each trial separately ('trial') or in a GLM featuring all trials of a given block ('block').
% .layout       = string, cap to use for topoplots.
% .lockSettings	= string, type of event-locked data to include, either 'stimlocked' or 'resplocked'.
% .trialdur 	= numeric scalar, trial duration to use when epoching upsampled BOLD data in seconds (recommended by Hauser et al. 2015: 8 seconds).
% .regNames 	= cell, one or more string(s) of all regressors (fMRI and behavior) to include in design matrix plus selected trials.
% dirs          = cell, directories where to save image file:
% .topoplot     = string, directory where to save topoplots.
% Tvalues       = Fieldtrip object with t-value for each
%               channel/frequency/time bin in .powspctrm
% iROI          = numeric scalar, index of selected ROI (to retrieve correct name).
% startTime     = numeric scalar, time from where to start plotting.
% endTime       = numeric scalar, time where to stop plotting.
% steps         = numeric scalar, time steps for which to make separate plot.
% zlim          = vector of two numerics, limit for color axis (-zlim zlim) in units of T.
% isSave        = Boolean, save plot (true) or not (false).
% 
% OUTPUTS:
% none, just create plot.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Adapted from Tobias Hauser (https://github.com/tuhauser/TAfT).
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/CueLockedAnalyses/TAfT/

%% Fill settings:

if nargin < 4
    iROI = 1;
    fprintf('No ROI specified -- use first ROI\n')
end

if nargin < 5
    if strcmp(job.lockSettings,'stimlocked')
        startTime   = 0.0;
        endTime     = 0.7;
        steps       = 0.1;
    elseif strcmp(job.lockSettings,'resplocked')
        startTime   = -1.0;
        endTime     = 0.4;
        steps       = 0.1;
    else
        error('Unknown lock settings\n')
    end    
    fprintf('No timings specified -- use default startTime = %.1f,endTime = %.1f, steps = %.1f\n',...
        startTime,endTime,steps)
end

if nargin < 8
    zlim = 3; % in units of T
    fprintf('No zlim specified -- use default zlim = %d\n',zlim)
end

if nargin < 9
    isSave = false;
    fprintf('Do not save by default\n')
end

%% Topoplot settings:

fprintf('ROI %s: Create topoplot over time\n',char(job.regNames(iROI)))

nCol            = 4; % number of columns in subplot
fontSize        = 18; 

if strcmp(job.lockSettings,'stimlocked')
    job.lockName= 'stim';
elseif strcmp(job.lockSettings,'resplocked')
    job.lockName= 'resp';
else
    error('Unknown lock settings\n')
end    

% Timing settings:
startTimeVec    = startTime:steps:endTime; % vector of start time for each subplot
endTimeVec      = (startTimeVec(1)+steps):steps:(startTimeVec(end)+steps);  % vector of end time for each subplot
nRow            = ceil(length(startTimeVec)/nCol); % number of rows following from bins

%% Create topoplot:

% Start figure:
figure('units','normalized','outerposition',[0 0 1 1]); hold on % fullscreen

for iPlot = 1:length(endTimeVec) % loop over given time windows

    subplot(nRow,ceil(length(endTimeVec)/nRow),iPlot)
    
    cfg = []; cfg.figure = gcf; cfg.zlim = [-1*zlim 1*zlim]; cfg.marker='on'; cfg.style='straight';
    cfg.layout = job.layout; cfg.comment = 'no'; cfg.xlim = [startTimeVec(iPlot) endTimeVec(iPlot)];
%    cfg.highlight = 'on'; cfg.highlightjob.channel = job.channels; cfg.highlightsymbol = 'o';
    cfg.colorbar = 'no'; % want 'no', i.e. do it yourself % --> add externally
    
    ft_topoplotER(cfg,Tvalues);
    
    title(sprintf('%.2f to %.2f sec', startTimeVec(iPlot),endTimeVec(iPlot)),'fontsize',fontSize) % 2 digits: 18
%     title(sprintf('%.1f to %.1f sec', startTimeVec(iPlot),endTimeVec(iPlot)),'fontsize',28) % 1 digit: 28

end % end iPlot

sgtitle(sprintf('Topoplot for ROI %s, %s-locked, time window %.02f - %.02f sec.',...
    char(job.regNames(iROI)),job.lockName,startTimeVec(1),endTimeVec(end)),'fontsize',fontSize)

fprintf('Finished topoplot\n');

%% Save:

if isSave
    saveas(gcf,fullfile(dirs.topoplot,sprintf('Topoplot_time_%s_%s_%s_%dsec_%s-%sms.png',...
        char(job.regNames(iROI)),job.lockSettings,job.HRFtype,job.trialdur,...
        num2str(1000*startTimeVec(1)),num2str(1000*endTimeVec(end)))))
end

end % end of function.
