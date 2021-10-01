function multiTopoPlot(job,data,zlim,nRows,startTime,endTime,steps)

% multiTopoPlot(job,data,zlim,nRows,startTime,endTime,steps)
%
% Plots multiple topoplots based on start and end time in steps contrasting
% two conditions.
% Works for both ERP and TFR data.
%
% INPUTS:
% job                   = cell with relevant fields for plotting settings, 
% required fields:
%   .layout             = string, cap to use for topoplots.
%   .lockSettings       = string, type of event-locking, 'stimlocked' or
%   'resplocked'.
%   .freq               = numeric vector of 2 elements, range of selected
%   frequencies.
%   .bandName           = string, name of selected frequency band
%   capitalized.
% data                  = cell with the following fields:
%   .topoplot           = averaged over conditions/ subjects, for
%   topoplots.
%
% zlim                  = numeric, positive/ negative limit for for color 
% scale (default: 1).
% nRows                 = number of rows of grid of plots (default: 2).
% startTime             = start time of first plot (default: based on
% job.lockSettings).
% endTime               = start time of last plot (default: based on
% job.lockSettings).
% steps                 = time steps between each plot (default: based on
% job.lockSettings).
%
% OUTPUTS:
% none, just plotting.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel/

%% Complete settings:

if nargin < 3
    zlim = 1;
end

if nargin < 4
    nRows = 2;
end

if nargin < 7
    if nargin > 4
        fprintf('Incomplete timing settings, will use defaults')
    end
    if strcmp(job.lockSettings,'stimlocked')
        startTime   = 0.0;
        endTime     = 1.3;
        steps       = 0.1;
    elseif strcmp(job.lockSettings,'resplocked')
        startTime   = -0.5;
        endTime     = 0.4;
        steps       = 0.1;
    else
        error('Invalid response settings')
    end
end

%% Close any open figure window:

close all

%% Fixed settings:

fontSize    = 12; % smaller because smaller panels

%% Timing settings:

startTimeVec = startTime:steps:endTime;
endTimeVec = (startTimeVec(1)+steps):steps:(startTimeVec(end)+steps);

%% Plot:

figure('units','normalized','outerposition',[0 0 1 1]); hold on % fullscreen

for iPlot = 1:length(endTimeVec)

    % Start subplot:
    subplot(nRows,ceil(length(endTimeVec)/nRows),iPlot)
    
    % Set config object for Fieldtrip:
    cfg = []; cfg.zlim = [-1*zlim 1*zlim]; cfg.marker='on'; cfg.style='straight';
    cfg.layout = job.layout; cfg.comment = 'no'; cfg.xlim = [startTimeVec(iPlot) endTimeVec(iPlot)];
%    cfg.highlight = 'on'; cfg.highlightchannel = job.channels; cfg.highlightsymbol = 'o'; % individual plots too small for highlights
    cfg.colorbar = 'no'; % want 'no', i.e. do it yourself % --> add externally
    cfg.figure = gcf;    
    
    % Determine y limits based on input data, plot, add title:
    if isfield(job,'freq') % frequency data
        cfg.ylim = job.freq; 
        p = ft_topoplotTFR(cfg,data.topo2plot);
        title(sprintf('%s %.2f to %.2f sec', job.bandName, startTimeVec(iPlot), endTimeVec(iPlot)), 'fontsize', fontSize)

    else % time domain data
        p = ft_topoplotER(cfg,data.topo2plot);
        title(sprintf('%.2f to %.2f sec', startTimeVec(iPlot),endTimeVec(iPlot)),'fontsize',fontSize)
    end
    

end

end
% END
