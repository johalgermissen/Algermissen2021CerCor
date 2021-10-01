function singleTopoPlot(job,data,zlim)

% singleTopoPlot(job,data,zlim)
%
% Plots one single topoplot based on job.sigTime contasting two conditions.
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
%   .channels           = vector of strings, selected channel names.
%   .chanName           = string, name of selected channel area
%   capitalized.
%   .bandName           = string, name of selected frequency band
%   capitalized.
%   .twoLineLabels      = vector of 2 strings, labels for conditions.
% data                  = cell with the following fields:
%   .topoplot           = averaged over conditions/ subjects, for
%   topoplots.
% zlim                  = numeric, positive/ negative limit for for color 
% scale (default: 1).
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

if ~isfield(job,'sigTime')
    error('job.sigTime not specified');
end

%% Close any open figure window:

close all

%% Fixed settings:

fontSize    = 24;
lineWidth   = 5;


%% Start figure:

figure('name',sprintf('%s TF power Valence x Actual action',job.chanName),'units','normalized','outerposition',[0 0 1 1]); hold on % fullscreen

% Set config object for Fieldtrip:
cfg = []; cfg.zlim = [-1*zlim 1*zlim]; cfg.marker='on'; cfg.style='straight';
cfg.layout = job.layout; cfg.comment = 'no'; cfg.xlim = job.sigTime;
cfg.highlight = 'on'; cfg.highlightchannel = job.channels; cfg.highlightsymbol = 'o'; cfg.highlightsymbol = 'o'; cfg.highlightsize = 12; cfg.highlightfontsize = 12;
cfg.colorbar = 'yes'; % if 'no', then add colorbar yourself later
cfg.figure = gcf;

% Fieldtrip plot:
% Determine y limits based on input data, plot, add title:
if isfield(job,'freq')
    cfg.ylim = job.freq; 
    ft_topoplotTFR(cfg,data.topo2plot);
    title(sprintf('%s vs. %s:\n %s (%d-%d Hz) power from %.3f to %.3f sec', job.twoLineLabels{1},job.twoLineLabels{end}, job.bandName, job.freq(1), job.freq(end), job.sigTime(1),job.sigTime(end)),'fontsize',fontSize)

else
    ft_topoplotER(cfg,data.topo2plot);
    title(sprintf('%s vs. %s:\n voltage from %.3f to %.3f sec', job.twoLineLabels{1},job.twoLineLabels{end}, job.sigTime(1),job.sigTime(end)),'fontsize',fontSize)
end

%% Adapt linewidth of highligted disks:
% 1) ring, 2) left ear, 3) right ear, 4) nose, 5) line around, 6) small crosses, 
% 7) ring, 8) left ear, 9) right ear, 10) nose, 11) line around, 12) highlighted crosses,
% 13) left ear, 14) right ear, 15) nose, 16) line around
% --> loop from 12 to 16 to get all relevant objects

p = gca;
for i = 12:16
    t = p.Children(i); 
    t.LineWidth = lineWidth;
end

end % end of function.
