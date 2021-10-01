function TFplot(job,data,zlim)

% TFplot(job,data,zlim)
%
% Creates time-frequency plot (2D heat map) contasting two conditions.
%
% INPUTS:
% job                   = cell with relevant fields for plotting settings, 
% required fields:
%   .TFtiming           = vector of 2 scalars, start and end timing of TF
%   plot (determined usually by lockSettings).
%   .channels           = vector of strings, selected channel names.
%   .chanName           = string, name of selected channel area
%   capitalized.
%   .twoLineLabels      = vector of 2 strings, labels for conditions.
% data                  = cell with the following fields:
%   .TF2plot            = averaged over channels/ conditions/ subjects, for
%   time-frequency plots.
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

%% Close any open figure window:

close all

%% Fixed settings:

fontSize    = 24;
lineWidth   = 3;

%% Plot figure:

figure('name',sprintf('%s vs. %s: Time-frequency power over %s (%s)',job.twoLineLabels{1},job.twoLineLabels{end},strjoin(job.channels,'/'),job.lockName),'units','normalized','outerposition',[0 0 1 1]); hold on % fullscreen

% Contour plot:
contourf(data.mu.time,data.mu.freq,data.TF2plot,40,'linestyle','none'); hold on

% Vertical lines:
plot([0 0],get(gca,'ylim'),':k','LineWidth',lineWidth); % cue onset
plot([1.3 1.3],get(gca,'ylim'),':k','LineWidth',lineWidth); % cue offset

% Settings:
set(gca,'xlim',[job.TFtiming],'ylim',[1 15],'clim',[-1*zlim 1*zlim],...
    'xtick',[-0.5 -0.25 0 0.25 .5 0.75 1],... 
    'ytick',2:2:14,'yscale','lin',...
    'fontsize',fontSize,'Linewidth',lineWidth) %

box off

% Labels:
ylabel('Frequency (Hz)','fontsize',fontSize);
xlabel('Time (s)','fontsize',fontSize);

% Title:
title(sprintf('%s vs. %s: Time-frequency power over %s (%s)',job.twoLineLabels{1},job.twoLineLabels{end},strjoin(job.channels,'/'),job.lockName),'fontsize',15)

% Colorbar:
colorbar

end % end of function.
