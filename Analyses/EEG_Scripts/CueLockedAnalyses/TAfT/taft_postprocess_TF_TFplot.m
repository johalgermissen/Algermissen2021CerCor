function [tg, corrp] = taft_postprocess_TF_TFplot(job,dirs,sortBetas,iROI,selChans,thresh,nP,zlim,isSave)

% [tg, corrp] = taft_postprocess_TF_TFplot(job,dirs,sortBetas,iROI,selChans,thresh,nP,zlim,isSave)
% 
% For fMRI-EEG regression weights for selected ROI, 
% select data for given channels/ time range, perform cluster-based 
% permutation test for mean signal across selected channels, output p-value,
% create TF-plot of T-values with with "significant clusters" highlighted.
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
% nP            = numeric scalar, number of permutations in cluster-based 
%               permutation test, default is 10000.
% zlim          = numeric scalar, limit for color axis (-zlim zlim) in units of T,
%               default is thresh.
% isSave        = Boolean, save plot (true) or not (false).
% 
% OUTPUTS:
% tg            = 2D-matrix with T-value for each time/frequency bin
%               (averaged over channels).
% corrp         = 2D-matrix with p-value for respective cluster to which
%               time/frequency bin belongs.
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
    fprintf('No ROI specified -- use ROI %d\n',iROI)
end

if nargin < 5
    selChans = {'FCz','Cz'}; % Go/NoGo contrast
    fprintf('No channels selected -- use %s\n',strjoin(selChans,'/'));
end

if nargin < 6
    thresh = 2.0;
    fprintf('No cluster-forming threshold specified -- use %.1f\n',thresh)
end

if nargin < 7
    nP = 10000;
    fprintf('No number of permutations specified -- use %d\n',nP)
end

if nargin < 8
    zlim = thresh;
    fprintf('Use threshold of permutation test as zlim for plotting\n')
end

if nargin < 9
    isSave = false;
    fprintf('Save plot by default\n')
end

%% Fixed settings:

fontSize    = 24;
pCrit       = 0.05;

%% Select channels, average over channels:

chanBetas = cell(length(sortBetas),1); % initialize objects for averaged channels (note: can have different time/ frequency axes!)

c   = []; % initialize empty 4D object

cfg = []; % empty cfg file for later data selection

% a) Select channels:
cfg.channel     = selChans;
cfg.avgoverchan = 'yes';

% b) Select latency:
if strcmp(job.lockSettings,'stimlocked')
    cfg.latency     = [0 1.3];
    job.lockName    = 'stim';
%     cfg.latency = [0 2]; % for extra long analyses
elseif strcmp(job.lockSettings,'resplocked')
    cfg.latency     = [-1.0 0.5]; % maximum
    job.lockName    = 'resp';
else
    error('Unknown lock settings')
end
cfg.avgovertime = 'no';

% c) Select frequency:
% Default: don't average over frequency
% cfg.frequency = [1 15];
% cfg.frequency = [4 8];
% cfg.frequency = [1 8];
% cfg.avgoverfreq = 'no';

fprintf('ROI %s: Average over channels %s, reshape\n',char(job.regNames(iROI)),strjoin(selChans,'/'));

%% Reshape to 4D object:

for iSub = 1:length(sortBetas) % iSub = 1;

    chanBetas{iSub} = ft_selectdata(cfg,sortBetas{iSub}); % average over selected channels
    c(iSub,:,:,:)   = chanBetas{iSub}.powspctrm; % bring into 4D, with subject as first dimension

end

%% Cluster-based permutation test:

[corrp,tg] = clustertf(c,thresh,nP); % default test statistic
% [corrp,tg] = clustertf(c,thresh,nP,'nExtent');
% [corrp,tg] = clustertf(c,thresh,nP,'maxT');
% [corrp,tg] = clustertf(c,thresh,nP,'sumT');

%% Two subplots: heat map of T-values, heat map of p-values:

figure('units','normalized','outerposition',[0 0 1 1]); hold on;

% First subplot: heat map of t-values:
subplot(1,2,1);
contourf(chanBetas{1}.time,chanBetas{1}.freq,real(squeeze(tg)),30,'linestyle','none');colorbar
% set(gca,'yscale','log','ylim',[1.5 15],'ytick',[2 4 8 16],'clim',[-1*zlim 1*zlim],'fontsize',fontSize);
set(gca,'clim',[-1*zlim 1*zlim],'fontsize',fontSize);
xlabel('Time (in s)','fontweight','bold','fontsize',fontSize); 
ylabel('Frequency (in Hz)','fontweight','bold','fontsize',fontSize)
title('T-values','fontsize',fontSize)

if sum(corrp(:) < pCrit) > 0 % if anything significant
    
   % Plot overlay:
   % see https://nl.mathworks.com/matlabcentral/answers/250279-how-to-make-contour-plots-transparent-in-matlab-r2015a#answer_211204
   
   pause(1)
      
   % Locate significant clusters:
   isSig = squeeze(double(corrp < pCrit));
   isSigTimeIdx = find(sum(isSig,1)>0);
   isSigFreqIdx = find(sum(isSig,2)>0);
   fprintf('significant cluster from %.03f - %.03f sec., %.03f - %.03f Hz\n',...
       chanBetas{1}.time(min(isSigTimeIdx)),chanBetas{1}.time(max(isSigTimeIdx)),...
       chanBetas{1}.freq(min(isSigFreqIdx)),chanBetas{1}.freq(max(isSigFreqIdx)));

   % New contour plot highlighting "significant clusters":
   hold on % on top of old plot
   [~, hContour]  = contourf(chanBetas{1}.time,chanBetas{1}.freq,squeeze(isSig),1); 
   hContour.LineWidth = 5;
   drawnow;  % this is important, to ensure that FacePrims is ready in the next line!
   hFills = hContour.FacePrims;  % array of TriangleStrip objects
   [hFills.ColorType] = deal('truecoloralpha');  % default = 'truecolor'
   for idx = 1:numel(hFills)
       hFills(idx).ColorData(4) = 1;   % default=255
   end
   hold off
end

% Second subplot: Heat map of p-values:
subplot(1,2,2);
if sum(corrp(:)<1) > 1
    contourf(chanBetas{1}.time,chanBetas{1}.freq,squeeze(corrp),'linestyle','none');colorbar
    set(gca,'clim',[0 1],'fontsize',fontSize)
    xlabel('Time (in s)','fontweight','bold','fontsize',fontSize); 
    ylabel('Frequency (in Hz)','fontweight','bold','fontsize',fontSize)
    title('p-values')
end

% Overall title:
sgtitle(sprintf('TF-plot for ROI %s, %s-locked, GLM-type %s, channels %s',...
    char(job.regNames(iROI)),job.lockName,job.HRFtype,strjoin(selChans,'/')),'fontsize',fontSize)

fprintf('Finished TF plot\n');

%% Save figure:

if isSave 
    saveas(gcf,fullfile(dirs.TFplot,sprintf('TFplot_%s_%s_%s_%dsec_%s.png',char(job.regNames(iROI)),job.lockSettings,job.HRFtype,job.trialdur,strjoin(selChans,''))))
end

end % end of function.
