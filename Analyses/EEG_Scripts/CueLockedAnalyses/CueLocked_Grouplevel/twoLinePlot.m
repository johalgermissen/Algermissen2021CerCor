function twoLinePlot(job,data,errorbar,within)

% twoLinePlot(job,data,errorbar,within)
%
% Creates lineplot for two conditions specified as contrast in job.
% Works for both ERP and TFR data.
%
% INPUTS:
% job                   = cell with relevant fields for plotting settings, 
% required fields:
%   .validSubs          = numeric vector, subjects to include.
%   .nValidSubs         = numeric, number of subjects to include.
%   .lockSettings       = string, type of event-locking, 'stimlocked' or
%   'resplocked'.
%   .lockName           = string, name of event-locking, fully written out.
%   .twoLineLabels      = vector of 2 strings, labels for conditions.
%   .freq               = numeric vector of 2 elements, range of selected
%   frequencies.
%   .channels           = vector of strings, selected channel names.
%   .sigTime            = numeric vector of 2 elements, timing for
%   significant difference to highlight (optional).
% data                  = cell with the following fields:
%   .mat1 and mat2      = averaged over channels/ frequencies/ conditions.
%   .matMean            = mean of mat1 and mat2 (over all conditions).
%   .matGrandMean       = matMean averaged over subjects.
% errorbar              = add errorbars (true) or not (false) to plot
% (default: true):
% within                = perform Cousineau-Morey correct on error bars
% (true) or not (default: true)
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
    errorbar = true;
end

if nargin < 4
    within = true;
end

%% Close any open figure windows:

close all

%% Fixed settings:

fontSize        = 32; % 24; % 12 
fontSizeTitle   = 15; % keep overall title at 15
lineWidth1      = 5; % actual lines in plot
lineWidth2      = 3; % all other lines

%% Until when data are completely available:

if isfield('job','freq') 
    if strcmp(job.lockSettings,'stimlocked')
            until = 110; % for stimlocked;
    elseif strcmp(job.lockSettings,'resplocked')% resplocked:
            until = 80; % for stimlocked;
    else
        error('Unknown lockSettings')
    end
else
    until = length(data.mu.time);
end

%% Determine baseline (i.e. lowest line at time point zero):

iBaseline   = find(round(data.mu.time,3) == 0); % find baseline via time in sec
baseline    = min(nanmean(data.mat1(:,iBaseline)),nanmean(data.mat2(:,iBaseline)));

%% Plot:

p = {}; % delete p

figure('units','normalized','outerposition',[0 0 1 1]); hold on % fullscreen

% Lines:
% Use Cousineau 2005 method for uncorrected-within-subjects CIs 
% (only correction, no z; i.e. approx SE with only 2 conditions):

if errorbar % with error bars:

    if within % within-subjects (use Cousineau/Morey method): subtract individual mean, add group level mean, correct for #conds bz x/(x-1)
        p{1} = boundedline(data.mu.time(1:until),...
            nanmean(data.mat1(:,1:until))-baseline,...
            2/(2-1)*nanstd(data.mat1(:,1:until)-data.matMean(:,1:until)+repmat(data.matGrandMean(1:until),job.nValidSubs,1))./sqrt(job.nValidSubs),...
            'cmap',job.twoLineColor(1,:),'alpha','transparency',0.2); %
        set(p{1},'linewidth',lineWidth1) % adjust linewidth (not possible within boundedline)
        p{2} = boundedline(data.mu.time(1:until),...
            nanmean(data.mat2(:,1:until))-baseline,...
            2/(2-1)*nanstd(data.mat2(:,1:until)-data.matMean(:,1:until)+repmat(data.matGrandMean(1:until),job.nValidSubs,1))./sqrt(job.nValidSubs),...
            'cmap',job.twoLineColor(2,:),'alpha','transparency',0.2); %
        set(p{2},'linewidth',lineWidth1) % adjust linewidth (not possible within boundedline)

    else % between-subjects (standard way of compute SEs):
        p{1} = boundedline(data.mu.time(1:until),...
            nanmean(data.mat1(:,1:until)),nanstd(data.mat1(:,1:until))./sqrt(length(job.validSubs)),...
            'cmap',job.twoLineColor(1,:),'alpha','transparency',0.2); %
        p{2} = boundedline(data.mu.time(1:until),...
            nanmean(data.mat2(:,1:until)),nanstd(data.mat2(:,1:until))./sqrt(length(job.validSubs)),...
            'cmap',job.twoLineColor(2,:),'alpha','transparency',0.2); %
    end
    
else % without error bars:
    p{1} = plot(data.mu.time,squeeze(nanmean(data.mat1)),'-','linewidth',lineWidth1,'Color',job.twoLineColor(1,:)); % data.mat1: blue
    p{2} = plot(data.mu.time,squeeze(nanmean(data.mat2)),'-','linewidth',lineWidth1,'Color',job.twoLineColor(2,:)); % data.mat2: red
end

%% y-axis labels:

xlabel('Time (s)','fontsize',fontSize,'fontweight','bold')

%% Retrieve y-axis limits, change later:

yLim    = get(gca,'ylim'); % retrieve automatic settings, change later
yMinLim = yLim(1); 
yMaxLim = yLim(2); 
if isfield(job,'freq'); yMinLim = yMinLim - 0.1; end

%% Add x-axis settings and labels, add vertical lines:

if strcmp(job.lockSettings,'stimlocked')
    set(gca,'xlim',[-0.25 1.3],'ylim',[yMinLim yMaxLim],....
        'xtick',[-1 -0.5 0 0.5 .815 1 1.3],'xtickLabel',{'-1','-0.5','Cue','0.5','AvgRT','1','+'},...
        'fontsize',fontSize,'linewidth',lineWidth2)
    % Vertical lines:
    plot([1.3 1.3], get(gca,'ylim'),'k') % cue offset
    plot([.815 .815], get(gca,'ylim'),':k') % average response time

elseif strcmp(job.lockSettings,'resplocked')
    set(gca,'xlim',[-1 1],'ylim',[yMinLim yMaxLim],'fontsize',fontSize);
    % Vertical lines:
    plot([0 0],get(gca,'ylim'),':k');
end

%% Add title and y-label:

if isfield(job,'freq')
    
    title(sprintf('%s vs. %s: Time-frequency power %d-%d Hz over %s (%s)',...
        job.twoLineLabels{1},job.twoLineLabels{end},job.freq(1),job.freq(2),strjoin(job.channels,'/'),job.lockName),...
        'fontsize',fontSizeTitle)
    ylabel('Power (dB)','fontsize',fontSize,'fontweight','bold')
    yMinLim = yMinLim - 0.1; % adjust a bit
    yMaxLim = yMaxLim + 0.1; % adjust a bit
    set(gca,'ylim',[yMinLim yMaxLim]) % set anew
    boost   = 0.05; % boost: negative offset from yMaxLim to draw significance lines
    
else  
  
    if isfield(job,'ROI2use')
        title(sprintf('Split per %s BOLD: Voltage over %s (%s)',job.type.ROI2use,strjoin(job.channels,'/'),job.lockName),...
	'fontsize',fontSizeTitle)
    
    else
        title(sprintf('%s vs. %s: Voltage over %s (%s)',...
	job.twoLineLabels{1},job.twoLineLabels{end},strjoin(job.channels,'/'),job.lockName),...
	'fontsize',fontSizeTitle)
    end
    
    ylabel('Amplitude (A/cm²)','fontweight','bold','fontsize',32); 
    boost   = 0.0015; % boost: negative offset from yMaxLim to draw significance lines

end

%% Add area of significance:

legendNames = job.twoLineLabels; % initialize labels

if isfield(job,'sigTime') && ~isempty(job.sigTime) % add line of significance
    
    nSig = length(job.sigTime)/2; % job.sigTime allows for pairs of number from where to where difference is significant

    for iSig = 1:nSig % loop over pairs
        
        fprintf('Found significant timing at %.03f - %.03f sec\n',job.sigTime(iSig*2-1),job.sigTime(iSig*2));
%         fprintf('Plotted at %0.3f\n',yMaxLim-boost); % give y-axis coordinate of added horizontal line
        p{end+1} = plot([job.sigTime(iSig*2-1) job.sigTime(iSig*2)],[yMaxLim-boost yMaxLim-boost],'k','linewidth',lineWidth1);

    end
    legendNames{end+1} = 'Significant difference'; % update if necessary
end

%% Add legend:

legend([p{:}],legendNames,'fontsize',fontSize); legend boxoff

end % end of function.
