function EEGfMRIPav_CueLocked_3_visualICs(rootdir,subVec)

% EEGfMRIPav_CueLocked_3_visualICs(rootdir,subVec)
% 
% Plot independent components obtained from ICA.
% 
% INPUTS:
% rootdir           = string, root directory of project
% subVec            = vector of integers, numbers of subjects to process
%
% OUTPUTS:
% saves plots of each component to disk under dirs.ICAplots.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% By J.C.Swart, 2016/ J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Preprocessing/

%% Complete input settings:

% clear all; close all; clc
% dbstop if error

if ~exist('rootdir','var')
    rootdir =  preprocessing_set_rootdir(); % '/project/3017042.02';
    fprintf('rootdir unspecified, assume %s\n',rootdir)
end

if ~exist('subVec','var')
    subVec = 1:36; 
    fprintf('subVec unspecified, assume %s\n',strjoin(string(subVec),', '));
end

%% Set directories:

% Add helper files to path:
fprintf('Add helper files to path\n');
addpath(fullfile(rootdir,'Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Helpers'));

% Set directories and parameters:
dirs    = set_dirs(rootdir);
par     = set_par();

%% Detect EEG raw data directories:

% Detect all available input folders (all subjects):
tmp         = dir(fullfile(dirs.ICA,'*.mat'));
fileList    = {tmp.name};
nInput      = length(fileList);
fprintf('Found data files from %d subjects\n',nInput);

% Extract subject numbers from folders available:
subNames    = str2double(extractBetween(fileList,8,10));

% Extract indices of selected subjects:
selIndices  = find(ismember(subNames,subVec));
nSub        = length(selIndices);

%% Loop over subjects:

for iSub = selIndices % iSub = 1;

    %% Create input file name:
    
    inputFile       = fullfile(dirs.ICA, fileList{iSub}); % go to subject directory
    fprintf('EEG input file is %s\n',inputFile);

    %% Load input data:
    
    fprintf('Subject %03d: Start loading ICs\n',iSub)
    load(fullfile(dirs.ICA,inputFile));

    % Provide subject ID on screen:
    fprintf('Subject %03d: Finished loading ICs\n',iSub)
    
    %% Retrieve component activations:
    
    compActivation = cat(3,comp.trial{:});
    nTrial 	   = length(com.trial);

    %% Compute frequency spectrum of components (adapted from ft_icabrowser.m). 
    
    % Takes a while -> 1 minute or so for 61 components :)
    fprintf('Subject %03d: Compute frequency spectrum of each component\n',iSub);
    fft_data = cat(2,comp.trial{1:5:end});
    smo = 50; % smoothing 
    steps = 10; % steps
    Fs = comp.fsample; % sampling rate, should be 1000 Hz --> for Nyquist frequency
    N = floor(size(fft_data,2)); % number of samples --> for Rayleigh frequency
    
    %% Loop over components:
    
    nComp       = numel(comp.label); % number of components
    smoothed    = cell(nComp,1); % initialize
    freq        = cell(nComp,1); % initialize
    strt        = cell(nComp,1); % initialize
    stp         = cell(nComp,1); % initialize

    for iComp = 1:nComp % Do FFT on each component to plot power spectrum
        
        fprintf('Subject %03d Component %02d: Start FFT\n',iSub,iComp);
        
        % Perform FFT:
        xdft = fft(fft_data(iComp,:)); % FFT on this component
        xdft = xdft(1:N/2+1); % select frequencies: 0 till Nyquist frequency
        psdx = (1/(Fs*N)).*abs(xdft).^2; % power or mean squares = sum of squares divided by length (number samples times sampling rate)
        psdx(2:end-1) = 2*psdx(2:end-1); % times 2 (except for first and last frequency)

        % Smooth components:
        j = 1; k = 1; % j = frequency; k = smoothing strength (size of moving window)
        while j < length(psdx)-smo % while j low enough
            smoothed{iComp}(k)=mean(psdx(j:j+smo)); % mean of moving window
            j = j + steps; k = k + 1;
        end

        % Template for frequencies to-be-displayed:
        freq{iComp} = linspace(0,Fs/2,size(smoothed{iComp},2)); % create template of evenly spaced frequencies
        strt{iComp} = find(freq{iComp} > 2,1,'first'); % index of first frequency above 2 Hz
        stp{iComp}  = find(freq{iComp} < 500,1,'last'); % index of last frequency below 500 Hz (Nyqvist)

    end

    %% Plot topoplots of several components (20 components per plot):
    
    fprintf('Subject %03d: Create overviews of topoplots of components\n',iSub); 

    % a) Components 1-20:
    figure('Position',[200 200 800 800],'visible','off');
    cfg                     = []
    cfg.figure              = gcf;
    cfg.layout              = par.chan.layoutfile;
    cfg.component           = 1:min(20,nComp);
    ft_topoplotIC(cfg,comp);
    saveas(gcf,fullfile(dirs.ICAplots,sprintf('MRIpav_%03d_ICs_allTopoplots_01-20.png',iSub)));
    
    % b) Components 21-40:
    if nComp > 20
        figure('Position',[200 200 800 800],'visible','off');
        cfg.component           = 21:min(40,nComp);
        cfg.figure              = gcf;
        ft_topoplotIC(cfg,comp);
        saveas(gcf,fullfile(dirs.ICAplots,sprintf('MRIpav_%03d_ICs_allTopoplots_21-40.png',iSub)));
    end
    
    % c) Components 41-60:
    if nComp > 40
        figure('Position',[200 200 800 800],'visible','off');
        cfg.component           = 41:min(60,nComp);
        cfg.figure              = gcf;
        ft_topoplotIC(cfg,comp);
        saveas(gcf,fullfile(dirs.ICAplots,sprintf('MRIpav_%03d_ICs_allTopoplots_41-60.png',iSub)));
    end
    
    % d) Components 61-80:
    if nComp > 60
        figure('Position',[200 200 800 800],'visible','off');
        cfg.component           = 61:min(80,nComp);
        cfg.figure              = gcf;
        ft_topoplotIC(cfg,comp);
        saveas(gcf,fullfile(dirs.ICAplots,sprintf('MRIpav_%03d_ICs_allTopoplots_60Plus.png',iSub)));
    end
    
    %% Plot for each component: topoplot, trial-by-time activation and power spectrum.
    
    fprintf('Subject %03d: Create overview plots for each component\n',iSub)

    for iComp = 1:nComp
        
        fprintf('Subject %03d Component %02d: Start plotting\n',iSub,iComp)
        
        % 1) Topoplot:
        figure('Position',[100 600 1800 300],'visible','off');
        subplot(1,5,1)
        cfg                     = [];
        cfg.figure              = gcf;
        cfg.layout              = par.chan.layoutfile;
        cfg.component           = iComp;
        ft_topoplotIC(cfg,comp);
        title(num2str(iComp,'Topoplot Component No. %02d')); 
        
        % 2) Trial-by-time activation (time x trial matrix):
        subplot(1,5,2)
        imagesc(comp.time{1},1:nTrial,squeeze(compActivation(iComp,:,:))')
        set(gca,'clim',[-.5*mean(abs(get(gca,'clim'))) .5*mean(abs(get(gca,'clim')))],...
            'ytick',0:50:nTrial,'ydir','normal')
        xlabel('Time (s)'); ylabel('Trial')
        title('Amplitude per trial per time point');
        
        % 3) Event-related activation (mean of all trials):
        subplot(1,5,3)
        plot(comp.time{1},squeeze(mean(compActivation(iComp,:,:),3)))
        set(gca,'xlim',comp.time{1}([1 end])); xlabel('Time (s)'); 
        title('Mean time course during trial');
        
        % 4) Activation of selected component for all trials (all trials in one plot superimposed):
        subplot(1,5,4); hold on
        cmap = [colormap;colormap];
        for iTrial = 1:nTrial
            plot(comp.time{1},squeeze(compActivation(iComp,:,iTrial)),'Color',cmap(iTrial,:));
        end
        set(gca,'xlim',comp.time{1}([1 end])); xlabel('Time (s)'); 
        title('Time course for each trial');
        
        % 5) Power spectrum:
        subplot(1,5,5)
        plot(freq{iComp}(strt{iComp}:stp{iComp}),...
            log10(smoothed{iComp}(strt{iComp}:stp{iComp})));
        set(gca,'TickDir','out','XTick',0:10:50,'xlim',[1 50])
        xlabel('Frequency (Hz)'); ylabel('(dB/Hz)');
        title('Power spectrum');

        % Save if necessary:
        set(gcf, 'PaperUnits', 'centimeters');
        set(gcf,'PaperPosition',[0 0 30 5]);
        saveas(gcf,fullfile(dirs.ICAplots,sprintf('MRIpav_%03d_ICs_Comp_%02d.png',iSub,iComp)));
       
    end % end iComp-loop.

    fprintf('Subject %03d: Finished plots\n',iSub)
    clear comp compActivation fft_data psdx stp strt xdft tmp freq
    close all
    fprintf('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% \n');
    
end % end iSub-loop.

fprintf('Done :-) \n');

end % end of function.
