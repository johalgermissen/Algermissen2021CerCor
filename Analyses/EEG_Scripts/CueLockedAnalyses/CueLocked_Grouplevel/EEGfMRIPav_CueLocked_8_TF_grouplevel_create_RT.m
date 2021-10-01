function EEGfMRIPav_CueLocked_8_TF_grouplevel_create_RT(rootdir,subVec,lockSettings,baselineSettings,baselineTimings,...
    accSettings,TFtype,nFreq,nBin,stimERPcor,respERPcor)

% EEGfMRIPav_CueLocked_8_TF_grouplevel_create_RT(rootdir,subVec,lockSettings,baselineSettings,baselineTimings,...
%     accSettings,TFtype,nFreq,ROI2use,nBin,stimERPcor,respERPcor)
%
% Load EEG, rejected trials, behavior; apply baseline correction; sort
% trials into bins based on RTs and cue valence, aggregate within bins,
% save. 
%
% INPUTS:
% rootdir           = string, root directory of project.
% subVec            = vector of integers, numbers of subjects to include
% (default: 1:36).
% lockSettings      = string, type of time-locking before TF decomposition,
% either 'stimlocked' (default) or 'resplocked'.
% baselineSettings  = string, type of baseline correction, 'trend'
% (regression-based), 'all' (grand mean oer subject), 'condition'
% (grand mean separately per condition), 'trial' (per trial, own
% implementation), 'ft_trial' (with Fieldtrip's ft_freqbaseline).
% baselineTimings   = vector of 1 or 2, timing of baseline correction
% for which to load data.
% accSettings       = string, whether certain trials are selected based on
% accuracy,either 'correct', 'incorrect', 'noAcc' (no distinction).
% TFtype            = type of TF decomposition, 'hanning' or 'morlet'.
% nFreq             = numeric, number of frequency bins decomposed
% (default: 15).
% nBin              = integer, number of bins to create based on BOLD
% signal (default: 3).
% stimERPcor        = Boolean, read data corrected for stimulus-locked 
% ERP (true) or not (false), default false.
% respERPcor        = Boolean, read data corrected for response-locked 
% ERP (true) or not (false), default false.
%
% OUTPUTS:
% saves data to disk under dirs.TFgroup.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% By J.C.Swart, 2016/ J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel/

%% Complete input settings:

% clear all; close all; clc
% dbstop if error

if ~exist('rootdir','var')
    rootdir = grouplevel_set_rootdir(); % '/project/3017042.02';
    fprintf('rootdir unspecified, assume %s\n',rootdir)
end

if ~exist('subVec','var')
    subVec = 1:36; 
    fprintf('subVec unspecified, assume %s\n',strjoin(string(subVec),', '));
end

if ~exist('lockSettings','var')
    lockSettings = 'stimlocked'; % stimlocked resplocked
    fprintf('lockSettings unspecified, assume %s\n',lockSettings);
end

if ~exist('baselineSettings','var')
    baselineSettings = 'trend'; % trend, all, condition, trial, ft_trial
    fprintf('baselineSettings unspecified, assume %s\n',baselineSettings);
end

if ~exist('baselineTimings','var')
    baselineTimings = [-0.250 -0.050]; % stimlocked resplocked
    fprintf('baselineTimings unspecified, assume %s\n',strjoin(string(baselineTimings),' - '));
end

if ~exist('accSettings','var')
    accSettings = 'correct'; % correct incorrect noACC
    fprintf('accSettings unspecified, assume %s\n',accSettings);
end

if ~exist('TFtype','var')
    TFtype = 'hanning'; % hanning, morlet
    fprintf('TFtype unspecified, assume %s\n',TFtype);
end

if ~exist('nFreq','var')
    nFreq = 15; 
    fprintf('nFreq unspecified, assume %d\n',nFreq);
end

if ~exist('nBin','var')
    nBin = 3;
    fprintf('nBin unspecified, assume %d\n',nBin);
end

if ~exist('stimERPcor','var')
    stimERPcor = false;
end
    
if ~exist('respERPcor','var')
    respERPcor = false;
end

%% Set directories:

% Add helper files to path:
fprintf('Add helper files to path\n');
addpath(fullfile(rootdir,'Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Helpers'));

% Set directories and parameters:
dirs        = set_dirs(rootdir);
par         = set_par();

% Add TF directory based on lockSettings:
dirs.TF     = fullfile(dirs.results,sprintf('TF_subjectlevel_%s_%s%d',lockSettings,par.TF.TFtype,par.TF.nFreq));

if stimERPcor 
    fprintf('Use data with stimulus-locked ERP subtracted before TF decomposition \n');
    dirs.TF = sprintf('%s_stimERPcor',dirs.TF); % overwrite target file name
end
if respERPcor
    fprintf('Use data with response-locked ERP subtracted before TF decomposition \n');
    dirs.TF = sprintf('%s_respERPcor',dirs.TF); % overwrite target file name
end

if strcmp(accSettings,'noACC')
    allAcc = [1 0];

elseif strcmp(accSettings,'correct')
    allAcc = 1;

elseif strcmp(accSettings,'incorrect')
    allAcc = 0;

else
    fprintf('Unknown accuracy setting')
end

%% Detect input files:

% Detect all available input folders (all subjects):
tmp         = dir(fullfile(dirs.TF,sprintf('*_TF_%s.mat',lockSettings)));
fileList    = {tmp.name};
nInput      = length(fileList);
fprintf('Found data files from %d subjects\n',nInput);

% Extract subject numbers from folders available:
subNames    = str2double(extractBetween(fileList,8,10));

% Extract indices of selected subjects:
selIndices  = find(ismember(subNames,subVec));
nSub        = length(selIndices);
nCond       = nBin*2; % number of bins, separately for each valence

%% Create final file name:

outputFile = sprintf('EEGfMRIPav_TF_%s_baseCor_%s_baseTime_%03dms_%03dms_%s_%s%d_RT%d', ...
    lockSettings,baselineSettings,baselineTimings(1)*-1000,baselineTimings(end)*-1000,...
    accSettings,TFtype,nFreq,nBin);

if stimERPcor
    outputFile = [outputFile '_stimERPcor'];
end

if respERPcor
    outputFile = [outputFile '_respERPcor'];
end

%% Create final file name:

outputFile = fullfile(dirs.TFgroup,sprintf('EEGfMRIPav_TF_%s_baseCor_%s_baseTime_%03dms_%03dms_%s%d_RT%d.mat', ...
    lockSettings,baselineSettings,baselineTimings(1)*-1000,baselineTimings(end)*-1000,TFtype,nFreq,nBin));

fprintf('Selected file is %s\n',outputFile)

%% Loop over subjects:

if ~exist(outputFile,'file')

    subCount = 1; % relative subject number, i.e. also works if certain numbers are missing

    TFall = cell(nSub,1);
    nTrialCond = nan(nSub, nCond);
    
    for iSub = selIndices % %iSub = 1;
                
        %% 1) Load EEG data:
        
        fprintf('Subject %03d: Start loading TF data\n', iSub)
        load(fullfile(dirs.TF,fileList{iSub}))
        fprintf('Subject %03d: Finished loading TF data\n', iSub)
        
        %% 2) Load rejected trials:

        fprintf('Subject %03d: Load rejected trials\n',iSub)
        filenameRejTrials = fullfile(dirs.rejTrials,sprintf('MRIpav_%03d_rejTrials.mat',iSub));
        load(filenameRejTrials);

        %% 3) Load and recode behavior:
        
        fprintf('Subject %03d: Load behavioral data\n',iSub)
        behav = load_recode_behav(rootdir,iSub,rejectedtrials); 
        
        %% 4) Sort RT quantiles:
        
        % RT quantiles:
        clear X Xsel G2Widx G2Aidx G2Wval G2Aval G2Wcutoff G2Acutoff G2WbinIdx G2AbinIdx
        fprintf('Subject %03d: Sort RTs into %d bins per G2W/ G2A \n',iSub, nBin) % start recoding

        X = behav.RT; % retrieve RTs
        G2Widx = (behav.val == 1 & behav.isgo == 1 & ismember(behav.accuracy,allAcc)); % indices of G2W trials
        G2Aidx = (behav.val == 2 & behav.isgo == 1 & ismember(behav.accuracy,allAcc)); % indices of G2A trials
        G2Wval = X(G2Widx); % G2W RT values
        G2Aval = X(G2Aidx); % G2A RT values
        
        % Cutoffs:
        G2Wcutoff = nan(1,nBin); 
        G2Acutoff = nan(1,nBin);
        for iBin = 1:nBin
                G2Wcutoff(iBin) = quantile(G2Wval,iBin/nBin); % lowest third cutoff
                G2Acutoff(iBin) = quantile(G2Aval,iBin/nBin); % lowest third cutoff
        end
        
        % Initialize logical indices of quantiles:
        G2WbinIdx = nan(length(X),nBin); 
        G2AbinIdx = nan(length(X),nBin);
        
        % First and last quantile:
        G2WbinIdx(:,1) = (G2Widx & X <= G2Wcutoff(1)); % lowest indices
        G2AbinIdx(:,1) = (G2Aidx & X <= G2Acutoff(1)); % lowest indices
        G2WbinIdx(:,end) = (G2Widx & X > G2Wcutoff(end-1)); % highest indices
        G2AbinIdx(:,end) = (G2Aidx & X > G2Acutoff(end-1)); % highest indices
        
        % Middle quantiles:
        if nBin > 2 
            for iBin = 2:(nBin-1)
                G2WbinIdx(:,iBin) = (G2Widx & X > G2Wcutoff(iBin-1) & X <= G2Wcutoff(iBin));
                G2AbinIdx(:,iBin) = (G2Aidx & X > G2Acutoff(iBin-1) & X <= G2Acutoff(iBin));
            end
        end
        
        % Create decision variable:
        Xsel = G2WbinIdx * [1:nBin]' + G2AbinIdx * [(nBin+1):(2*nBin)]';
        % tabulate(Xsel)
        
        if max(Xsel) > nBin*2
            error('Too many conditions, check Xsel')
        end        
        
        %% 5) Baseline-correction per subject (and transformation into dezibel). 

        % Retrieve baseline timings:
        timeIdx = dsearchn(freq.time',baselineTimings'); % select indices for baseline time window
        
        if length(timeIdx) > 1; timeIdx = timeIdx(1):timeIdx(end); end % extend to interval if start and end given
        fprintf('Subject %03d: Compute baseline, setting %s\n', iSub, baselineSettings);
        
        if strcmp(baselineSettings,'all')
            baseline        = nanmean(nanmean(freq.powspctrm(:,:,:,timeIdx),1),4); % first average over trials, then over time bins
            freq.powspctrm  = 10*log10(freq.powspctrm ./ repmat(baseline,size(freq.powspctrm,1),1,1,size(freq.powspctrm,4)));
        
        elseif strcmp(baselineSettings,"condition")
            
            for iCondi = 1:max(Xsel)

                trials = find(Xsel == iCondi); % trials in this bin
                   
                % Compute baseline averaged over trials:
                baseline  = nanmean(nanmean(freq.powspctrm(trials,:,:,timeIdx),1),4); % first average over time bins, then over trials     

                % Correct each trial by that baseline:
                freq.powspctrm(trials,:,:,:)    = 10*log10(freq.powspctrm(trials,:,:,:) ./ repmat(baseline,length(trials),1,1,size(freq.powspctrm,4)));
                   
            end % end of iCondi
            
        elseif strcmp(baselineSettings,"trial")
            fprintf("Subject %03d: Compute baseline per trial\n",iSub);
            baseline            = nanmean(freq.powspctrm(:,:,:,timeIdx),4); % keep trial, channel, freq bin, average over time bin
            freq.powspctrm      = 10*log10(freq.powspctrm ./ repmat(baseline,1,1,1,size(freq.powspctrm,4)));
            
        elseif strcmp(baselineSettings,"ft_trial")
            fprintf("Subject %03d: Compute baseline per trial using Fieldtrip\n",iSub);
            cfg                 = [];
            cfg.baseline        = par.TF.baselinetime;
            cfg.baselinetype    = 'db'; % use db; try out relative
            freq                = ft_freqbaseline(cfg, freq);        
            
        elseif strcmp(baselineSettings,"trend")
            fprintf("Subject %03d: Compute baseline per trial by detrending\n",iSub);
            % freq.powspctrm contains are trials, channels, frequencies, time bins
            
            for iChan = 1:length(freq.label) % iChan = 1;
                for iFreq = 1:length(freq.freq) % iFreq = 20;
                    
                    x = freq.trialinfo; % trial numbers
                    y = nanmean(freq.powspctrm(:,iChan,iFreq,timeIdx),4);  % average in baselineTimings
                    validIdx = ~isnan(y); % non-zero
                    p = polyfit(x(validIdx),y(validIdx),1); % fit regression
                    f = polyval(p,x); % predict values from regression for x
                    baseline(1:length(freq.trialinfo),iChan,iFreq) = f; % save for all trials
                    
                end % end for-loop for iFreq
            end % end for-loop for iChan 
            
            % Divide by baseline, convert to decibel (sometimes imaginary part of 0.000i):
            freq.powspctrm      = real(10*log10(freq.powspctrm ./ repmat(baseline,1,1,1,length(freq.time))));
            
        end % end baseline correction

        %% 6) Select trials into conditions:
        cfg = [];
        
        for iCondi = 1:max(Xsel)
           
            cfg.trials = find(Xsel == iCondi); % based on exact response
            nTrialCond(iSub, iCondi) = length(cfg.trials);
            
            fprintf('Subject %03d condition %d: Found %d trials\n', iSub, iCondi, length(cfg.trials))
            
            % 3) Average over selected trials:
            if ~isempty(cfg.trials) % if any trials in this condition
               TFall{subCount}.ValxAct{iCondi} = ft_freqdescriptives(cfg,freq); % define for this condition

            else % otherwise fill in NaNs
               fprintf("Subject %03d condition %d: Found no trials, store NaNs !\n",iSub, iCondi);
               TFall{subCount}.ValxAct{iCondi}.powspctrm = nan(size(freq.powspctrm,2),size(freq.powspctrm,3),size(freq.powspctrm,4));
               
            end % end of check for empty trials
        end % end of for loop for iCond
        
        fprintf('Store subject %03d\n', iSub)

        clear freq behav baselineData baseline
        
        % Update effective subject number
        subCount = subCount+1;
        
    end % end of for loop for subject   
    
    % Save data for all subjects.
    fprintf('Save all files\n')
    save(filenameTFall,'TFall','nTrialCond','-v7.3')

end % end if data does not yet exist

end % end of function.
