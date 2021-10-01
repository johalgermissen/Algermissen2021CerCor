function EEGfMRIPav_CueLocked_8_TF_grouplevel_create(rootdir,subVec,lockSettings,baselineSettings,baselineTimings,...
    responseSettings,actionSettings,accSettings,...
    TFtype,nFreq,stimERPcor,respERPcor)

% EEGfMRIPav_CueLocked_8_TF_grouplevel_create(rootdir,subVec,lockSettings,baselineSettings,baselineTimings,...
%     responseSettings,actionSettings,accSettings,...
%     TFtype,nFreq,stimERPcor,respERPcor)
% 
% Load EEG, rejected trials, behavior; apply baseline correction; sort
% trials into specified conditions, aggregate within conditions, save.
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
% which to use.
% responseSettings  = string, type of response setting for which
% conditions are split, 'Go' (Go/ NoGo) or 'Hand' (Left Go/ Right Go/
% NoGo.)
% actionSettings    = string, type of action setting for which
% conditions are split, 'exeAct' (executed action) or 'reqAct' (required
% action).
% accSettings       = string, type of accuracy setting for which 
% conditions are split, 'correct', 'incorrect', 'bothAcc' (both
% accuracies separately), 'noAcc' (no distinction), 'reqAct' (required
% action).
% TFtype            = type of TF decomposition, 'hanning' or 'morlet'.
% nFreq             = numeric, number of frequency bins decomposed
% (default: 15).
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
    baselineTimings = [-0.250 -0.050]; % 
    fprintf('baselineTimings unspecified, assume %s\n',strjoin(string(baselineTimings),' - '));
end

if ~exist('responseSettings','var')
    responseSettings = 'Go'; % Go Hand
    fprintf('responseSettings unspecified, assume %s\n',responseSettings);
end

if ~exist('actionSettings','var')
    actionSettings = 'exeAct'; % exeAct reqAct
    fprintf('actionSettings unspecified, assume %s\n',actionSettings);
end

if ~exist('accSettings','var')
    accSettings = 'correct'; % correct incorrect bothAcc noAcc reqAct
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

%% Automatically fill in downstream settings:

if strcmp(responseSettings,'Go')
    allResp = [1 0]; % Go/NoGo.
    nCond = 4;

elseif strcmp(responseSettings,'Hand')
    allResp = [101 97 0]; % left/right/nogo.
    nCond = 6;

else
    fprintf('Unknown response setting')
end

if strcmp(accSettings,'bothAcc') || strcmp(accSettings,'reqAct')
    allAcc = [1 0];

elseif strcmp(accSettings,'correct')
    allAcc = 1;

elseif strcmp(accSettings,'incorrect')
    allAcc = 0;

else
    fprintf('Unknown accuracy setting')
end

nCond = nCond * length(allAcc); % times accuracy conditions

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

%% Create final file name:

outputFile = sprintf('EEGfMRIPav_TF_%s_baseCor_%s_baseTime_%03dms_%03dms_resp_%s_%s_%s_%s%d', ...
    lockSettings,baselineSettings,baselineTimings(1)*-1000,baselineTimings(end)*-1000,responseSettings,actionSettings,accSettings,TFtype,nFreq);

if stimERPcor
    outputFile = [outputFile '_stimERPcor'];
end

if respERPcor
    outputFile = [outputFile '_respERPcor'];
end

% Add directory and .mat ending:
outputFile = fullfile(dirs.TFgroup,sprintf('%s.mat',outputFile));

fprintf('Selected file is %s\n',outputFile)

%% Loop over subjects:

if ~exist(outputFile,'file')

    subCount = 1; % relative subject number, i.e. also works if certain numbers are missing

    TFall = cell(nSub,1);
    nTrialCond = nan(nSub,nCond);
    
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
               
        %% 4) Baseline-correction per subject (and transformation into decibel): 
        
        % Load stimulus-locked baselines for response-locked data:
        if strcmp(lockSettings,'stimlocked')
            baselineData = freq; % use same data

        elseif strcmp(lockSettings,'resplocked')
            fprintf('Subject %03d: Additionally load stimlocked data for baseline correction\n', iSub)
            tmp = load(fullfile(dirs.results,[sprintf('/TF_subjectlevel_stimlocked_%s%d/',TFtype,nFreq) fileList{iSub}(1:end-14) 'stimlocked.mat']),'freq'); % use stimulus-locked data
            baselineData = tmp.freq; 
            clear tmp
        
        else
            error('Unknown lock settings')
        end
        
        % Retrieve indices of baseline timings:
        timeIdx = dsearchn(round(freq.time,3)',baselineTimings'); % select indices for baseline time window
        
        if length(timeIdx) > 1; timeIdx = timeIdx(1):timeIdx(end); end % extend to interval if start and end given
        fprintf('Subject %03d: Compute baseline, setting %s\n', iSub, baselineSettings);
        
        if strcmp(baselineSettings,'all')
            baseline        = nanmean(nanmean(baselineData.powspctrm(:,:,:,timeIdx),1),4); % first average over trials, then over time bins; keep 4D format
            freq.powspctrm  = 10*log10(freq.powspctrm ./ repmat(baseline,size(freq.powspctrm,1),1,1,size(freq.powspctrm,4)));
        
        elseif strcmp(baselineSettings,'condition')
                        
            iCondi = 0; % iCondi is stimulus number (Response x Valence, i.e. 6)
           
            for iAcc = allAcc % Correct/ incorrect: iAcc = 1; iAcc = 0;
                for iResp = allResp % Go/NoGo action: iResp = 1; iResp = 0;
                    for iVal = 1:2 % Win/Avoid cue: iVal = 1; iVal = 2;

                        iCondi = iCondi + 1; % update condition count.

                        fprintf('Subject %03d condition %d: Compute baseline subject per val and %s for iResp %d, iVal %d\n',iSub,iCondi,responseSettings,iResp,iVal)                   
                        if strcmp(responseSettings,'Go') && strcmp(actionSettings,'exeAct') && strcmp(accSettings,'reqAct')
                            selTrials = find(behav.val == iVal & behav.isgo == iResp & behav.reqgo == iAcc); % based on exact response

                        elseif strcmp(responseSettings,'Go') && strcmp(actionSettings,'reqAct')
                            selTrials = find(behav.val == iVal & behav.reqgo == iResp & behav.accuracy == iAcc); % based on exact response

                        elseif strcmp(responseSettings,'Go') && strcmp(actionSettings,'exeAct')
                            selTrials = find(behav.val == iVal & behav.isgo == iResp & behav.accuracy == iAcc); % based on exact response

                        elseif strcmp(responseSettings,'Hand') && strcmp(actionSettings,'reqAct')
                            selTrials = find(behav.val == iVal & behav.reqresp == iResp & behav.accuracy == iAcc); % based on exact response

                        elseif strcmp(responseSettings,'Hand') && strcmp(actionSettings,'exeAct')
                            selTrials = find(behav.val == iVal & behav.resp == iResp & behav.accuracy == iAcc); % based on exact response

                        else
                            fprintf('Invalid response specification: %s',responseSettings)
                        end % end of trial selection
                    
                        % Compute baseline averaged over trials:
                        baseline     = nanmean(nanmean(baselineData.powspctrm(selTrials,:,:,timeIdx),1),4); % first average over time bins, then over trials     

                        % Correct each trial by that baseline:
                        freq.powspctrm(selTrials,:,:,:)    = 10*log10(freq.powspctrm(selTrials,:,:,:) ./ repmat(baseline,length(selTrials),1,1,size(freq.powspctrm,4)));
                    
                    end % end of iVal
                end % end of iResp
            end % end of iACC            
            
        elseif strcmp(baselineSettings,'trial')
            fprintf('Subject %03d: Compute baseline per trial\n',iSub)
            baseline   = nanmean(baselineData.powspctrm(:,:,:,timeIdx),4); % keep trial, channel, freq bin, average over time bin
            freq.powspctrm      = 10*log10(freq.powspctrm ./ repmat(baseline,1,1,1,size(freq.powspctrm,4)));

        elseif strcmp(baselineSettings,'ft_trial')
            fprintf('Subject %03d: Compute baseline per trial using Fieldtrip\n',iSub)
            % no adaptation possible for response-locked data
            cfg                 = [];
            cfg.baseline        = baselineTimings;
            cfg.baselinetype    = 'db'; % use db;
            freq                = ft_freqbaseline(cfg, freq);        
        
        elseif strcmp(baselineSettings,'trend')
            fprintf('Subject %03d: Compute baseline per trial by detrending\n',iSub)
            
            for iChan = 1:length(freq.label) % iChan = 1;
                for iFreq = 1:length(freq.freq) % iFreq = 20;
                    
                    x = freq.trialinfo; % trial numbers
                    y = nanmean(baselineData.powspctrm(:,iChan,iFreq,timeIdx),4);  % average in baselineTimings
                    validIdx = ~isnan(y); % non-zero
                    p = polyfit(x(validIdx),y(validIdx),1); % fit regression
                    f = polyval(p,x); % predict values from regression for x
                    TFall{subCount}.power.baseline(1:length(freq.trialinfo),iChan,iFreq) = f; % save for all trials
                    
                end % end for-loop for iFreq
            end % end for-loop for iChan 
            
           % Divide by baseline, convert to decibel (sometimes imaginary part of 0.000i):
           freq.powspctrm      = real(10*log10(freq.powspctrm ./ repmat(TFall{subCount}.power.baseline,1,1,1,length(freq.time))));
            
        end % end baseline correction
        
        %% 5) Select trials into conditions:
        % Go: Store data in TFall.ValxAct (G2W, G2A, NG2W, NG2A).
        % Hand: Store data in TFall.ValxAct (G2WL, G2AL, G2WR, G2AR, NG2W, NG2A).
        
        cfg = []; iCondi = 0; % iCondi is stimulus number (Response x Valence, i.e. 6)
        
        for iAcc = allAcc % Correct/ incorrect: iAcc = 1; iAcc = 0;
            for iResp = allResp % Go/NoGo action: iResp = 1; iResp = 0;
                for iVal = 1:2 % Win/Avoid cue: iVal = 1; iVal = 2;
                   
                   iCondi = iCondi + 1; % update condition count
        
                   % 1) Select trial numbers of condition:
                   fprintf('Subject %03d condition %d: Start TF extraction per val and %s for iAcc = %d, iResp = %d, iVal = %d\n', iSub, iCondi, responseSettings, iAcc, iResp, iVal);
                   if strcmp(responseSettings,'Go') && strcmp(actionSettings,'exeAct') && strcmp(accSettings,'reqAct')
                       cfg.trials = find(behav.val == iVal & behav.isgo == iResp & behav.reqgo == iAcc); % based on exact response
                   
                   elseif strcmp(responseSettings,'Go') && strcmp(actionSettings,'reqAct')
                       cfg.trials = find(behav.val == iVal & behav.reqgo == iResp & behav.accuracy == iAcc); % based on exact response
                   
                   elseif strcmp(responseSettings,'Go') && strcmp(actionSettings,'exeAct')
                       cfg.trials = find(behav.val == iVal & behav.isgo == iResp & behav.accuracy == iAcc); % based on exact response
                   
                   elseif strcmp(responseSettings,'Hand') && strcmp(actionSettings,'reqAct')
                       cfg.trials = find(behav.val == iVal & behav.reqresp == iResp & behav.accuracy == iAcc); % based on exact response
                   
                   elseif strcmp(responseSettings,'Hand') && strcmp(actionSettings,'exeAct')
                       cfg.trials = find(behav.val == iVal & behav.resp == iResp & behav.accuracy == iAcc); % based on exact response
                   
                   else
                       fprintf('Invalid response specification: %s',responseSettings)
                   end % end of trial selection

                  nTrialCond(iSub, iCondi) = length(cfg.trials); % cell size for this condition

                   % 2) Average over selected trials:
                   if ~isempty(cfg.trials) % if any trials in this condition
                       TFall{subCount}.ValxAct{iCondi} = ft_freqdescriptives(cfg,freq); % define for this condition
                   
                   else % otherwise fill in NaNs
                       fprintf('Subject %03d condition %d: Found no trials in condition iResp = %d, iVal=%d, store NaNs! \n',iSub, iCondi, iResp, iVal);
                       TFall{subCount}.ValxAct{iCondi}.powspctrm = nan(size(freq.powspctrm,2),size(freq.powspctrm,3),size(freq.powspctrm,4));
                       
                   end % end of check for empty trials
                end % end of for loop for iVal
            end % end of for loop for iResp
        end % end of for loop for iAcc
        
        %% Store data, update subject count:
        
        fprintf('Store subject %03d\n', iSub)
        
        clear freq behav baselineData baseline
        
        % Update effective subject number
        subCount = subCount+1;
        
    end % end of for loop for subject   
    
    %% Save data for all subjects:
    
    fprintf('Save file under %s\n',outputFile)
    save(outputFile,'TFall','nTrialCond','-v7.3')
    fprintf('Finished :-)\n')

end % end if ~exist(outputFile)

fprintf('Done :-)\n');

end % end of function.
