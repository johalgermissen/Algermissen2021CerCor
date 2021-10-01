function EEGfMRIPav_CueLocked_8_extract_singleTrial(rootdir,subVec,lockSettings,contrast,channels,band,GLMID)

% EEGfMRIPav_CueLocked_8_extact_singleTrial(rootdir,subVec)
% 
% Compute trial-by-trial mean power within mask defined previously in
% two-condition permutation test (from EEGfMRIPav_grouplevel_TF_test.m).
% 
% INPUTS:
% rootdir           = string, root directory of project.
% subVec            = vector of integers, numbers of subjects to process.
% lockSettings      = string, type of time-locking before TF decomposition,
% either 'stimlocked' (default) or 'resplocked'.
% contrast          = string, contrast from which mask was created.
% (featured in mask name), either 'Congruency' (default) or 'Go'.
% channels          = string, list of channels (without separators) for
% which mask was created (featured in mask name), default 'CzFCzCz'.
% band              = string, band for which mask was created (featured in
% mask name), either 'alpha' (default) or 'theta'.
% GLMID             = string, name of GLM in which timingfiles to save
% regressors (default: 3).

%% Complete input settings:

% clear all; close all; clc
% dbstop if error

if ~exist('rootdir','var')
    rootdir = grouplevel_set_rootdir(); % '/project/3017042.02';
    fprintf('rootdir unspecified, assume %s\n',rootdir)
end

if ~exist('iSub','var')
    subVec = 1:36; 
    fprintf('subVec unspecified, assume %s\n',strjoin(string(subVec),', '));
end

if ~exist('lockSettings','var')
    lockSettings = 'stimlocked'; % stimlocked resplocked
end

if ~exist('contrast','var')
    contrast = 'Congruency'; % Go Congruency
end

if ~exist('channels','var')
    channels = 'CzFCzFz'; % CzFCzFz
end

if ~exist('band','var')
    band = 'alpha'; % alpha, theta
end

if ~exist('GLMID','var')
    GLMID = '3';
end

nSub            = length(subVec);

%% Set directories:

% Add helper files to path:
fprintf('Add helper files to path\n');
addpath(fullfile(rootdir,'Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Helpers'));

% Set directories and parameters:
dirs        = set_dirs(rootdir);
par         = set_par();

% Add TF directory based on TF settings:
dirs.TF     = fullfile(dirs.results,sprintf('TF_%s_%s%d',...
    lockSettings,par.TF.TFtype,par.TF.nFreqs));

% fMRI (output):
dirs.fMRI       = fullfile(rootdir,'Log','fMRI'); % subject-specific directory created later

%% Create names of mask:

maskName  = sprintf('Mask_%s_%s_%s_%s.mat', contrast, channels, band, lockSettings);
fprintf('Selected mask is %s\n',maskName);

%% Load mask:

fprintf('Load mask %s \n', maskName);    
load(fullfile(dirs.mask,maskName));

%% Locate single-subject EEG files:

tmp         = dir(fullfile(dirs.TF,sprintf('*_%s.mat',lockSettings)));
fileList    = {tmp.name};
nInput      = length(fileList);
fprintf('Found data files from %d subjects\n',nInput);

% Extract subject numbers from folders available:
subNames    = str2double(extractBetween(fileList,8,10));

% Extract indices of selected subjects:
selIndices  = find(ismember(subNames,subVec));
nSub        = length(selIndices);

%% Loop over subjects: 

allImputeMean       = nan(nSub,par.nCond); % initialize mean per condition per subject
allImputeCount      = nan(nSub,par.nCond); % initialize count of NaN per condition per subject

for iSub = selIndices

    %% Create output file name:
    
    dirs.target  = fullfile(dirs.fMRI,sprintf('sub-%03d/GLM%s/timings_regressors',iSub,GLMID));
    if ~exist(dirs.target,'dir'); mkdir(dirs.target); end
    
    outputFile  = fullfile(dirs.target,sprintf('EEG_%s.txt',targetFile));
    fprintf('EEG output file is %s\n',outputFile);
    
    % Check if output file already exists:
    if exist(outputFile,'file')
        warning('Subject %03d: File %s already exists, skip subject\n', iSub, outputFile);
        continue;
    end

    %% Create input file name:
    
    inputFile  = fullfile(dirs.TF, fileList{iSub}); % go to subject directory
    fprintf('EEG input file is %s\n',inputFile);

    %% Load EEG:
    
    fprintf('Subject %03d: Start loading TF data\n',iSub)
    tmp     = load(inputFile);
    freq    = tmp.freq;
    hist    = tmp.hist;
    fprintf('Subject %03d: Finished loading TF data\n',iSub)

    %% Load rejected trials: 

    fprintf('Subject %03d: Load rejected trials\n',iSub)
    
    filenameRejTrials = fullfile(dirs.rejTrials,sprintf('MRIpav_%03d_rejTrials.mat',iSub));
    load(filenameRejTrials);

    %% Load and store relevant behavioral data.
    
    fprintf('Subject %03d: Load behavior\n',iSub)

    behav = load_recode_behav(rootdir,iSub,rejectedtrials); 
      
    %% Loop over trials and extract summary measure of mask:
       
    nNonRejTrial = length(freq.trialinfo); % number non-rejected trials
    powSumRaw  = nan(nNonRejTrial,1); % initialize 
    
    for iTrial = 1:nNonRejTrial % iTrial = 1;
        trlPow      = squeeze(freq.powspctrm(iTrial,:,:,:)); % power for this trial
        if sum(size(trlPow) ~= size(freqMask))>0; error('Data and mask have different dimensions'); end
        trlPowMask  = trlPow .* freqMask; % multiply power of this trial with mask
        powSumRaw(iTrial) = sum(trlPowMask(:)); % sum across all dimensions
    end
    
    %% Compute mean per action x valence condition:
    
    if(nNonRejTrial ~= length(behav.go)); error('EEG and behavior have different numbers of trials'); end

    fprintf('Subject %03d: Compute condition means \n',iSub);
    
    imputeVec   = nan(par.nCond,1); % initialize vector of condition means
    iCond = 0;
    for iResp = [1 0] % iResp = 1; loop over Go/NoGo
        for iVal = [1 2] % iVal = 1; % loop over Win/Avoid
            iCond = iCond + 1;
            selTrials = behav.go == iResp & behav.val == iVal; % all trials in this condition
            imputeVec(iCond)    = mean(powSumRaw(selTrials)); % condition mean
        end
    end
    
    %% Interpolate rejected & invalid trials:
    
    fprintf('Subject %03d: Interpolate trials \n',iSub);
    
    powSumImp   = nan(par.nTrial,1); % initialize output power per trial
    valTrlCount = 0; % initialize valid trial count
    imputeCount = zeros(par.nCond,1); % count # imputations per condition

    for iTrial  = 1:par.nTrial
        
        if ismember(iTrial,behav.recoded.trlidx) % if valid trial: insert
            valTrlCount             = valTrlCount + 1; % increase valid trial count
            powSumImp(iTrial)       = powSumRaw(valTrlCount); % retrieve valid power value
            
        else % retrieve condition of trial, impute:
            iCond                   = 2*(1-behav.go(iTrial)) + behav.val(iTrial); % determine condition of this trial
            powSumImp(iTrial)       = imputeVec(iCond); % retrieve value used for imputation
            imputeCount(iCond)      = imputeCount(iCond) + 1; % increase imputation count
        end
    end
    
    allImputeMean(iSub,:)   = imputeVec; % save mean imputed values per subject
    allImputeCount(iSub,:)  = imputeCount; % save counts of imputed values per subject
    fprintf('Subject %03d: Imputed %s trials \n',iSub,strjoin(string(imputeCount)));

    %% Check if still NaN left:

    if sum(isnan(powSumImp))>0; error('Found NaN'); end
    
    %% Demean:
    
    fprintf('Subject %03d: Demean power values \n',iSub);
    powSumImp      = powSumImp - mean(powSumImp);
    
    %% Save:
    
    save(outputFile,'powSumImp','-ascii')

end
fprintf('Finished :-)\n');

end % end of function.