function EEGfMRIPav_CueLocked_7_TF(rootdir,subVec,lockSettings,stimERPcor,respERPcor)

% EEGfMRIPav_CueLocked_7_TF(rootdir,subVec,lockSettings,stimERPcor,respERPcor)
% 
% Reject trials, compute surface LaPlacian filter (and at the same time
% interpolate bad channels). Final steps of pre-processing.
% 
% INPUTS:
% rootdir           = string, root directory of project
% subVec            = vector of integers, numbers of subjects to process
% lockSettings      = string, type of time-locking before TF decomposition,
% either 'stimlocked' (default) or 'resplocked'
% stimERPcor        = Boolean, whether to trial-wise subtract condition-wise 
% ERP before TF decompositions or not (default: false)
% respERPcor        = Boolean, whether to trial-wise subtract condition-wise 
% ERP before TF decompositions or not (default: false)
%
% OUTPUTS:
% saves data to disk under dirs.TF.
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
    rootdir = preprocessing_set_rootdir(); % '/project/3017042.02';
    fprintf('rootdir unspecified, assume %s\n',rootdir)
end

if ~exist('subVec','var')
    subVec = 1:36; 
    fprintf('subVec unspecified, assume %s\n',strjoin(string(subVec),', '));
end

if ~exist('lockSettings','var')
    lockSettings = 'stimlocked'; % stimlocked resplocked
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
dirs.TF     = fullfile(dirs.results,sprintf('TF_%s_%s%d',lockSettings,par.TF.TFtype,par.TF.nFreqs));

if stimERPcor 
    fprintf('Subtract stimulus-locked ERP before TF decomposition \n');
    dirs.TF = sprintf('%s_stimERPcor',dirs.TF); % overwrite target file name
end
if respERPcor
    fprintf('Subtract response-locked ERP before TF decomposition \n');
    dirs.TF = sprintf('%s_respERPcor',dirs.TF); % overwrite target file name
end

if ~exist(dirs.TF,'dir'); mkdir(dirs.TF); end

% Time-points to decompose based on lockSettings:
if strcmp(lockSettings,'stimlocked')
    par.TF.toi4tf          = -1:.025:2; % from 1 sec before till 2 sec after cue presentation; windows of 0.025 sec
elseif strcmp(lockSettings,'resplocked')
    par.TF.toi4tf          = -1.5:.025:1;
else
    error('Invalid lock settings')
end

fprintf('Lock setting is %s \n',lockSettings);
fprintf('Perform TF decomposition from %.03f to %.03f sec. \n',...
    par.TF.toi4tf(1),par.TF.toi4tf(end));

%% Detect input files:

% Detect all available input folders (all subjects):
tmp         = dir(fullfile(dirs.finalPP,'*.mat'));
fileList    = {tmp.name};
nInput      = length(fileList);
fprintf('Found data files from %d subjects\n',nInput);

% Extract subject numbers from folders available:
subNames    = str2double(extractBetween(fileList,8,10));

% Extract indices of selected subjects:
selIndices  = find(ismember(subNames,subVec));

%% Loop over subjects:

for iSub = selIndices % %iSub = 1;
    
    %% Create output file name:
    
    outputFile      = fullfile(dirs.TF, sprintf('MRIpav_%03d_TF_%s.mat',iSub,lockSettings));
    fprintf('EEG output file is %s\n',outputFile);
    
    % Check if output file already exists:
    if exist(outputFile,'file')
        warning('Subject %03d: File %s already exists, skip subject\n', iSub, outputFile);
        continue;
    end

    %% Create input file name:
    
    inputFile       = fullfile(dirs.finalPP, fileList{iSub}); % go to subject directory
    fprintf('EEG input file is %s\n',inputFile);
       
    %% Load data from previous step:
    
    fprintf('Subject %03d: Start loading time-domain data\n',iSub)
    tmp     = load(inputFile);
    scd     = tmp.scd;
    hist    = tmp.hist;
    fprintf('Subject %03d: Finished loading time-domain data\n',iSub)

    % Filter out non-EEG channels:
    scd.elec.elecpos    = scd.elec.elecpos(par.chan.EEG,:); 
   
    %% Load rejected trials:
    
    fprintf('Subject %03d: Load rejected trials\n',iSub)
    
    filenameRejTrials = fullfile(dirs.rejTrials,sprintf('MRIpav_%03d_rejTrials.mat',iSub));
    load(filenameRejTrials);

    %% Load and store relevant behavioral data.
    
    fprintf('Subject %03d: Load behavior\n',iSub)

    behav = load_recode_behav(rootdir,iSub,rejectedtrials); 
        
    %% Subtract stimulus-locked ERP:
    
    if stimERPcor
        
        % --------------------------------------------------------------- %
        % 1) Convert trials stored as cells into matrix:
        cfg = [];
        cfg.keeptrials = 'yes';
        scd = ft_timelockanalysis(cfg,scd);
        
        % --------------------------------------------------------------- %
        % 2) Baseline correction:
        % a) With Fieldtrip:
%         cfg = [];
%         cfg.baseline = par.TF.baselineTimings; % -250 - 50 ms
%         base_scd = ft_timelockbaseline(cfg,scd);

        % b) Regression based:
        startIdx    = dsearchn(scd.time',par.TF.baselinetime(1)); % sample index for start of baseline period
        stopIdx     = dsearchn(scd.time',par.TF.baselinetime(end)); % sample index for stop of baseline period
        scd.trial   = regress_baseline(scd.trial,scd.trialinfo,startIdx,stopIdx,'subtraction');
        
        % --------------------------------------------------------------- %
        % 3) Compute and subtract ERP per condition:
        for iResp = [1 0] % Go/NoGo action: iResp = 1; iResp = 0;
            
            for iVal = [1 2] % Win/Avoid cue: iVal = 1; iVal = 2;
                
               fprintf('Subject %03d: Correct for stimulus-locked ERP in iResp = %d, iVal = %d\n', iSub, iResp, iVal);
                
                % a) Select trials per condition:
                selTrials = find(behav.val == iVal & behav.go == iResp); % 
               
                % b) Compute ERP for selected trials, subtract from selected trials:
                scd_cor = subtract_ERP(scd,scd,selTrials,selTrials,'Fieldtrip'); % entire data with Fieldtrip
                
                % c) Insert selectively corrected trials into original object:
                scd.trial(selTrials,:,:) = scd_cor.trial(selTrials,:,:);
                
                % d) Do baseline-correction again:
                startIdx    = dsearchn(scd.time',par.TF.baselinetime(1));
                stopIdx     = dsearchn(scd.time',par.TF.baselinetime(end));
                scd.trial   = regress_baseline(scd.trial,scd.trialinfo,startIdx,stopIdx,'subtraction');
                
            end % end iResp
        end % end iVal
    end % end stimERPcor
    
    %% Use response-locked data:
    
    if strcmp(lockSettings,'resplocked')
        
        scd = redefine_resplock(par,behav,scd);

    end
    
    %% Subtract response-locked ERP:
    
    if respERPcor
        
        % Use response-locked data for computing ERPs:
        respscd = redefine_resplock(iSub,par.sRate,par.chan.EEG,behav,scd);
        
        % Selected trials: Subtract ERP only from Go trials:
        goIdx   = find(behav.go); % only on Go trials
        
        % Timing of subtraction:
        RT2use  = behav.RT; % retrieve RTs
        RT2use(RT2use > 1.3035) = NaN; % remove too long RTs
        RT2use  = RT2use(~isnan(RT2use)); % remove NaNs
        
        % Compute ERP for Go trials (response-locked), subtract from Go trials (stimulus-locked):
        scd = subtract_ERP(scd,respscd,goIdx,goIdx,RT2use,'Manual',RT2use,0.7);
        
    end
    
    %% Perform Morlet wavelet convolution.
    
    if strcmp(par.TF.TFtype,'morlet')

        fprintf('Subject %03d: Start TF decomposition using Morlet wavelets\n',iSub)
        
        cfg                 = [];
        cfg.method          = 'wavelet';
        cfg.output          = 'pow'; %'fourier';% to also keep the phase.
        cfg.keeptrials      = 'yes';
%         cfg.foilim          = [parTF.minFreq parTF.maxFreq]; % either frequency min and max
        cfg.foi             = par.TF.freq4tf;      % or directly frequencies of interest.
        cfg.width           = par.TF.width4tf;     % number of cycles.
        cfg.toi             = par.TF.toi4tf;       % times (s) to centre the analysis window.
        cfg.channel         = 'all';
        freq                = ft_freqanalysis(cfg,scd);
        
    %% Perform multi-taper-method convolution with Hanning tapers:
    
    elseif strmcp(par.TF.TFtype,'hanning')
        
        fprintf('Subject %03d: Start TF decomposition using Hanning tapers\n',iSub)
        
        cfg                 = [];
        cfg.method          = 'mtmconvol'; % mtmfft mtmconvol
        cfg.taper           = 'hanning';
%         cfg.foilim          = [parTF.minFreq parTF.maxFreq]; % either frequency min and max        
        cfg.foi             = par.TF.freq4tf; % or directly frequencies of interest.
        cfg.t_ftimwin       = repmat(par.TF.ftimwin,1,length(cfg.foi)); % length of time window for each frequency
        cfg.toi             = par.TF.toi4tf; % time bins for which to compute frequency
        cfg.pad             = par.TF.pad; % pad up to 8 sec. to get well-behaved frequency bins
        cfg.output          = 'pow'; %'fourier' to also keep the phase.
        cfg.keeptrials      = 'yes';
        cfg.channel         = 'all';
        freq                = ft_freqanalysis(cfg,scd);

    end
    
    %% Save data.
    
    fprintf('Subject %03d: Start saving data\n',iSub)
    
    hist.TF.dirs       = dirs;
    hist.TF.par        = par;

    save(outputFile,'freq','hist','-v7.3')
    clear scd freq par behav 
    
    fprintf('Subject %03d: Finished saving data :-)\n',iSub)  
        
end % iSub-loop.

fprintf('Done :-)\n');

end % end of function.
