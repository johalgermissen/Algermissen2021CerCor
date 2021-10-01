function job = taft_preprocess_initialize_job(EEGdomain,iSub,ROIs2use,behav2use,selTrials)

% job = taft_preprocess_initialize_job(iSub,ROIs2use,behav2use,selTrials)
% 
% Defines a job for each subject give a set of ROIs (file names to use), 
% behavioral variables to use, selected trials. 
% Use this function to adjust various settings of the computation.
% Mind setting the root directory in taft_set_dirs.root().
% Mind to adjust paths to SPM (for pinv) and Fieldtrip.
%
% INPUTS:
% EEGdomain 	= string, either 'TF' (data in TF domain) or 'time' (data in time domain).
% iSub          = numeric scalar, subject identifier.
% ROIs2use      = cell, one or more string(s) specifying the file names
%               of fMRI ROIs data files to include in design matrix. For any volume-by-volume regressor (also nuisance regressor).
% behav2use       cell, one or more string(s) specifying the behavioral variables 
%               to include in design matrix (optional). For any trial-by-trial regressor.
% selTrials     = string, sub-selection of trials to use (optional).
%
% OUTPUTS:
% job           = structure with settings for creating TAfT object, specifically:
% .add_intercept= Boolean, whether to include intercept to design matrix or not. 
% .behavFile    = string, full file path for behavioral raw data file of respective subject.
% .behavReg 	= cell of strings, behavioral regressors to include as loaded and specified in taft_preprocess_load_behavior (optional), provided as input via behav2use.
% .dBconversion = Boolean, whether to convert outcome variable (EEG data) to decibel (10*log10(Y)) or not.
% .demeanX      = Boolean, whether to demean design matrix X or not.
% .demeanX      = Boolean, whether to demean outcome variable (EEG data) or not.
% .description  = string, contains important information about job to be saved in output file name, includes EEGdomain, lockSettings, regNames.
% .EEGdomain 	= string, either 'TF' (data in TF domain) or 'time' (data in time domain).
% .EEGfile      = string, full file path of subject-specific EEG data set to load.
% .foi          = vector of two numerics, frequency range of interest (depends on nFreq).
% .HRFtype      = string, whether to perform estimation of HRF amplitude for each trial separately ('trial') or in a GLM featuring all trials of a given block ('block').
% .lockSettings	= string, type of event-locked data to include, either 'stimlocked' or 'resplocked'.
% .nFreq        = numeric scalar, number of frequencies used in TF decomposition.
% .ons_unit 	= Boolean, whether onsets are given in seconds or not.
% .outputfile 	= string, complete path of output file (description with .mat at the end, with respective directory concatenated).
% .regNames 	= cell, one or more string(s) of all regressors (fMRI and behavior) to include in design matrix plus selected trials.
% .regType      = string, type of regression over trials performed, either 'pinv', 'fitlm', 'robustfit'
% .ROIs 		= further settings/file names to use for particular ROI for particular block of particular subject, specifically:
%   .ROIname 	= string, name of particular ROI, taken from ROI2use. 
%   .ROIdef.rawfMRIfile = string, file from which to load volume-by-volume MRI data.
%   .ROIdef.onsets  	= string, file from which to obtain trial onsets.
%   .ROIdef.nuisance 	= string, file from which to load volume-by-volume nuisance regressors from.
%   .ROIdef.betafMRIfile = string, file to which to save trial-by-trial HRF amplitude estimates.
%   .ROIdef.fitType  	= string, type of turning volume-by-volume data into trial-by-trial data, either 'HRF' (default) or 'avg' (for mp_Fellner regressor).
% .save         = Boolean, whether to save estimates of trial-by-trial HRF amplitude to disk.
% .selTrialsFile= string, name of file containing indices of selected trials to load. 
% .subID        = numeric scalar, subject ID.
% .TFtype       = string, typoe of TF decomposition used in TF decomposition, either 'hanning' or 'morlet'.
% .toi          = vector of two numerics, time range of interest (depends on lockSettings). 
% .TR           = numeric scalar, repetition time of MRI sequence (used for upsampling volumes) in seconds.
% .trialdur 	= numeric scalar, trial duration to use when epoching upsampled BOLD data in seconds (recommended by Hauser et al. 2015: 8 seconds).
% .ups          = numeric scalar, factor by which to upsample MRI sequence.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Adapted from Tobias Hauser (https://github.com/tuhauser/TAfT).
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/CueLockedAnalyses/TAfT/

%% Complete settings:

if nargin < 4
	behav2use = {};
end

if nargin < 5
	selTrials = 'all';
end

fprintf('Subject %03d: Initialize job\n',iSub);

job         = [];

job.subID   = iSub; % iSub
nBlocks     = 6;    % number of blocks per subject

%% Directories:

% Data directories:
dirs.root   	= taft_set_rootdir(); % /project/3017042.02
dirs.log        = fullfile(dirs.root,'Log');
dirs.behav      = fullfile(dirs.log,'Behavior/Data_beh_mat');    
dirs.fMRI       = fullfile(dirs.log,'fMRI');    
dirs.EEG        = fullfile(dirs.log,'EEG','CueLockedResults');    
dirs.TAfT       = fullfile(dirs.EEG,'TAfT_Betas');
if ~exist(dirs.TAfT,'dir'); mkdir(dirs.TAfT); end

% Software package directories:
dirs.fieldtrip  = '/home/common/matlab/fieldtrip'; % needs to be adapted to users' folder structure
dirs.SPM        = '/home/common/matlab/spm12'; % needs to be adapted to users' folder structure

%% Load necessary software packages:

if isempty(strfind(path,dirs.fieldtrip))
    addpath(dirs.fieldtrip);
    ft_defaults;
end
if isempty(strfind(path,dirs.SPM))
    addpath(dirs.SPM);
end

%% General behavior settings:

job.behavFile       = fullfile(dirs.root,'/Log/Behavior/Data_beh_mat',...
    sprintf('3017042.02_emmvdij_%03d_001_results.mat',iSub));    

%% General EEG settings:

% Time- or time-frequency domain:
if ~strcmp(EEGdomain,'TF') && ~strcmp(EEGdomain,'time'); error('Invalid EEGdomain specified'); end
job.EEGdomain       = EEGdomain; % time TF

job.nFreqs          = 15;   % number frequencies to be expected in EEG data
job.TFtype          = 'hanning'; % type of TF decomposition, either hanning or morlet

job.foi             = [1 job.nFreqs];
job.dBconversion    = false; % perform decibel conversion before regression or not

% Stimulus- or response-locked:
job.lockSettings    = 'stimlocked'; % stim or resp

% Select time of EEG for which to perform regression:
if strcmp(job.lockSettings,'stimlocked')
    job.toi = [0 2.0]; 
elseif strcmp(job.lockSettings,'resplocked')
    job.toi = [-1 0.5]; 
else
    error('Invalid lock settings')
end

% Directory from which to load EEG data:
if strcmp(job.EEGdomain, 'time')
    job.EEGfile    = fullfile(dirs.EEG,sprintf('finalPP/MRIpav_%03d_finalPP.mat',job.subID)); % where EEG data are
    job.EEGname    = 'scd';
elseif strcmp(job.EEGdomain, 'TF')
    job.EEGfile    = fullfile(dirs.EEG,sprintf('TF_%s_%s%d/MRIpav_%03d_TF_%s.mat',job.lockSettings,job.TFtype,job.nFreqs,job.subID,job.lockSettings)); % where EEG data are
    job.EEGname    = 'freq';
else
    error('unknown data format')
end

%% General fMRI settings:

% Fit betas based on trial-wise or block-wise GLM:
job.HRFtype     = 'trial'; % trial or block

% Settings for BOLD deconvolution:
job.TR          = 1.4;  % Original repetition time
job.hpass       = 128;  % High-pass filter about 128 sec.
job.ups         = 10;   % Factor by which to upsample original RT, new RT is job.TR/job.ups (Hauser, Hunt et al., 2015 achieved new RT of 185 ms)
job.ons_unit    = 1;    % Boolean that onsets are given in seconds
job.trialdur    = 8;    % Trial durations when epoching upsampled BOLD data into trials (8 sec. as Hauser, Hunt et al. 2015)

% For motion regressors: Average over shorter window (1-2 volumes):
if strcmp(ROIs2use,'mp_Fellner')
    job.trialdur = 2;
end

% Settings for regression:
job.demeanX         = true;     % demean design matrix (fMRI and behavior) before regression
job.demeanY         = true;     % demean EEG data for each channel/frequency/time bin before regression
job.add_intercept   = false;    % add intercept for GLMs
job.save            = false;    % save trial-by-trial HRF amplitudes as separate file  
job.regType         = 'pinv';   % type of regression: pinv (fast), ols, robust (quite slow)

if ischar(ROIs2use) % if input only string, not cell
    ROIs2use = {ROIs2use};
end

%% Regressors:

% 1) BOLD ROI regressors:
job.regNames    = ROIs2use;

% 2) Behavioral regressors:
if ~isempty(behav2use) % if specified: add to regNames 
    job.behavReg = behav2use;
    job.regNames = [job.regNames job.behavReg];
end

% 3) Selected trials to use:
if ~strcmp(selTrials,'all') % if specified: add to regNames
    job.selTrials  	= selTrials;
    job.selTrialsFile   = fullfile(dirs.TAfT,sprintf('TAfT_selTrials/selTrials_%s_Sub%03d.mat',job.selTrials,iSub));
    job.regNames        = [job.regNames job.selTrials];
end

%% Name of output file:

% Basic information: EEG domain, lock setting, all regressors, HRF
% deconvolution type, trial duration, rergression type
job.description         = sprintf('TAfT_%s_%s_%s_%s_%dsec_%s',...
    job.EEGdomain,job.lockSettings,strjoin(job.regNames,'_'),job.HRFtype,job.trialdur,job.regType);

% Add type of TF decomposition and # frequencies:
if strcmp(job.EEGdomain,'TF')
    job.description         = sprintf('%s_freq%02d',job.description,job.nFreqs); % add number of frequencies
end

% Add if decibel conversion:
if job.dBconversion
    job.description     = sprintf('%s_dB',job.description);
end

% Specify full output directory:
job.outputFile          = fullfile(dirs.TAfT,sprintf('/%s.mat',job.description));

% Warn of overwriting file:
if exist(job.outputFile,'file')
    warning('Output file already exists');
end

%% Initialize directories for each ROI, for each block:

for iROI = 1:length(ROIs2use) % loop through ROIs
    
   for iBlock = 1:nBlocks % loop through ROIs
    
        job.ROIs(iROI).ROIname = ROIs2use{iROI}; %  name of ROIs just carried forward
        
        % Where BOLD time series, onsets and nuisance regressors are expected to be:
        fMRIsubdir  = fullfile(dirs.fMRI,sprintf('sub-%03d/FEAT_Block%d.feat/AROMA',iSub,iBlock));

        % BOLD time series extracted from ROI:
        job.ROIs(iROI).ROIdef(iBlock).rawfMRIfile = {fullfile(fMRIsubdir,sprintf('%s.txt',job.ROIs(iROI).ROIname))}; % where to load volumes from

        % Trial onsets (text file):
        job.ROIs(iROI).ROIdef(iBlock).onsets = {fullfile(fMRIsubdir,'trialOnsets.txt')}; % trial onsets

        % Potential nuisance regressors:
        if ~strcmp(job.ROIs(iROI).ROIname,'mp_Fellner')
            job.ROIs(iROI).ROIdef(iBlock).nuisance = {fullfile(fMRIsubdir,'NuisanceRegressors.txt')}; % where to load volumes from
        end

        % Output directory:
        job.ROIs(iROI).ROIdef(iBlock).betafMRIfile = {fullfile(dirs.TAfT,sprintf('betas_ROI_%s_sub%03d_block%d.txt', job.ROIs(iROI).ROIname,iSub,iBlock))};
        
        % Deconvolution type:
        job.ROIs(iROI).ROIdef(iBlock).fitType = 'HRF';        
        % Or just average:
        if strcmp(job.ROIs(iROI).ROIname,'mp_Fellner')
            job.ROIs(iROI).ROIdef(iBlock).fitType = 'avg';
        end
        
   end % end iBlock
end % end iROI

% Print all settings to screen:
fprintf('Create job in %s domain, %s-locked, %s-wise GLM for fitting HRF, trial duration %d sec., %s regression, \nRegressors: \n%s\n',job.EEGdomain,job.lockSettings,job.HRFtype,job.trialdur,job.regType,strjoin(job.regNames,'\n'));

% END
