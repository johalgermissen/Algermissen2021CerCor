function [job, dirs, betas] = taft_postprocess_load_job(EEGdomain,ROIs2use,behav2use,selTrials)

% taft_preprocess_main()
% 
% Interactive script.
% Load beta-weights from multiple linear regression of EEG on fMRI and
% behavior per subject.
% Only regressors and selected trials as input; all other settings have to
% be changed within this file as fields of job object.
% Mind setting the root directory in taft_set_rootdir().
% Mind adjusting the paths to SPM and Fieldtrip.
% 
% INPUTS:
% EEGdomain 	= string, either 'TF' (data in TF domain) or 'time' (data in time domain).
% ROIs2use      = cell, one or more string(s) specifying the file names
%               of fMRI ROIs data files to include in design matrix. For any volume-by-volume regressor (also nuisance regressor).
% behav2use     = cell, one or more string(s) specifying the behavioral variables 
%               to include in design matrix (optional). For any
%               trial-by-trial regressor. Default empty ({}).
% selTrials     = string, sub-selection of trials to use (optional).
%               Default 'all' (all trials, no selection). 
%
% OUTPUTS:
% job           = structure with settings for loaded TAfT object,
%               specifically:
% .dBconversion = Boolean, whether to convert outcome variable (EEG data) to decibel (10*log10(Y)) or not.
% .EEGdomain 	= string, either 'TF' (data in TF domain) or 'time' (data in time domain).
% .HRFtype      = string, whether to perform estimation of HRF amplitude for each trial separately ('trial') or in a GLM featuring all trials of a given block ('block').
% .layout       = string, cap to use for topoplots.
% .lockSettings	= string, type of event-locked data to include, either 'stimlocked' or 'resplocked'.
% .nFreq        = numeric scalar, number of frequencies used in TF decomposition.
% .regNames 	= cell, one or more string(s) of all regressors (fMRI and behavior) to include in design matrix plus selected trials.
% .regType      = string, type of regression over trials performed, either 'pinv', 'fitlm', 'robustfit'
% .TFtype       = string, typoe of TF decomposition used in TF decomposition, either 'hanning' or 'morlet'.
% .trialdur 	= numeric scalar, trial duration to use when epoching upsampled BOLD data in seconds (recommended by Hauser et al. 2015: 8 seconds).
% dirs          = structure with relevant directory paths, especially for
%               topoplots and TFplots.
% betas         = cell, Fieldtrip object per subject, with 3D (regressor/ 
% channel/ time) or 4D (regressor/ channel/ frequency/ time) matrix with 
% regression weights of given fMRI/ behavioral regressor on EEG data.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Adapted from Tobias Hauser (https://github.com/tuhauser/TAfT).
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/CueLockedAnalyses/TAfT/

if nargin < 3
	behav2use = {};
end

if nargin < 4
	selTrials = 'all';
end

%% Load betas:

set(0, 'DefaultFigureColormap', jet(64)) % set color bar  

% Directories:
dirs.root   	= taft_set_rootdir(); % /project/3017042.02
dirs.log        = fullfile(dirs.root,'Log');
dirs.TAfT       = fullfile(dirs.log,'EEG/CueLockedResults/TAfT_Betas');
dirs.topoplot   = fullfile(dirs.TAfT,'TAfT_currentPlots');
dirs.TFplot     = fullfile(dirs.TAfT,'TAfT_currentPlots');

%% Add Fieldtrip and SPM:

addpath /home/common/matlab/fieldtrip;
rmpath /home/common/matlab/spm12;
ft_defaults;

%% Job settings:

if ~strcmp(EEGdomain,'TF') && ~strcmp(EEGdomain,'time'); error('Invalid EEGdomain specified'); end

% EEG:
job.EEGdomain       = EEGdomain; % time TF
job.nFreqs          = 15;   % number frequencies to be expected in EEG data
job.TFtype          = 'hanning'; % type of TF decomposition, either hanning or morlet
job.dBconversion    = false; % decibel conversion or not

job.lockSettings    = 'stimlocked'; % stimlocked resplocked

% fMRI:
job.HRFtype         = 'trial'; % trial or block
job.trialdur        = 8;

% Regression:
job.regType         = 'pinv'; % 'fitlm', 'robustfit'

% Update names of all regressors:
if ischar(ROIs2use) % if input only string, not cell
    ROIs2use = {ROIs2use};
end
job.regNames        = ROIs2use;
if exist('behav2use','var')
    job.regNames = [job.regNames behav2use];
end
if ~strcmp(selTrials,'all')
    job.regNames = [job.regNames selTrials];
end
% End of variable selection

% Enrich file name:
fileName         = sprintf('%s_%s_%s_%s_%dsec_%s',...
    job.EEGdomain,job.lockSettings,strjoin(job.regNames,'_'),job.HRFtype,job.trialdur,job.regType);

% Add type of TF decomposition and # frequencies:
if strcmp(job.EEGdomain,'TF')
    fileName         = sprintf('%s_freq%02d',fileName,job.nFreqs); % add number of frequencies
end

% Add if decibel conversion:
if job.dBconversion
    fileName     = sprintf('%s_dB',fileName);
end

% Specify full output directory:
fprintf('No existing job in workspace, load file \nTAfT_%s.mat\n',fileName)
fileNameFull = fullfile(dirs.TAfT,sprintf('TAfT_%s.mat',fileName));

%% Load file:

if exist(fileNameFull,'file')
    load(fileNameFull)

    % Special settings for plotting:
    job.layout   = 'easycapM11.mat';

    % Update on older files:
    if ~isfield(job,'regNames') 
        job.regNames = ROIs2use;
    end
    if ~isfield(job,'lockSettings') 
        job.lockSettings = sprintf('%slocked',job.lock);
    end
    if ~isfield(job,'EEGdomain')
        job.EEGdomain = job.MEEGdomain;
        job.EEGfile = job.MEEGfile;
    end
    fprintf('Regressors are \n%s\n',strjoin(job.regNames,'\n'))
    fprintf('Done :-) \n')

else
    error('Unable to load file: File does not exist');
end
    
end % end of function.
