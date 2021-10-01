function [job,data,data_ups,i_ons] = taft_preprocess_filter_upsample_epoch(job,iROI,iBlock)

% [job,data,data_ups,i_ons] = taft_preprocess_filter_upsample_epoch(job,iROI,iBlock)
% 
% For given ROI, for given block, upsample volume-by-volume time series, 
% epoch it into trials.
% 
% INPUTS:
% job           = structure with settings for creating TAfT object, specifically:
% .hpass        = numeric scalar, highpass filter cutoff (in seconds), applied if ~=0.
% .ons_unit     = numeric scalar, whether trial onsets are given in seconds (==1) or just as volume indices (==2).
% .ROIs 		= further settings/file names to use for particular block of particular subject, specifically:
%   .ROIdef.rawfMRIfile = string, file from which to load volume-by-volume MRI data.
%   .ROIdef.onsets  	= string, file from which to obtain trial onsets.
%   .ROIdef.nuisance 	= string, file from which to load volume-by-volume nuisance regressors from.
% .TR           = numeric scalar, repetition time of MRI sequence (used for upsampling volumes) in seconds.
% .trialdur     = numeric scalar, the desired trial duration (in seconds), used for epoching.
% .ups          = numeric scalar, upsampling factor, ratio between repetition time before and after upsampling.
% iROI          = numeric scalar, index of selected ROI.
% iBlock        = numeric scalar, index of selected block.
%
% OUTPUTS:
% job           = structure with settings for creating TAfT object, new field added:
% .TR_ups       = new repetition time after upsampling (in seconds).
% data          = 2D (trials x samples within trial) matrix with upsampled data.
% data_ups      = numeric vector, long (# data points after upsampling) vector of data after upsampling.
% i_ons         = numeric vectors, trial onsets (in seconds).
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Adapted from Tobias Hauser (https://github.com/tuhauser/TAfT).
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/CueLockedAnalyses/TAfT/

%% Load volume-by-volume fMRI data (txt-file):

fprintf('Load fMRI data\n');
dat = load(job.ROIs(iROI).ROIdef(iBlock).rawfMRIfile{:}); % volume-by-volume data in column vector

%% High-pass filter:

if job.hpass ~= 0
    
    fprintf('Perform high-pass filtering at %d s\n',job.hpass);
    [b,a]   = butter(1,(1/job.hpass/(1/job.TR/2)),'high');
    dat     = filter(b,a,dat);
    
end

%% Regress out nuissance regressors:

if isfield(job.ROIs(iROI).ROIdef(iBlock),'nuisance')
    
    fprintf('Found nuisance regressor directory, regress out \n')
    dm          = load(job.ROIs(iROI).ROIdef(iBlock).nuisance{:}); % load nuisance regressors
    dm(:,end+1) = 1; % add intercept
    pdm         = taft_pinv(dm); % compute pseudo-inverse
    betas       = pdm*dat;  % compute betas
    pred_dat    = dm*betas; % get predicted signal
    res         = dat - pred_dat; % residuals (actual minus predicted signal)
    dat         = res; % save residuals instead of data

else 
    warning('No nuisance regressors found, no correction\n')
end

%% Compute onsets of trials for epoching:

% Load onsets:
rawOnsets = load(job.ROIs(iROI).ROIdef(iBlock).onsets{:});

if job.ons_unit==1          % if volume onsets already given in seconds:
    t_ons = rawOnsets;

elseif job.ons_unit==2      % if just volume indices given:
    t_ons = rawOnsets .* job.TR; % convert number of repetition times into seconds

else
    error('Unknown onset units');
end

%% Upsample data:

job.TR_ups = job.TR/job.ups; % length of data samples after interpolation (in s)

fprintf('Upsample data by factor %d using spline: new TR = %.03f sec.\n',job.ups,job.TR_ups);
data_ups     = spline(1:length(dat),dat,1:(1/job.ups):length(dat));

% Provide interpolated data points at data sides from first until final volume, in
% steps of 1/upsampling factor (how many new samples between original samples)

%% Compute trial timing and duration in samples:

t_ups_ons   = 0:job.TR_ups:(length(data_ups)-1)*job.TR_ups; % timings of upsampled data samples (in s)
i_ons       = taft_findc(t_ups_ons,t_ons); % find timing of upsampled samples closest to trial onsets, i.e. upsampled trial onsets
dur_hrf_ups = round(job.trialdur/job.TR_ups); % number upsampled samples per trial (given trialdur)

%% Epoch trials:

fprintf('Epoch data into trials of %.2f sec. duration\n',job.trialdur);

% Create t_ind matrix (# trials x # time points) that contains start and
% end point of each trial; apply to data_ups:

t_ind = nan(length(i_ons),(dur_hrf_ups+1)); % initialize t_ind

for t = 1:length(i_ons) % for each trial: retrieve indices of data points to be included
    t_ind(t,:) = i_ons(t):(i_ons(t)+dur_hrf_ups); % from new trial onset until trial onset + # samples per trial
end

data = data_ups(t_ind); % retrieve from upsampled data: rows are trials, columns time within trial

end % end of function.