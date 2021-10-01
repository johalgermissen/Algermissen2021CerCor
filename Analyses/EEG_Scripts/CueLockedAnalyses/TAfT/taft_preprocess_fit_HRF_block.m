function betas = taft_preprocess_fit_HRF_block(job,data_ups,i_ons)

% betas = taft_preprocess_fit_HRF_trial(job,data_ups,i_ons)
% 
% Take upsampled data of entire block (one long vector), fit HRF amplitude
% for every trial within block with one single design matrix for entire
% block.
%
% INPUTS:
% job           = structure with settings for creating TAfT object, specifically:
% .add_intercept= whether to add intercept to design matrix (true) or not (false).
% .demeanX      = whether to demean design matrix (true) or not (false).
% .demeanY      = whether to demean trial-wise BOLD data (true) or not (false).
% .TR_ups       = new repetition time after upsampling (in seconds).
% data_ups      = numeric vector, long (# data points after upsampling) vector of data after upsampling.
% i_ons         = numeric vectors, trial onsets (in seconds).
%
% Output:
% betas         = numeric vector, HRF amplitude for each trial.
%
% OUTPUTS:
% betas         = numeric vector, HRF amplitude for each trial.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Adapted from Tobias Hauser (https://github.com/tuhauser/TAfT).
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/CueLockedAnalyses/TAfT/

%% Obtain canonical HRF for specific RT:

hrf     = spm_hrf(job.TR_ups); % retrieve SPM HRF template given resolution (i.e. upsampled repetition time).

%% Prepare design matrix based on HRF:

% Initialize # samples x # trials matrix:
dm      = zeros(length(data_ups),length(i_ons)); % initialize empty design matrix

for iTrial = 1:length(i_ons) % for each trial:
    dm(i_ons(iTrial):(i_ons(iTrial)+length(hrf)-1),iTrial) = hrf; % insert HRF template
end

dm      = dm(1:length(data_ups),:); % crop if HRF overlength at the end

% Demean design matrix:
if job.demeanX
    dm  = dm - nanmean(dm,1);
end

% Demean data:
if job.demeanY
    Y   = data_ups - nanmean(data_ups);
else
    Y   = data_ups;
end

% Add intercept if demanded:
if job.add_intercept  
    dm(:,end+1) = 1; % add intercept to remove offset-effects
end

%% Fit GLM across all trials of entire block:

% Compute pseudo-inverse:
pdm     = taft_pinv(dm);

% Compute beta estimates:
betas   = (pdm*Y')';  % apply inverted design matrix to data
betas   = betas(1:length(i_ons)); % drop intercept at the end if necessary

end % end of function.