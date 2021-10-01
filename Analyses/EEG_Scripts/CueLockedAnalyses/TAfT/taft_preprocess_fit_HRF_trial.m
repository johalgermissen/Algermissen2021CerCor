function betas = taft_preprocess_fit_HRF_trial(job,data)

% betas = taft_preprocess_fit_HRF_trial(job,data)
% 
% Take upsampled data epoched into trials, fit separate HRF per trial.
%
% INPUTS:
% job           = structure with settings for creating TAfT object, specifically:
% .add_intercept= whether to add intercept to design matrix (true) or not (false).
% .demeanX      = whether to demean design matrix (true) or not (false).
% .demeanY      = whether to demean trial-wise BOLD data (true) or not (false).
% .TR_ups       = new repetition time after upsampling (in seconds).
% data          =  2D (trials x samples within trial) matrix with upsampled data.
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

% GLM has only one predictor (HRF) (and optional intercept), thus only use
% respective part of SPM HRF template:

% Initialize # samples x 1 matrix:
dm      = hrf(1:size(data,2));

% Demean design matrix:
if job.demeanX
    dm = dm - nanmean(dm,1);
end

% Demean data:
if job.demeanY
    data = data - nanmean(data,2);
end

% Add intercept if requested:
if job.add_intercept  
    dm(:,end+1) = 1; % add intercept to remove offset-effects
end

%% Fit GLM for each trial separately:

% Pseudo inverse:
pdm     = taft_pinv(dm);

% Beta estimates:
coefs   = pdm*data';    % apply inverted design matrix to data; receive beta weight for each trial.
betas   = coefs(1,:);   % ignore intercept (2nd column)

end % end of function.