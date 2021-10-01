function pX = taft_pinv(X)

% pX = taft_pinv(X)
%
% Compute pseudo-inverse using SPM's spm_sp function.
% Check for rank deficiency.
%
% INPUTS:
% X         = 2D (trials x regressors) design matrix.
% OUTPUTS:
% pX        = 2D pseudo-inverse calculated with SPM.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Adapted from Tobias Hauser (https://github.com/tuhauser/TAfT).
% Should work in Matlab 2018b.

%% Initialize design matrix:

nR = size(X,2); % # columns of design matrix 
x  = spm_sp('Set',X); % Initialize design matrix as SPM object

%% Check for rank deficiency:

if x.rk~=nR % Compare # columns of manually created and SPM design matrix
    error('Rank deficient design matrix (after removal of bad trials)');
end

%% Actually calculate design matrix using the SPM function:

pX = spm_sp('x-',x); % Compute pseudoinverse of design matrix

end % end of function.