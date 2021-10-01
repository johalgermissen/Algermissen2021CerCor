function [sortBetas,Tvalues] = taft_postprocess_TF_selectData(job,betas,iROI)

% Extracts time-frequency data for selected ROI,
% sort channels across subjects, selects valid subjects,
% combine subjects into 4D data, perform one-sample t-test for each
% time-frequency-channel bin across subjects.
%
% Use as 
%   [sortBetas,Tvalues] = taft_postprocess_TF_selectData(job,betas,iROI)
%
% INPUTS:
% job           = cell, necessary settings for subject selection:
% .invalidSubs 	= numeric vectors, subjects to exclude from analyses (optional, will otherwise fill in subjects excluded in analyses reported in paper).
% .regNames 	= cell, one or more string(s) of all regressors (fMRI and behavior) to include in design matrix plus selected trials.
% betas         = cell, Fieldtrip object for each subject, contains
% regression weights for each ROI/channel/frequency/time bin in .powspctrm.
% iROI          = numeric scalar, number of ROI to be extracted.
%
% OUTPUTS:
% sortBetas     = cell, Fieldtrip object for each subject, contains 
% regression weights for each ROI/channel/frequency/time bin in .powspctrm, 
% channels sorted for each subject, for selected subjects.
% Tvalues       = Fieldtrip object with t-value for each
% channel/frequency/time bin in .powspctrm.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Adapted from Tobias Hauser (https://github.com/tuhauser/TAfT).
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/CueLockedAnalyses/TAfT/

%% Complete settings:

if nargin < 3
    iROI = 1;
    fprintf('No ROI specified -- use first ROI\n')
end

fprintf('Extract betas for selected ROI: %s\n',char(job.regNames(iROI)))

%% Initialize:

nSub = length(betas); % Number of subjects
ROIBetas = cell(1,nSub); % initialize cell for only one ROI

%% Loop over subjects, extract only selected ROI:

for iSub = 1:nSub % iSub = 1;
    
    ROIBetas{iSub}              = betas{iSub}; % reassign to new object
    ROIBetas{iSub}.powspctrm    = squeeze(ROIBetas{iSub}.powspctrm(iROI,:,:,:)); % extract ROI (in first dimension)

end

%% Sort channels across subjects:

sortBetas       = cell(1,nSub);
labelTemplate   = ROIBetas{1}.label; % take subject 1 as reference

for iSub = 1:nSub % iSub = 1
    
    sortBetas{iSub}             = ROIBetas{iSub}; % copy template data
    [~,labIdx2]                 = match_str(labelTemplate,ROIBetas{iSub}.label); % indices of labels
    sortBetas{iSub}.label       = ROIBetas{iSub}.label(labIdx2); % sort labels
    sortBetas{iSub}.powspctrm   = sortBetas{iSub}.powspctrm(labIdx2,:,:); % sort data 
    
end

%% Check label consistency across subjects:

for iSub = 1:nSub
    % Return warnings if labels of certain subjects don't match template
    if sum(strcmp(sortBetas{iSub}.label,labelTemplate)) < length(labelTemplate)
        warning(sprintf('Mismatch in labels: Inspect subject %d\n',iSub)) %
    end

end

%% Select subjects: 

if ~isfield(job,'invalidSubs')
%     job.invalidSubs = []; % no invalid subs (for plotting)
    % job.invalidSubs = [15 25]; % subjects with bad co-registrations to start with
    job.invalidSubs = [1 11 15 19 21 25 26]; % invalid co-registrations and outliers
end

job.validSubs   = setdiff([1:nSub],job.invalidSubs);
fprintf('Exclude subjects %s\n',strjoin(string(job.invalidSubs),', '))

% Exclude subjects from data object:
sortBetas       = sortBetas(job.validSubs);

%% Reshape into 4D (subject as first dimension):

fprintf('Reshape into 4D\n')

c = zeros(length(sortBetas),length(sortBetas{1}.label),length(sortBetas{1}.freq),length(sortBetas{1}.time)); %
for iSub = 1:length(sortBetas) % iSub = 1;
    c(iSub,:,:,:) = sortBetas{iSub}.powspctrm; % bring into 4D, with subject as first dimension
end

%% Perform T-test across all channels:

fprintf('Compute T-test across subjects\n')

% Reshape into 2D:
nS = size(c,1); % number of 'subjects' (recording sites)
nR = size(c,2); % number of regressors (e.g. channels -- anything spatially unclustered)
nF = size(c,3); % number of frequencies
nT = size(c,4); % number of timebins
cr = reshape(c,nS,nR*nF*nT); % concatenate regressors/frequencies/timebins into 2D

% One-sample t-test:
dm          = ones(nS,1); % design matirx: mean of all subjects
[~,~,tg]    = ols(cr,dm,eye(size(dm,2))); % compute t-value for each data point
tg          = reshape(tg,size(dm,2),nR,nF,nT); % reshape back into 4D

% Reshape back into Fieldtrip object:
cfg                 = [];
Tvalues             = ft_freqgrandaverage(cfg,sortBetas{:}); % template
Tvalues.powspctrm   = squeeze(tg); % assign back

end % end of function.