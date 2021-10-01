function difBetas = taft_postprocess_subtract(betas1, betas2)

% difBetas = taft_postprocess_subtract(betas1, betas2)
%
% Loop over subjects, subtract one ROI from the other.
%
% INPUTS:
% betas1,betas2 = cell with Fieldtrip object for each subject, contains
%               regression weights for each ROI/channel(/frequency)/time bin
%               in .avg or .powspctrm.
% OUTPUTS:
% difBetas      = cell with Fieldtrip object for each subjects, difference
%               between both ROIs.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Adapted from Tobias Hauser (https://github.com/tuhauser/TAfT).
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/CueLockedAnalyses/TAfT/

fprintf('Subtract inputs (1st input - 2nd input) from each other\n');

%% Initialize:

difBetas    = cell(1,length(betas1)); % initialize empty objects
nSubs       = length(betas1);

%% Loop over subjects:

for iSub = 1:nSubs % loop over subjects:
    
    if isfield(betas1{iSub},'powspctrm') && isfield(betas2{iSub},'powspctrm') % if time-frequency domain
        difBetas{iSub}              = betas1{iSub}; % initialize object
        difBetas{iSub}.powspctrm    = betas1{iSub}.powspctrm - betas2{iSub}.powspctrm; % subtract second from first input

    elseif isfield(betas1{iSub},'avg') && isfield(betas2{iSub},'avg') % if time domain
        difBetas{iSub}              = betas1{iSub}; % initialize object
        difBetas{iSub}.avg          = betas1{iSub}.avg - betas2{iSub}.avg; % subtract second from first input
    
    else
        error('Dimensions of both inputs unknown or not-matching')
    end
end

end % end of function.
