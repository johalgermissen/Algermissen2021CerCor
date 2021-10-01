function output = taft_postprocess_FisherZ(input)

% output = taft_postprocess_FisherZ(input)
% 
% Convert betas stored in Fieldtrip object for multiple subjects into 
% z-values using the Fisher z-transformation (Matlab's atanh)
% In my experience, this transformation hardly matters for numerically
% small values...
%
% INPUTS:
% input     = cell with Fieldtrip object for each subject, either in 
%           time (.avg) or frequency (.freq) domain.
% OUTPUTS:
% output    = same object, but .avg or .freq transformed.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Adapted from Tobias Hauser (https://github.com/tuhauser/TAfT).
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/CueLockedAnalyses/TAfT/

fprintf('Perform Fisher z-transformation on betas \n');

%% Initialize:

output  = cell(1,length(input)); % initialize empty objects
nSubs   = length(input); % number subjects

%% Loop over subjects: 

for iSub = 1:nSubs
    
    if isfield(input{iSub},'powspctrm') % if time-frequency domain
        output{iSub}            = input{iSub}; % initialize object
        output{iSub}.powspctrm  = atanh(input{iSub}.powspctrm); % Fisher z-transformation, which is inverse hyperbolic tangent
        
    elseif isfield(input{iSub},'avg') % if time domain
        output{iSub}            = input{iSub}; % initialize object
        output{iSub}.avg        = atanh(input{iSub}.avg); % Fisher z-transformation, which is inverse hyperbolic tangent
        
    else
        error('Dimensions of inputs unknown')
    end
end

end % end of function.