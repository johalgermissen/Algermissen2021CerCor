function out = redefine_resplock(iSub,sRate,validChan,behav,data)

% out = redefine_resplock(iSub,sRate,validChan,,behav,data)
% 
% Re-epochs stimulus-locked EEG data (time domain) based on a vector of
% reaction times.
%
% INPUT:
% iSub          = scalar integer, subject number
% sRate         = scalar integer, sampling rate
% validChan     = vector of integers, indices of valid channels
% behav         = cell with stimulus identity (.stim, .winCue) and RTs (.RT)
% data          = Fieldtrip structure with EEG data (time-domain)
%
% OUTPUT:
% out           response-locked behavioral data
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

%% Extract RTs from input:

RT2use  = behav.RT'; % RTs for each trial (all 640 trials) in EEG data.

%% Delete too long RTs:

RTthreshold = 1.3035; % maximal possible RT for response "in time"
fprintf('Subject %03d: Found %d trials with too long RTs, set to NaN\n',...
    iSub,sum(RT2use>RTthreshold));
RT2use(RT2use > RTthreshold) = NaN;

%% Impute RTs for NoGos:

% RTs for NoGo trials are NaNs, so impute them by mean of Go trials per cue valence:

% Retrieve mean RT for Win and Avoid cues:
RTwin       = nanmean(behav.RT(ismember(behav.stim,behav.winCue))); % mean RT for Win cues.
RTavoid     = nanmean(behav.RT(~ismember(behav.stim,behav.winCue))); % mean RT for Avoid cues.

RTwin       = round(RTwin,3); % round to EEG sampling rate (1 ms)
RTavoid     = round(RTwin,3); % round to EEG sampling rate (1 ms)

% Impute missing RT values per Win and Avoid with mean values:
RT2use(ismember(behav.stim,behav.winCue) & isnan(RT2use))  = RTwin; % for NoGo trials, use mean RT for Win cues.
RT2use(~ismember(behav.stim,behav.winCue) & isnan(RT2use)) = RTavoid; % for NoGo trials, use mean RT for Win cues.

% Check:
if sum(isnan(RT2use)) > 0
    error('Still remaining NaNs in RT even after imputation')
end

%% Convert seconds to ms samples (how many ms samples relative to current t=0):

RT2useSamples  = RT2use .* sRate;

%% Redefine trials to response-locked:

% Delete non-EEG channels:
if isfield(data,'elec') && length(data.elec.elecpos) > length(validChan)
    data.elec.elecpos    = data.elec.elecpos(validChan,:); % filter out non-EEG channels
end

% Redefine trial (RT becomes 0):
fprintf('Subject %03d: Redefine trial \n',iSub)
cfg                 = [];
cfg.offset          = -RT2useSamples; % vector with begin samples relative to t=0
out                 = ft_redefinetrial(cfg,data);

%% Check new start and end points of trials: too late/ too early?

startVec    = nan(1,length(out.trial));
endVec      = nan(1,length(out.trial));

for iTrial = 1:length(out.trial)

    startVec(iTrial) = out.time{iTrial}(1); % Retrieve start of each trial
    endVec(iTrial) = out.time{iTrial}(end); % Retrieve end of each trial
%        fprintf('Trial %d: start %.03f, end %.03f\n',iTrial,startVec(iTrial),endVec(iTrial))

end

% Cutoffs hard coded, but seem reasonable:
if sum(startVec > -1) > 0
    error('Found trials starting later than -1 sec.')
elseif sum(endVec < 0.2) > 0
    error('Found trials stopping earlier than 0.2 sec.')
end

end % end of function.