function cordata = subtract_ERP(cordata,ERPdata,corselTrials,ERPselTrials,mode,subtractTime,ERPdur)

% cordata = subtract_ERP(cordata,ERPdata,corselTrials,ERPselTrials,mode,subtractTime,ERPdur)
% 
% Computes the mean ERP given input data ERPdata for selected trials
% ERPselTrials, correct input data cordata for this mean ERP for selected
% trials corselTrials.
% ERP computation and ERP correction allow for separate data sets and
% trials. For example, it is possible to compute a mean response-locked ERP
% and subtract it from each stimulus-locked trials at the respective RT.
% ERP computation and ERP correction can be applied when trials are
% misaligned in time (mode 'Manual') or alligned in time (mode
% 'Fieldtrip').
%
% INPUT:
% cordata       = trial-by-trial data to be corrected. Trials can very in
% timing (e.g. after response-locking data that were previously
% stimulus-locked).
% ERPdata       = trial-by-trial data used for computing ERP (e.g. subset
% of trials). It is assumed that trials are aligned in time.
% corselTrials  = vector with trials to be corrected (default: all trials)
% ERPselTrials  = vector with trials used for ERP computation (default: all trials)
% mode          = mode of correction, either 'Manual' or 'Fieldtrip'
% subtractTime  = in mode 'Manual': vector with start of correction window for each trial in cordata (default: each trial starts at 0)
% ERPdur        = in mode 'Manual': scalar duration of meanERP to be computed based on ERPdata (start at 0, default duration 0.7 sec.)
%
% OUTPUT:
% cordata       = trial-by-trial data corrected for ERP
%
% Recommendation:
% 'Manual': Useful if timing of to-be-computed ERP differs per trial,
% e.g. when data are response-locked and Fieldtrip objects still contains
% old trial timings.
% 'Fieldtrip': Useful as much simpler version if all trials are aligned
% in time.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

%% Fill in settings:

% Indices of selected trials:
if nargin < 2
    ERPdata = cordata;
    fprintf('Use same data for computing ERP and correction\n')
end

if nargin < 3
    corselTrials = 1:length(cordata.trial);
    fprintf('No selected trials to-be-corrected specified -- use all trials\n')
end

if nargin < 4
    ERPselTrials = 1:length(ERPdata.trial);
    fprintf('No selected trials for ERP computation specified -- use all trials\n')
end

if nargin < 5
    mode = 'Manual';
    fprintf('No mode specified---use manual correction')
end

if nargin < 6 && strcmp(mode,'Manual') 
    subtractTime = zeros(1,length(cordata.trial));
    fprintf('No trialwise timings for start of correction window specified -- use 0 for all trials\n')
end

if nargin < 7 && strcmp(mode,'Manual') % no ERPdur needed for Fieldtrip mode
    ERPdur = 0.7;
    fprintf('No duration of ERP computation specified -- use %0.1f sec.\n',ERPdur);
end

%% Part I: Compute mean ERP

fprintf('Compute ERP based on %d trials \n',length(ERPselTrials))

%% Alternative A: Extract ERP manually in specified window:

if strcmp(mode,'Manual')

    ERPstart = 0; % assume that all trials in ERPdata are aligned, so start at 0 to compute ERP
    
    % ------------------------------------------------------------------- %
    % Initialize channel x time x trial matrix as template for ERP:
    if iscell(ERPdata.time)
        meanERP = nan(length(ERPdata.label),length(dsearchn(ERPdata.time{1}',ERPstart):dsearchn(ERPdata.time{1}',ERPdur)),length(ERPselTrials));
    else
        meanERP = nan(length(ERPdata.label),length(dsearchn(ERPdata.time',ERPstart):dsearchn(ERPdata.time',ERPdur)),length(ERPselTrials));
    end
    
    % ------------------------------------------------------------------- %
    % Loop over trials, extract relevant data to compute average ERP:    
    for iTrial = 1:length(ERPselTrials) % loop through selected trials
        
        iTrialIdx = ERPselTrials(iTrial); % relative index of this trial within vector of all trials
        
        % Each trial in one cell, time line can differ between trials (e.g.
        % because response-locked after previously stimulus-locked)
        if iscell(ERPdata.time) 
            startIdx = dsearchn(round(ERPdata.time{iTrialIdx},3)',ERPstart); % start point of ERP computation for this trial
            stopIdx = dsearchn(round(ERPdata.time{iTrialIdx},3)',ERPstart + ERPdur); % end point of ERP computation for this trial
            meanERP(:,:,iTrial) = cordata.trial{iTrialIdx}(:,startIdx:stopIdx); % extract data in specified window
            
        % All trials in one single 3D matrix (rather than different fields
        % in cell), more convenient, requires that all trials are
        % time-aligned (e.g. if stimulus-locked)
        else 
            startIdx = dsearchn(round(ERPdata.time{iTrialIdx},3)',ERPstart); % start point of ERP computation for this trial
            stopIdx = dsearchn(round(ERPdata.time{iTrialIdx},3)',ERPstart + ERPdur); % end point of ERP computation for this trial
            meanERP(:,:,iTrial) = cordata.trial(iTrialIdx,:,startIdx:stopIdx); % extract data in specified window

        end
    end
    
    % Average over trials to get ERP:
    meanERP = nanmean(meanERP,3); % average over trials
    meanERP = meanERP - meanERP(:,1); % demean vertically (ERPs starts at 0 at time index 1)

%% Alternative B: average with Fieldtrip ft_timelockanalysis
% (assumes that all trials are temporally aligned/ timing information
% provided)

elseif strcmp(mode,'Fieldtrip')
    
    cfg = [];
    cfg.trials = ERPselTrials;
    tmp = ft_timelockanalysis(cfg,ERPdata); % define for this condition
    meanERP = tmp.avg;
    meanERP = meanERP - meanERP(:,1); % demean vertically (ERPs starts at 0 at time index 1)

else
    error('Unknown mode')
    
end

%% Part II: Subtract ERP:

fprintf('Compute ERP from selected trials based on %s mode\n',mode)

%% Alternative A: Manually create trial-by-trial template (potentially different trials timings) and subtract

if strcmp(mode,'Manual')
    
    for iTrial = 1:length(corselTrials) % loop through selected trials only (leave others unaltered); iTrial = 1
        
        iTrialIdx = ERPselTrials(iTrial); % relative index of this trial within vector of all trials
        
        % Timings:
        startTime   = subtractTime(iTrial); % retrieve start of correction window
        stopTime    = startTime + ERPdur;   % end of correction window
        
        % Each trial in one cell, time line can differ between trials (e.g.
        % because response-locked after previously stimulus-locked)
        if iscell(cordata.time) 
            
            % Start and stop of ERP per trial:
            startIdx = dsearchn(round(cordata.time{iTrialIdx},3)',startTime); % start of time window to subtract from
            stopIdx = dsearchn(round(cordata.time{iTrialIdx},3)',stopTime); % end of time window to subtract from
            
            % Initialize empty template based on start and stop time:
            trialERP = zeros(size(cordata.trial{iTrialIdx},1),size(cordata.trial{iTrialIdx},2)); % empty trial template
            
            % Fill template in specified time window with mean ERP:
            trialERP(:,startIdx:stopIdx) = meanERP + cordata.trial{iTrialIdx}(:,startIdx); % insert ERP; vertically align template with trial (intercept at startIdx)
            
            % Subtract mean ERP template from cordata of this trial:
            cordata.trial{iTrialIdx} = cordata.trial{iTrialIdx} - trialERP; % subtract template
            
        % All trials in one single 3D matrix (rather than different fields
        % in cell), more convenient, requires that all trials are
        % time-aligned (e.g. if stimulus-locked)
        else 
            
            % Start and stop of ERP per trial:
            startIdx = dsearchn(round(cordata.time,3)',startTime); % start of time window to subtract from
            stopIdx = dsearchn(round(cordata.time,3)',stopTime); % end of time window to subtract from
            
            % Initiate empty template based on start and stop time:
            trialERP = zeros(size(cordata.trial,2),size(cordata.trial,3)); % empty trial template

            % Fill template in specified time window with mean ERP:
            trialERP(:,startIdx:stopIdx) = meanERP + cordata.trial(iTrialIdx,:,startIdx)'; % insert ERP; vertically align template with trial (intercept at startIdx)
            
            % Subtract mean ERP template from cordata of this trial:
            cordata.trial(iTrialIdx,:,:) = squeeze(cordata.trial(iTrialIdx,:,:)) - trialERP; % subtract template

        end
    end
    
%% Alternative B: if mode Fieldtrip (so trials temporally alligned): subtract from entire matrix
    
elseif strcmp(mode,'Fieldtrip')

    for iTrial = 1:length(corselTrials) % loop through selected trials only (leave others unaltered); iTrial = 1

        iTrialIdx = corselTrials(iTrial); % retrieve relative trial index
        
        % Vertical alignment unnecessary because no restricted trial window used
        cordata.trial(iTrialIdx,:,:) = squeeze(cordata.trial(iTrialIdx,:,:)) - meanERP; % subtract from each trial - meanERP already in right format

    end
    
else
    error('Unknown mode')
end

end % end of function.
