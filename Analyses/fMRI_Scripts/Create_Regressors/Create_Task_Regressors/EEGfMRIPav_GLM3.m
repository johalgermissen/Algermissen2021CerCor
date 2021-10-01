function EEGfMRIPav_GLM3()

% EEGfMRIPav_GLM3
% 
% Create regressors for subject-level GLM (concatenated blocks).
% Allows for contrasts valence, required action, performed action, 
% conflict (8 regressors), stimulus-locked midfrontal alpha power,
% response-locked midfrontal theta power.
% Use three-column (onset, duration, factor level) format.
% Saves each regressor as separate .txt file in respective
% timings_regressor folder of respective GLM of respective subject.
%
% Regressors of interest are:
% Valence (Win/Avoid) x Required Action (Go/NoGo) x Actual Action
% (Go/NoGo):
% 1) Win_ReqGo_ActGo
% 2) Win_ReqGo_ActNoGo
% 3) Win_ReqNoGo_ActGo
% 4) Win_ReqNoGo_ActNoGo
% 5) Avoid_ReqGo_ActGo
% 6) Avoid_ReqGo_ActNoGo
% 7) Avoid_ReqNoGo_ActGo
% 8) Avoid_ReqNoGo_ActNoGo
%
% Covariates of no interest are:
% 9) Handedness of response: hand (left, NoGo, right)
% 10) Performance: Error (correct/incorrect)
% 11) Outcome Onset
% 12) Outcome Valence
% 13) Invalid Outcome
%
% EEG:
% 14) EEG alpha (conflict contrast, stimulus-locked) over Fz/FCz/Cz
% 15) EEG theta (action contrast, response-locked) over Fz/FCz/Cz
%
% Mind setting root directory dirs.root
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

% clear all; clc

%% Set directories:

dirs.root       = '/project/3017042.02'; % root directory--needs to be adapted to users' folder structure
dirs.behav      = fullfile(dirs.root,'Log/Behavior');
dirs.raw        = fullfile(dirs.behav,'Data_beh_mat');

%% Fixed settings:

nSub = 36;

GLMID       = '3';
nRegressors = 15; % for specification of empty regressors

%% Detect which files are in behavioral raw data directory:

input = dir(dirs.raw); % store details of all files there
allfiles = {input.name}; % extract names and store as cell array
allfiles = allfiles'; % tranpose
allfiles{1,1} = []; % delete '.'
allfiles{2,1} = []; % delete '..'   
dirs.subs = allfiles(~cellfun('isempty',allfiles)); % remove deleted entries

% Object to store all valid trials:
validTrialsAll  = nan(length(dirs.subs),640);

%% Loop over subjects, create and save regressors:

for iSub = 1:nSub % length(dirs.subs) % specify for which subjects to create regressors
    
    % Print start and end of each subject in case something goes wrong:
    fprintf('Subject %03d: Start \n', iSub)
 
    % Load event and timing data:
    fileName = strcat('3017042.02_emmvdij_',sprintf('%03d',iSub),'_001_results.mat'); % raw data of specific subject (mind that subject ID can be 2 or 3 digits)
    load(fullfile(dirs.raw,fileName));
    
    % Set directory where regressors will be printed:
    dirs.timings = fullfile(dirs.root,sprintf('Log/fMRI/sub-%03d/GLM%s/timings_regressors',iSub,GLMID));

    % Create directory if it doesn't exist yet: 
    if ~exist(dirs.timings ,'dir'); mkdir(dirs.timings); end
    
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    %% Retrieve data for regressors:
    
    % Just concatenate both sessions:
    stimIDAll      = [prep.seq.learn.stim{1}; prep.seq.learn.stim{2}]; % cue condition (1-8)
    respAll        = [results.learn{1}.response; results.learn{2}.response]; % button ID pressed: 69,70,71,72,101,102,103,104 for left; 65,66,67,68,97,98,99,100 for right
    valAll         = double(ismember(stimIDAll,[1 2 5 6])); % stimulus valence
    respNumAll     = NaN(640,1); % responses recoded to 1,2,3
    isGoAll        = NaN(640,1); % respoonses recoded to 1,0
    motorNumAll    = NaN(640,1); % responses recoded to 1,-1,0
    RTAll          = [results.learn{1}.RT; results.learn{2}.RT]; % RT for determining whether response was too late (if too late, always wrong)
    correspAll     = [prep.seq.learn.resp{1}; prep.seq.learn.resp{2}]; % correct response
    reqActAll      = double(ismember(correspAll,[101,97])); % correct vigor type (Go/NoGo)
    outcomeAll     = [results.learn{1}.outcome; results.learn{2}.outcome]; % outcome obtained: -1, 0, 1
    correctAll     = repelem(NaN,length(stimIDAll))'; % incorrect responses; created after responses recoded
    incorrectAll   = repelem(NaN,length(stimIDAll))'; % incorrect responses; created after responses recoded

    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    %% Recode button presses, accuracy, RTs:
    % Recode uninstructed button presses on correct button box as instructed button presses
    % Recode accuracy based on these recoded responses; too late responses always as incorrect

    for k = 1:length(stimIDAll)
        if ismember(respAll(k),[69,70,71,72,101,102,103,104]) % left Go response
                respAll(k) = prep.par.key.left; % actually instructed left Go response is 101
                respNumAll(k) = 1; % recode numerically to 1
                isGoAll(k) = 1;
                motorNumAll(k) = 1; % left hand
        elseif ismember(respAll(k),[65,66,67,68,97,98,99,100]) % right Go response
                respAll(k) = prep.par.key.right; % actually instructed right Go response is 97
                respNumAll(k) = 2; % recode numerically to 1
                isGoAll(k) = 1;
                motorNumAll(k) = -1; % right hand
        else % NoGo response
                respAll(k) = 0; % no response = NoGo; implicitly also sets NAs to NoGo
                respNumAll(k) = 3; % recode numerically to 1
                isGoAll(k) = 0;
                motorNumAll(k) = 0; % NoGo
        end
        % Accuracy based on recoded responses
        if respAll(k) ~= correspAll(k) || RTAll(k) > 1.3035 % wrong key or too late (1.3035 latest RT) --> incorrect
            incorrectAll(k) = 1; % 1 if incorrect
            correctAll(k) = 0;
        else
            incorrectAll(k) = 0; % 0 if correct
            correctAll(k) = 1;
        end
        % Replace NAs in outcome with neutral outcome (in line with stan model):
        % if isnan(outcomeAll(k)) == 1;  
        %    outcomeAll(k) = 0;        
        % end
    end
    % Don't recode too short or too long RTs because motor response still given
    % Don't recode trials with too long RTs as NoGos (and accuracy respectively) because people might have thought they made a Go response

    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    %% Retrieve whether valence outcome experienced or not:
        
    % Initialize other task settings:
    nStim           = 16;
    nTrial          = 640;
    stimAll         = [stimIDAll(1:320);8+stimIDAll(321:640)];              % add 8 to stimulusID in second session of task
    cueDetected     = zeros(1,nStim);                                       % valenced outcome for cue already experienced or not
    cueDetectAll    = nan(nTrial,1);                                        % cue valence known or not per trial
    validOutcome    = nan(nTrial,1);                                       % store trials where outcome was NaN to later delete row in regressor

    for iTrial = 1:nTrial
        % Retrieve stimulus ID
        stim = stimAll(iTrial);
        % Store whether cue valence experienced in this trial or not
        cueDetectAll(iTrial) = cueDetected(stim);        
        % Retrieve actually received outcome:
        r = outcomeAll(iTrial);
        % "Unmute" stimulus if valence seen for the first time
        if abs(r) == 1; cueDetected(stim) = 1; end
        % Retrieve whether valid or invalid outcome received
        if isnan(r)==1; validOutcome(iTrial) = 0; else; validOutcome(iTrial) = 1; end
    end
    
    % Store valid trials (just for check afterwards):
    validTrialsAll(iSub,:) = validOutcome;

    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    %% Load EEG regressors:
    
    % a) Retrieve single-trial EEG time-course:
    EEGRegressor1  = load(fullfile(dirs.timings,'EEG_CongruencyCzFCzFzalphastimlocked'));
    EEGRegressor2  = load(fullfile(dirs.timings,'EEG_GoCzFCzFzthetaresplocked'));
    
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    %% Create timing for regressors:

    % Strategy for concatenated blocks:
    % a) Retrieve number of volumes per block, compute effective fMRITime in GLM as TE * position of volume
    % b) Extract when in fMRITime each blocks starts,
    % c) extract for when block started in behavioral time (set to 0),
    % d) extract timing of trials from behavioral data files and correct for when blocks started in behavioral time,
    % e) add when block started in fMRITime    
    
    % a) retrieve block lengths in volumes, create fMRI time:
    dirs.sub            = sprintf('/project/3017042.02/Log/fMRI/sub-%03d',iSub);
    cd(dirs.sub)
    nVolumes            = load(sprintf('Sub%03d_Nvolumes.txt',iSub)); % number of volumes for all 6 blocks
    TE                  = 1.4; % echo time
    fMRITime            = [0:(sum(nVolumes)-1)]*TE; % when each volume starts
    % how many volumes have occurred till end of each block, minus the block itself (so rather at the start of each block), +1 to start at index 1:
    fMRIBlockStart      = cumsum(nVolumes)-nVolumes+1; 
    
    % b) extract when each block started in fMRI time:
    trialBlockLength    = [110 110 100 110 110 100];
    trialTime0          = nan(0,0); % time of first volume per block
    for i = 1:length(trialBlockLength)
        trialTime0   = [trialTime0; repmat(fMRITime(fMRIBlockStart(i)),1,trialBlockLength(i))']; % repeat for each volume in block   
    end    
       
    % c) extract when blocks started in behavioral time:
    % extract timing at start of blocks, repeat 110 (or 100 times), 
    % put together into t0 that will be subtracted from behavioral time:
    % retrieve timing of start of first 3 blocks, repeat for each of the 110 trials:
    t0part1 = reshape(repmat(tm.learn{1}.stim([1 111 221]),1,110)',[330 1]); % one long vector with 330 elements
    t0part1 = t0part1(1:320); % drop last 10 entries as block 3 was 10 trials shorter
    % retrieve timing of start of second 3 blocks, repeat for each of the 110 trials:
    t0part2 = reshape(repmat(tm.learn{2}.stim([1 111 221]),1,110)',[330 1]);
    t0part2 = t0part2(1:320); % drop last 10 entries as block 6 was 10 trials shorter
    t0      = [t0part1; t0part2] - 10; % concatenate; now subtract 10 seconds that were used for getting MRI signal steady-state (dummy scans)
   
    % d) extract trial timing in behavioral time, correct for when block started:
    % 1) retrieve timings of cues corrected for t0 (onset of block):
    tCueAllBlock        = [tm.learn{1}.stim; tm.learn{2}.stim] - t0; % take cue onset; subtract start of each block
    % 2) retrieve timings of responses corrected for t0 (onset of block):
    tRespAllBlock       = [tm.learn{1}.response; tm.learn{2}.response] - t0; % note that in case of nogo, the cue offset time was recorded.
    % 3) retrieve timing of outcomes (warning), correct for NaNs when people got error message, corrected for t0:
    tFBPrel        = [tm.learn{1}.outcome; tm.learn{2}.outcome]; % take onset of learning, so far uncorrected
    tWarningAll    = [tm.learn{1}.warning; tm.learn{2}.warning]; % take timing of warning messages, so far uncorrected
    tFBPrel(isnan(tFBPrel)) = tWarningAll; % replace NaNs with timing of warnings
    tFBAllBlock    = tFBPrel - t0; % subtract start of each block --> now corrected
    % sum(isnan(tFBAllBlock)) % check whether NaN due to uninstructed key presses (NAs left)?
    
    % e) add start of blocks in fMRI time:
    tCueAll         = tCueAllBlock+trialTime0;
    tRespAll        = tRespAllBlock+trialTime0;
    tFBAll          = tFBAllBlock+trialTime0;
    
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    %% Create behavioral aspects of regressors:
    
    % 1) Cues:
    % Split up trials into those with Win vs. Avoid cues and Go vs. NoGo responses:
    % mind cueDetectAll: don't model cues for trials where outcome not yet experienced
    Win_ReqGo_ActGo_Trials         = (ismember(stimIDAll,[1 2 5 6]) & ismember(correspAll,[97 101]) & ismember(respNumAll,[1 2]) & cueDetectAll==1); % indicator whether trial was Win trial with Go response
    Win_ReqGo_ActNoGo_Trials       = (ismember(stimIDAll,[1 2 5 6]) & ismember(correspAll,[97 101]) & ismember(respNumAll,[3]) & cueDetectAll==1); % indicator whether trial was Win trial with Go response
    Win_ReqNoGo_ActGo_Trials       = (ismember(stimIDAll,[1 2 5 6]) & ismember(correspAll,[0]) & ismember(respNumAll,[1 2]) & cueDetectAll==1); % indicator whether trial was Win trial with Go response
    Win_ReqNoGo_ActNoGo_Trials     = (ismember(stimIDAll,[1 2 5 6]) & ismember(correspAll,[0]) & ismember(respNumAll,[3]) & cueDetectAll==1); % indicator whether trial was Win trial with Go response
    Avoid_ReqGo_ActGo_Trials       = (ismember(stimIDAll,[3 4 7 8]) & ismember(correspAll,[97 101]) & ismember(respNumAll,[1 2]) & cueDetectAll==1); % indicator whether trial was Win trial with Go response
    Avoid_ReqGo_ActNoGo_Trials     = (ismember(stimIDAll,[3 4 7 8]) & ismember(correspAll,[97 101]) & ismember(respNumAll,[3]) & cueDetectAll==1); % indicator whether trial was Win trial with Go response
    Avoid_ReqNoGo_ActGo_Trials     = (ismember(stimIDAll,[3 4 7 8]) & ismember(correspAll,[0]) & ismember(respNumAll,[1 2]) & cueDetectAll==1); % indicator whether trial was Win trial with Go response
    Avoid_ReqNoGo_ActNoGo_Trials   = (ismember(stimIDAll,[3 4 7 8]) & ismember(correspAll,[0]) & ismember(respNumAll,[3]) & cueDetectAll==1); % indicator whether trial was Win trial with Go response

    % 2) Responses:
    ErrorTrials     = incorrectAll==1;
    
    % 3) Feedback:
    OutcomeTrials = validOutcome==1; 
    InvalidTrials   = validOutcome==0;  
    
    emptyRegressorMatrix = zeros(1,nRegressors); % initialize empty regressors (default: non-empty, so zero)
            
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    %% Create three-column regressors:

    % 1) Win ReqGo ActGo response:
    if sum(Win_ReqGo_ActGo_Trials) == 0
        tWin_ReqGo_ActGo = [0 0 0];
        emptyRegressorMatrix(1) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
        fprintf('Note: Subject %d has empty regressors tWin_ReqGo_ActGo! \n', iSub)
    else
        tWin_ReqGo_ActGo  = [tCueAll(Win_ReqGo_ActGo_Trials ==1) ones(sum(Win_ReqGo_ActGo_Trials),2)];
    end
    
    % 2) Win ReqGo ActNoGo response:
    if sum(Win_ReqGo_ActNoGo_Trials) == 0
        tWin_ReqGo_ActNoGo = [0 0 0];
        emptyRegressorMatrix(2) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
        fprintf('Note: Subject %d has empty regressors tWin_ReqGo_ActNoGo! \n', iSub)
    else
        tWin_ReqGo_ActNoGo = [tCueAll(Win_ReqGo_ActNoGo_Trials ==1) ones(sum(Win_ReqGo_ActNoGo_Trials),2)];
    end
    
    % 3) Win ReqGo ActGo response:
    if sum(Win_ReqNoGo_ActGo_Trials) == 0
        tWin_ReqNoGo_ActGo = [0 0 0];
        emptyRegressorMatrix(3) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
        fprintf('Note: Subject %d has empty regressors tWin_ReqNoGo_ActGo! \n', iSub)
    else
        tWin_ReqNoGo_ActGo  = [tCueAll(Win_ReqNoGo_ActGo_Trials ==1) ones(sum(Win_ReqNoGo_ActGo_Trials),2)];
    end
    
    % 4) Win ReqGo ActNoGo response:
    if sum(Win_ReqNoGo_ActNoGo_Trials) == 0
        tWin_ReqNoGo_ActNoGo = [0 0 0];
        emptyRegressorMatrix(4) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
        fprintf('Note: Subject %d has empty regressors tWin_ReqNoGo_ActNoGo! \n', iSub)
    else
        tWin_ReqNoGo_ActNoGo = [tCueAll(Win_ReqNoGo_ActNoGo_Trials ==1) ones(sum(Win_ReqNoGo_ActNoGo_Trials),2)];
    end
    
    % 5) Avoid ReqGo ActGo response:
    if sum(Avoid_ReqGo_ActGo_Trials) == 0
        tAvoid_ReqGo_ActGo = [0 0 0];
        emptyRegressorMatrix(5) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
        fprintf('Note: Subject %d has empty regressors tAvoid_ReqGo_ActGo! \n', iSub)
    else
        tAvoid_ReqGo_ActGo  = [tCueAll(Avoid_ReqGo_ActGo_Trials ==1) ones(sum(Avoid_ReqGo_ActGo_Trials),2)];
    end
    
    % 6) Avoid ReqGo ActNoGo response:
    if sum(Avoid_ReqGo_ActNoGo_Trials) == 0
        tAvoid_ReqGo_ActNoGo = [0 0 0];
        emptyRegressorMatrix(6) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
        fprintf('Note: Subject %d has empty regressors tAvoid_ReqGo_ActNoGo! \n', iSub)
    else
        tAvoid_ReqGo_ActNoGo = [tCueAll(Avoid_ReqGo_ActNoGo_Trials ==1) ones(sum(Avoid_ReqGo_ActNoGo_Trials),2)];
    end
    
    % 7) Avoid ReqGo ActGo response:
    if sum(Avoid_ReqNoGo_ActGo_Trials) == 0
        tAvoid_ReqNoGo_ActGo = [0 0 0];
        emptyRegressorMatrix(7) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
        fprintf('Note: Subject %d has empty regressors tAvoid_ReqNoGo_ActGo! \n', iSub)
    else
        tAvoid_ReqNoGo_ActGo  = [tCueAll(Avoid_ReqNoGo_ActGo_Trials ==1) ones(sum(Avoid_ReqNoGo_ActGo_Trials),2)];
    end
    
    % 8) Avoid ReqGo ActNoGo response:
    if sum(Avoid_ReqNoGo_ActNoGo_Trials) == 0
        tAvoid_ReqNoGo_ActNoGo = [0 0 0];
        emptyRegressorMatrix(8) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
        fprintf('Note: Subject %d has empty regressors tAvoid_ReqNoGo_ActNoGo! \n', iSub)
    else
        tAvoid_ReqNoGo_ActNoGo = [tCueAll(Avoid_ReqNoGo_ActNoGo_Trials ==1) ones(sum(Avoid_ReqNoGo_ActNoGo_Trials),2)];
    end    

    % Check for NAs:
    if sum(isnan(tWin_ReqGo_ActGo(:,3))) > 0; fprintf("Sub%d tWin_ReqGo_ActGo has NAs!",iSub); end
    if sum(isnan(tWin_ReqGo_ActNoGo(:,3))) > 0; fprintf("Sub%d tWin_ReqGo_ActNoGo has NAs!",iSub); end
    if sum(isnan(tWin_ReqNoGo_ActGo(:,3))) > 0; fprintf("Sub%d tWin_ReqNoGo_ActGo has NAs!",iSub); end
    if sum(isnan(tWin_ReqNoGo_ActNoGo(:,3))) > 0; fprintf("Sub%d tWin_ReqNoGo_ActNoGo has NAs!",iSub); end
    if sum(isnan(tAvoid_ReqGo_ActGo(:,3))) > 0; fprintf("Sub%d tAvoid_ReqGo_ActGo has NAs!",iSub); end
    if sum(isnan(tAvoid_ReqGo_ActNoGo(:,3))) > 0; fprintf("Sub%d tAvoid_ReqGo_ActNoGo has NAs!",iSub); end
    if sum(isnan(tAvoid_ReqNoGo_ActGo(:,3))) > 0; fprintf("Sub%d tAvoid_ReqNoGo_ActGo has NAs!",iSub); end
    if sum(isnan(tAvoid_ReqNoGo_ActNoGo(:,3))) > 0; fprintf("Sub%d tAvoid_ReqNoGo_ActNoGo has NAs!",iSub); end
    
    % Handedness:
    % 9) Motor response. Note that we take cue onset as timing here: Regressor 5
    tHand       = [tCueAll ones(length(tCueAll),1) motorNumAll]; % sum command checks how many incorrect responses

    % 10) Error trials. Note that we take cue onset as timing here.
    if sum(ErrorTrials==1) == 0 % if in fact no errors        
        tError = [0 0 0];
        fprintf('Note: Subject %d has empty regressor Error! \n', iSub)
        emptyRegressorMatrix(10) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
    else
        tError  = [tCueAll(ErrorTrials) ones(sum(ErrorTrials),2)]; % sum command checks how many incorrect responses
    end
     
    % 11-12) Outcome onset: neutral, reward after Go/NoGo, punishment after Go/NoGo.
    tOutcomeOnset       = [tFBAll(OutcomeTrials) ones(sum(OutcomeTrials),2)];
    tOutcomeValence     = [tFBAll(OutcomeTrials) ones(sum(OutcomeTrials),1) outcomeAll(OutcomeTrials)];

    % 13) Invalid trials. Uninstructed key pressed, thus no outcome, but error message.
    if sum(InvalidTrials) == 0 % if in fact no uninstructed key presses
        tInvalid = [0 0 0];
        fprintf('Note: Subject %d has empty regressor Invalid! \n', iSub)
        emptyRegressorMatrix(13) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
    else
        tInvalid = [tFBAll(InvalidTrials) ones(sum(InvalidTrials),2)]; % sum command checks how many incorrect responses
    end
    
    % 14) EEG conflict signal.
    tEEGConflict = [tCueAll ones(length(tCueAll),1) EEGRegressor1'];
    % 15) EEG conflict signal.
    % consider whether to use tCueAll (starts before response) or tRespAll (starts with response)
    tEEGGo = [tCueAll ones(length(tCueAll),1) EEGRegressor2'];

    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    %% Save regressor matrices:
    
    % Onsets:
    save(fullfile(dirs.timings,'tWin_ReqGo_ActGo.txt'),'tWin_ReqGo_ActGo','-ascii')
    save(fullfile(dirs.timings,'tWin_ReqGo_ActNoGo.txt'),'tWin_ReqGo_ActNoGo','-ascii')
    save(fullfile(dirs.timings,'tWin_ReqNoGo_ActGo.txt'),'tWin_ReqNoGo_ActGo','-ascii')
    save(fullfile(dirs.timings,'tWin_ReqNoGo_ActNoGo.txt'),'tWin_ReqNoGo_ActNoGo','-ascii')
    save(fullfile(dirs.timings,'tAvoid_ReqGo_ActGo.txt'),'tAvoid_ReqGo_ActGo','-ascii')
    save(fullfile(dirs.timings,'tAvoid_ReqGo_ActNoGo.txt'),'tAvoid_ReqGo_ActNoGo','-ascii')
    save(fullfile(dirs.timings,'tAvoid_ReqNoGo_ActGo.txt'),'tAvoid_ReqNoGo_ActGo','-ascii')
    save(fullfile(dirs.timings,'tAvoid_ReqNoGo_ActNoGo.txt'),'tAvoid_ReqNoGo_ActNoGo','-ascii')

    % Handedness, errors and outcomes:
    save(fullfile(dirs.timings,'tHand.txt'),'tHand','-ascii')
    save(fullfile(dirs.timings,'tError.txt'),'tError','-ascii')
    save(fullfile(dirs.timings,'tOutcomeOnset.txt'),'tOutcomeOnset','-ascii')
    save(fullfile(dirs.timings,'tOutcomeValence.txt'),'tOutcomeValence','-ascii')
    save(fullfile(dirs.timings,'tInvalid.txt'),'tInvalid','-ascii')
    
    % Save EEG regressors:
    save(fullfile(dirs.timings,'tEEGConflict.txt'),'tEEGConflict','-ascii')        
    save(fullfile(dirs.timings,'tEEGGo.txt'),'tEEGGo','-ascii')        
    
    % Save empty regressors:
    csvwrite(fullfile(dirs.timings,'emptyregressors.txt'),emptyRegressorMatrix)        
end % end iSub-loop.

% Check how many invalid trials per subject:
% [(1:36)' 640-sum(validTrialsAll,2)]

% Empty regressors:
% Starting with subject 1
% Note: Subject 1 has empty regressor Invalid! 
% Starting with subject 2
% Note: Subject 2 has empty regressor Invalid! 
% Starting with subject 3
% Starting with subject 4
% Note: Subject 4 has empty regressor Invalid! 
% Starting with subject 5
% Note: Subject 5 has empty regressor Invalid! 
% Starting with subject 6
% Note: Subject 6 has empty regressor Invalid! 
% Starting with subject 7
% Note: Subject 7 has empty regressor Invalid! 
% Starting with subject 8
% Starting with subject 9
% Note: Subject 9 has empty regressor Invalid! 
% Starting with subject 10
% Starting with subject 11
% Note: Subject 11 has empty regressors tWin_ReqGo_ActNoGo! 
% Starting with subject 12
% Note: Subject 12 has empty regressor Invalid! 
% Starting with subject 13
% Note: Subject 13 has empty regressor Invalid! 
% Starting with subject 14
% Note: Subject 14 has empty regressor Invalid! 
% Starting with subject 15
% Note: Subject 15 has empty regressor Invalid! 
% Starting with subject 16
% Starting with subject 17
% Note: Subject 17 has empty regressor Invalid! 
% Starting with subject 18
% Note: Subject 18 has empty regressor Invalid! 
% Starting with subject 19
% Note: Subject 19 has empty regressor Invalid! 
% Starting with subject 20
% Note: Subject 20 has empty regressor Invalid! 
% Starting with subject 21
% Note: Subject 21 has empty regressor Invalid! 
% Starting with subject 22
% Note: Subject 22 has empty regressor Invalid! 
% Starting with subject 23
% Note: Subject 23 has empty regressor Invalid! 
% Starting with subject 24
% Note: Subject 24 has empty regressors tWin_ReqGo_ActNoGo! 
% Note: Subject 24 has empty regressor Invalid! 
% Starting with subject 25
% Starting with subject 26
% Starting with subject 27
% Note: Subject 27 has empty regressor Invalid! 
% Starting with subject 28
% Note: Subject 28 has empty regressor Invalid! 
% Starting with subject 29
% Starting with subject 30
% Starting with subject 31
% Starting with subject 32
% Note: Subject 32 has empty regressor Invalid! 
% Starting with subject 33
% Note: Subject 33 has empty regressor Invalid! 
% Starting with subject 34
% Note: Subject 34 has empty regressor Invalid! 
% Starting with subject 35
% Starting with subject 36
% Note: Subject 36 has empty regressor Invalid! 

end % end of function.
