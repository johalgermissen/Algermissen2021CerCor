function EEGfMRIPav_GLM1()

% EEGfMRIPav_GLM1
% 
% Create regressors for subject-level GLM (concatenated blocks).
% Allows for contrasts valence, performed action, conflict (4 regressors).
% Use three-column (onset, duration, factor level) format.
% Saves each regressor as separate .txt file in respective
% timings_regressor folder of respective GLM of respective subject.
%
% Regressors of interest are:
% Valence (Win/Avoid) x Actual Action (Go/NoGo)
% (Go/NoGo):
% 1) Win_ActGo
% 2) Win_ActNoGo
% 3) Avoid_ActGo
% 4) Avoid_ActNoGo
%
% Covariates of no interest are:
% 5) Handedness of response: hand (left, NoGo, right)
% 6) Performance: Error (correct/incorrect)
% 7) Outcome Onset
% 8) Outcome Valence
% 9) Invalid Outcome
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

nSub        = 36;

% GLM-specific settings:
GLMID       = '1';
nRegressors = 9; % for specification of empty regressors

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
    RTCorAll       = repelem(NaN,length(RTAll))'; % RTs

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
        stim                    = stimAll(iTrial);
        % Store whether cue valence experienced in this trial or not
        cueDetectAll(iTrial)    = cueDetected(stim);        
        % Retrieve actually received outcome:
        r                       = outcomeAll(iTrial);
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
    %% Create timing for regressors:

    % Strategy for concatenated blocks:
    % a) Retrieve number of volumes per block, compute effective fMRITime in GLM as TE * position of volume
    % b) Extract when in fMRITime each blocks starts,
    % c) extract for when block started in behavioral time (set to 0),
    % d) extract timing of trials from behavioral data files and correct for when blocks started in behavioral time,
    % e) add when block started in fMRITime    
    
    % a) retrieve block lengths in volumes, create fMRI time:
    dirs.sub            = fullfile(dirs.root,sprintf('Log/fMRI/sub-%03d',iSub));
    nVolumes            = load(fullfile(dirs.sub,sprintf('Sub%03d_Nvolumes.txt',iSub))); % number of volumes for all 6 blocks
    TE                  = 1.4; % echo time
    fMRITime            = [0:(sum(nVolumes)-1)]*TE; % when each volume starts
    % How many volumes have occurred till end of each block, minus the block itself (so rather at the start of each block), +1 to start at index 1:
    fMRIBlockStart      = cumsum(nVolumes)-nVolumes+1; 
    
    % b) extract when each block started in fMRI time:
    trialBlockLength    = [110 110 100 110 110 100];
    trialTime0          = nan(0,0); % time of first volume per block
    for i = 1:length(trialBlockLength)
        trialTime0   = [trialTime0; repmat(fMRITime(fMRIBlockStart(i)),1,trialBlockLength(i))']; % repeat for each volume in block   
    end    
       
    % c) extract when blocks started in behavioral time:
    % extract timing at start of blocks, repeat 110 (or 100) times, 
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
    tFBPrel             = [tm.learn{1}.outcome; tm.learn{2}.outcome]; % take onset of learning, so far uncorrected
    tWarningAll         = [tm.learn{1}.warning; tm.learn{2}.warning]; % take timing of warning messages, so far uncorrected
    tFBPrel(isnan(tFBPrel)) = tWarningAll; % replace NaNs with timing of warnings
    tFBAllBlock         = tFBPrel - t0; % subtract start of each block --> now corrected
    % sum(isnan(tFBAllBlock)) % check whether NaN due to uninstructed key presses (NAs left)?
    
    % e) add start of blocks in fMRI time:
    tCueAll             = tCueAllBlock+trialTime0;
    tRespAll            = tRespAllBlock+trialTime0;
    tFBAll              = tFBAllBlock+trialTime0;

    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    %% Create behavioral aspects of regressors:

    % 1) Cues:
    % Split up trials into those with Win vs. Avoid cues and Go vs. NoGo responses:
    % mind cueDetectAll: don't model cues for trials where outcome not yet experienced
    Win_ActGo_Trials         = (ismember(stimIDAll,[1 2 5 6]) & ismember(respNumAll,[1 2]) & cueDetectAll==1); % indicator whether trial was Win trial with Go response
    Win_ActNoGo_Trials       = (ismember(stimIDAll,[1 2 5 6]) & ismember(respNumAll,[3]) & cueDetectAll==1); % indicator whether trial was Win trial with Go response
    Avoid_ActGo_Trials       = (ismember(stimIDAll,[3 4 7 8]) & ismember(respNumAll,[1 2]) & cueDetectAll==1); % indicator whether trial was Win trial with Go response
    Avoid_ActNoGo_Trials     = (ismember(stimIDAll,[3 4 7 8]) & ismember(respNumAll,[3]) & cueDetectAll==1); % indicator whether trial was Win trial with Go response

    % Handedness: motorNumAll created above
    
    % 2) Responses:
    ErrorTrials = incorrectAll==1;
    
    % 3) Feedback:
    OutcomeTrials = validOutcome==1; 
    InvalidTrials   = validOutcome==0;  
    
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    %% Create three-column regressors:

    emptyRegressorMatrix = zeros(1,nRegressors); % initialize empty regressors (default: non-empty, so zero)
            
    % 1) Win ActGo response:
    if sum(Win_ActGo_Trials) == 0
        tWin_ActGo = [0 0 0];
        emptyRegressorMatrix(1) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
        fprintf('Note: Subject %d has empty regressors tWin_ActGo! \n', iSub)
    else
        tWin_ActGo  = [tCueAll(Win_ActGo_Trials ==1) ones(sum(Win_ActGo_Trials),2)];
    end
    
    % 2) Win ActNoGo response:
    if sum(Win_ActNoGo_Trials) == 0
        tWin_ActNoGo = [0 0 0];
        emptyRegressorMatrix(2) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
        fprintf('Note: Subject %d has empty regressors tWin_ActNoGo! \n', iSub)
    else
        tWin_ActNoGo = [tCueAll(Win_ActNoGo_Trials ==1) ones(sum(Win_ActNoGo_Trials),2)];
    end
    
    % 3) Avoid ActGo response:
    if sum(Avoid_ActGo_Trials) == 0
        tAvoid_ActGo = [0 0 0];
        emptyRegressorMatrix(3) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
        fprintf('Note: Subject %d has empty regressors tAvoid_ActGo! \n', iSub)
    else
        tAvoid_ActGo  = [tCueAll(Avoid_ActGo_Trials ==1) ones(sum(Avoid_ActGo_Trials),2)];
    end
    
    % 4) Avoid ActNoGo response:
    if sum(Avoid_ActNoGo_Trials) == 0
        tAvoid_ActNoGo = [0 0 0];
        emptyRegressorMatrix(4) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
        fprintf('Note: Subject %d has empty regressors tAvoid_ActNoGo! \n', iSub)
    else
        tAvoid_ActNoGo = [tCueAll(Avoid_ActNoGo_Trials ==1) ones(sum(Avoid_ActNoGo_Trials),2)];
    end

    % Check for NAs:
    if sum(isnan(tWin_ActGo(:,3))) > 0; fprintf("Sub%d tWin_ActGo has NAs!",iSub); end
    if sum(isnan(tWin_ActNoGo(:,3))) > 0; fprintf("Sub%d tWin_ActNoGo has NAs!",iSub); end
    if sum(isnan(tAvoid_ActGo(:,3))) > 0; fprintf("Sub%d tAvoid_ActGo has NAs!",iSub); end
    if sum(isnan(tAvoid_ActNoGo(:,3))) > 0; fprintf("Sub%d tAvoid_ActNoGo has NAs!",iSub); end
    
    % Handedness:
    % 5) Motor response. Note that we take cue onset as timing here: Regressor 5
    tHand       = [tCueAll ones(length(tCueAll),1) motorNumAll]; % sum command checks how many incorrect responses

    % 6) Error trials. Note that we take cue onset as timing here.
    if sum(ErrorTrials==1) == 0 % if in fact no errors        
        tError = [0 0 0];
        fprintf('Note: Subject %d has empty regressor Error! \n', iSub)
        emptyRegressorMatrix(6) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
    else
        tError  = [tCueAll(ErrorTrials) ones(sum(ErrorTrials),2)]; % sum command checks how many incorrect responses
    end
        
    % 7-8) Outcome onset: neutral, reward after Go/NoGo, punishment after Go/NoGo.
    tOutcomeOnset       = [tFBAll(OutcomeTrials) ones(sum(OutcomeTrials),2)];
    tOutcomeValence     = [tFBAll(OutcomeTrials) ones(sum(OutcomeTrials),1) outcomeAll(OutcomeTrials)];

    % 9) Invalid trials. Uninstructed key pressed, thus no outcome, but error message.
    if sum(InvalidTrials) == 0 % if in fact no uninstructed key presses
        tInvalid = [0 0 0];
        fprintf('Note: Subject %d has empty regressor Invalid! \n', iSub)
        emptyRegressorMatrix(9) = 1; % specify that regressor is empty (regressor position hard-coded given that GLM setup is known)
    else
        tInvalid = [tFBAll(InvalidTrials) ones(sum(InvalidTrials),2)]; % sum command checks how many incorrect responses
    end
    
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    % ------------------------------------------------------------------- %
    %% Save regressor matrices:
    
    % Cues:
    save(fullfile(dirs.timings,'tWin_ActGo.txt'),'tWin_ActGo','-ascii')
    save(fullfile(dirs.timings,'tWin_ActNoGo.txt'),'tWin_ActNoGo','-ascii')
    save(fullfile(dirs.timings,'tAvoid_ActGo.txt'),'tAvoid_ActGo','-ascii')
    save(fullfile(dirs.timings,'tAvoid_ActNoGo.txt'),'tAvoid_ActNoGo','-ascii')
    
    % Handedness, errors, and outcomes:
    save(fullfile(dirs.timings,'tHand.txt'),'tHand','-ascii')
    save(fullfile(dirs.timings,'tError.txt'),'tError','-ascii')
    save(fullfile(dirs.timings,'tOutcomeOnset.txt'),'tOutcomeOnset','-ascii')
    save(fullfile(dirs.timings,'tOutcomeValence.txt'),'tOutcomeValence','-ascii')
    save(fullfile(dirs.timings,'tInvalid.txt'),'tInvalid','-ascii')
    
    % Save empty regressors:
    csvwrite(fullfile(dirs.timings,'emptyregressors.txt'),emptyRegressorMatrix)  
    
end % end iSub-loop.

% Check how many invalid trials per subject:
% [(1:36)' 640-sum(validTrialsAll,2)]

% Empty regressors:
% Subject 001: Start 
% Note: Subject 1 has empty regressor Invalid! 
% Subject 002: Start 
% Note: Subject 2 has empty regressor Invalid! 
% Subject 003: Start 
% Subject 004: Start 
% Note: Subject 4 has empty regressor Invalid! 
% Subject 005: Start 
% Note: Subject 5 has empty regressor Invalid! 
% Subject 006: Start 
% Note: Subject 6 has empty regressor Invalid! 
% Subject 007: Start 
% Note: Subject 7 has empty regressor Invalid! 
% Subject 008: Start 
% Subject 009: Start 
% Note: Subject 9 has empty regressor Invalid! 
% Subject 010: Start 
% Subject 011: Start 
% Subject 012: Start 
% Note: Subject 12 has empty regressor Invalid! 
% Subject 013: Start 
% Note: Subject 13 has empty regressor Invalid! 
% Subject 014: Start 
% Note: Subject 14 has empty regressor Invalid! 
% Subject 015: Start 
% Note: Subject 15 has empty regressor Invalid! 
% Subject 016: Start 
% Subject 017: Start 
% Note: Subject 17 has empty regressor Invalid! 
% Subject 018: Start 
% Note: Subject 18 has empty regressor Invalid! 
% Subject 019: Start 
% Note: Subject 19 has empty regressor Invalid! 
% Subject 020: Start 
% Note: Subject 20 has empty regressor Invalid! 
% Subject 021: Start 
% Note: Subject 21 has empty regressor Invalid! 
% Subject 022: Start 
% Note: Subject 22 has empty regressor Invalid! 
% Subject 023: Start 
% Note: Subject 23 has empty regressor Invalid! 
% Subject 024: Start 
% Note: Subject 24 has empty regressor Invalid! 
% Subject 025: Start 
% Subject 026: Start 
% Subject 027: Start 
% Note: Subject 27 has empty regressor Invalid! 
% Subject 028: Start 
% Note: Subject 28 has empty regressor Invalid! 
% Subject 029: Start 
% Subject 030: Start 
% Subject 031: Start 
% Subject 032: Start 
% Note: Subject 32 has empty regressor Invalid! 
% Subject 033: Start 
% Note: Subject 33 has empty regressor Invalid! 
% Subject 034: Start 
% Note: Subject 34 has empty regressor Invalid! 
% Subject 035: Start 
% Subject 036: Start 
% Note: Subject 36 has empty regressor Invalid!  

end % end of function
