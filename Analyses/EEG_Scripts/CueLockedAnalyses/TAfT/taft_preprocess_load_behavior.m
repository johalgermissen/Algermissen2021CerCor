function [out] = taft_preprocess_load_behavior(job)

% [out] = taft_preprocess_load_behavior(iSub)
%
% Load behavioral data, 
% recode response and accuracy variables (uninstructed keys to
% instructed keys; delete Gos with too long RTs, accuracy based on recoded
% keys presses),
% 
% INPUTS:
% job               = structure with settings for creating TAfT object, specifically:
% .behavFile        = string, full file path for behavioral raw data file of respective subject.
%
% OUTPUTS:
% out               = pre-processed behavioral data, with following fields:
% .stim             = numeric, identifiers of stimulus identify (1-16)
% .goCue            = numeric, identifiers of Go cues
% .reqaction        = numeric, index whether trial is Go cue (1) or not (0)
% .winCue           = numeric, identifiers of Win cues
% .valence          = numeric, index whether trial is Win cue (1) or not (2).
% .iswin            = numeric, index whether trial is Win cue (1) or not (0).
% .splithalf        = numeric, whether first (1) or second (2) session.
% .reqresp          = numeric, left Go (101) or right Go (97) or NoGo (0)
% required.
% .resp             = numeric, left Go (101) or right Go (97) or NoGo (0)
% performed (recoded: uninstructed keys to instructed keys, remove too long
% RTs).
% .reqgo            = numeric, Go (1) or NoGo (1) required.
% .isgo             = numeric, Go (1) or NoGo (1) performed (recoded: 
% uninstructed keys to instructed keys, remove too long RTs).
% .accuracy         = numeric, correct (1) or incorrect (0) based on
% reqresp and resp.
% .RT              	= numeric, RTs, NaN if too fast (< 0.2 sec.) or too
% slow (> 1.3035) or NoGo.
% .isW2G            = numeric, whether Go action to Win cue (1) or not (0).
% .isW2NG           = numeric, whether NoGo action to Win cue (1) or not (0).
% .isA2G            = numeric, whether Go action to Avoid cue (1) or not (0).
% .isA2NG           = numeric, whether NoGo action to Avoid cue (1) or not (0).
% 
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

%% Retrieve trialinfo from behavioural file.

load(job.behavFile) % load behavior

%% Combine data from both sessions (blocks 1-3 and 4-6) into one vector:

seq.stim        = [prep.seq.learn.stim{1}; prep.seq.learn.stim{2}];
seq.resp        = [results.learn{1}.response; results.learn{2}.response];
seq.RT          = [results.learn{1}.RT; results.learn{2}.RT];
seq.outcome     = [results.learn{1}.outcome; results.learn{2}.outcome];
seq.accuracy    = [results.learn{1}.acc; results.learn{2}.acc];
seq.go          = [results.learn{1}.go; results.learn{2}.go];
seq.splithalf   = [ones(160,1); 2*ones(160,1);ones(160,1); 2*ones(160,1)];

%% Carry forward stimulus properties: 

clear out
out.stim        = seq.stim;
out.goCue       = 1:4;
out.reqaction   = ~ismember(out.stim,out.goCue) + 1; % 1 if yes, 2 if no
out.winCue      = [1 2 5 6];
out.valence     = ~ismember(out.stim,out.winCue) + 1; % 1 if yes, 2 if no
out.iswin       = ismember(out.stim,out.winCue); % 1 is Win, 0 if Avoid

%% Initialize objects to recode: 

out.reqresp     = nan(size(out.stim)); % initialize in same format 
out.resp        = nan(size(out.stim)); % initialize in same format 
out.reqgo       = nan(size(out.stim)); % initialize in same format 
out.isgo        = nan(size(out.stim)); % initialize in same format 
out.isnextgo    = nan(size(out.stim)); % initialize in same format 
out.accuracy    = nan(size(out.stim)); % initialize in same format 
out.RT          = nan(size(out.stim)); % initialize in same format 

%% Loop over trials:

nTrials = length(seq.resp);

for iTrial = 1:nTrials % k = 1

    %% a) Recode recorded response (also for uninstructed key presses):
    
    if ismember(seq.resp(iTrial),[69,70,71,72,101,102,103,104]) % left Go response
            out.resp(iTrial) = 101; % actually instructed left Go response is 101
            out.isgo(iTrial) = 1; % any Go
            
    elseif ismember(seq.resp(iTrial),[65,66,67,68,97,98,99,100]) % right Go response
            out.resp(iTrial) = 97; % actually instructed right Go response is 97
            out.isgo(iTrial) = 1; % any Go
   
    else % NoGo response
            out.resp(iTrial) = 0; % no response = NoGo; implicitly also sets NAs to NoGo
            out.isgo(iTrial) = 0; % any NoGo
    end
    
    %% b) Recode required response (based on stimulus number):
    
    if ismember(out.stim(iTrial),[1, 3]) % left Go response
            out.reqresp(iTrial) = 101; % actually instructed left Go response is 101
            out.reqgo(iTrial)   = 1;
            
    elseif ismember(out.stim(iTrial),[2, 4]) % right Go response
            out.reqresp(iTrial) = 97; % actually instructed right Go response is 97
            out.reqgo(iTrial)   = 1;
            
    else % NoGo response
            out.reqresp(iTrial) = 0; % no response = NoGo; implicitly also sets NAs to NoGo
            out.reqgo(iTrial)   = 0;
    end
    
    %% c) Recode RTs:
    
    if out.isgo(iTrial) == 0
        out.RT(iTrial)      = NaN; % NoGo response

    elseif seq.RT(iTrial) < 0.2 % if too short: delete
        out.RT(iTrial)      = NaN;
        
    elseif seq.RT(iTrial) > 1.3035 % if too long: delete, set to NoGo
        out.RT(iTrial)      = NaN;
        out.resp(iTrial)    = 0; % late response = NoGo
        out.isgo(iTrial)    = 0; % 
        
    else
        out.RT(iTrial)      = round(seq.RT(iTrial),3); % round to
        % ms be consistent with EEG sampling rate (in ms)
    end
    
    %% d) Recode accuracy based on recorded responses:
    
    if out.reqresp(iTrial) ~= out.resp(iTrial) || seq.RT(iTrial) > 1.3035 % wrong response (key) or after time
        out.accuracy(iTrial) = 0; % 0 if incorrect
        
    else
        out.accuracy(iTrial) = 1; % 1 if correct
    end
    
end % end of iTrial loop.

%% Executed action x valence:

out.isW2G   = out.iswin == 1 & out.isgo == 1;
out.isW2NG  = out.iswin == 1 & out.isgo == 0;
out.isA2G   = out.iswin == 0 & out.isgo == 1;
out.isA2NG  = out.iswin == 0 & out.isgo == 0;

end % end of function.