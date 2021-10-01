function [out] = EEGfMRIPav_aggr_data()

% [out] = EEGfMRIPav_aggr_data()
% 
% Loads behavioral data and 
% - recodes if necessary (uninstructed key presses, too short/long RTs,
% updated accuracy)
% - recodes outcomes
% - sorts data per stimulus
% - sorts data per condition
% - creates Win-stay/lose-shift variable
% - save as csv; return as function argument under /project/3017042.02/Log/Behavior/Matlab_behav_recode.csv
%
% INPUT:
% none
%
% Output:
% out           = cell with relevant behavioral data objects.
% Additionally saved as .csv argument under 
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

%% Set directories;

dirs.root = '/project/3017042.02';

%% General settings;

nSub                = 36;
nStim               = 16;
nCond               = 4;
nTrial              = 640;
nRep                = nTrial/nStim;

%% Initialize output fields as NaN:

out.stim            = nan(nSub,nTrial);
out.isgo            = nan(nSub,nTrial);
out.outcome         = nan(nSub,nTrial);
out.pGoStim         = nan(nSub,nStim,nRep);
out.pCorrectStim    = nan(nSub,nStim,nRep);
out.RTStim          = nan(nSub,nStim,nRep);
out.pGoCond         = nan(nSub,nCond,nRep);
out.pCorrectCond    = nan(nSub,nCond,nRep);
out.RTCond          = nan(nSub,nCond,nRep);
out.RTAcc           = nan(nSub,6);
saveData            = nan(nSub*nTrial,9); % subject, trialnr, valence, required action, action, response, outcome, stay

%% Loop over subjects:

for iSub = 1:nSub % iSub = 1;
    
    %% Load data:
    
    fprintf('Start loading subject %02d\n',iSub);
    load(fullfile(dirs.root,'Log/Behavior/Data_beh_mat',sprintf('3017042.02_emmvdij_%03d_001_results.mat',iSub)));
    fprintf('Finished loading subject %02d\n',iSub);
    
    %% Concatenate both sessions (blocks 1-3 and 4-6):

    fprintf('Concatenate both sessions\n')
    
    seq.stim        = [prep.seq.learn.stim{1}; prep.seq.learn.stim{2}];
    seq.stimID      = [prep.seq.learn.stim{1}; prep.seq.learn.stim{2}+8];
    seq.resp        = [results.learn{1}.response; results.learn{2}.response];
    seq.RT          = [results.learn{1}.RT; results.learn{2}.RT];
    seq.outcome     = [results.learn{1}.outcome; results.learn{2}.outcome];
    seq.accuracy    = [results.learn{1}.acc; results.learn{2}.acc];
    seq.go          = [results.learn{1}.go; results.learn{2}.go];
    seq.splithalf   = [ones(160,1); 2*ones(160,1);ones(160,1); 2*ones(160,1)];

    %% A) Carry forward conditions and outcomes:

    fprintf('Carry forward cue conditions and outcomes \n')

    clear dat
    dat.stim        = seq.stim;
    dat.stimID      = seq.stimID;
    dat.reqaction   =  ~ismember(dat.stim,1:4)+1; % 1 if Go, 2 if NoGo
    dat.valence     = ~ismember(dat.stim,[1 2 5 6])+1; % 1 if Win, 2 if Avoid
    dat.iswin       = ismember(dat.stim,[1 2 5 6]); % Win cue or not
    dat.outcome     = seq.outcome; % outcome
    dat.outVal      = seq.outcome + (dat.valence==2); % for Avoid: add 1, to -1 becomes 0, 0 becomes +1
    % [dat.valence dat.outcome dat.outVal]

    %% RECODE RESPONSES AND ACCURACY:
    
    fprintf('Recode responses and accuracy \n')
    
    % Initialize data objects:
    dat.reqresp     = nan(nTrial,1); % initialize in same format 
    dat.resp        = nan(nTrial,1); % initialize in same format 
    dat.reqgo       = nan(nTrial,1); % initialize in same format 
    dat.isgo        = nan(nTrial,1); % initialize in same format 
    dat.isnextgo    = nan(nTrial,1); % initialize in same format 
    dat.accuracy    = nan(nTrial,1); % initialize in same format 
    dat.RT          = nan(nTrial,1); % initialize in same format 
    dat.stimRep     = nan(nTrial,1); % initialize in same format 
    stimCount       = zeros(nStim,1); % count how often stimulus has already appeared

    % Loop over trials:
    for iTrial = 1:nTrial % iTrial = 1
        
        %% a) Count how often stimulus has appeared so far
        
        thisStim            = dat.stimID(iTrial); % retrieve stimulus ID
        stimCount(thisStim) = stimCount(thisStim) + 1; % increment count for this stimulus ID
        dat.stimRep(iTrial) = stimCount(thisStim); % store stimulus count
        
        %% b) Recode recorded response (also for uninstructed key presses):
        
        if ismember(seq.resp(iTrial),[69,70,71,72,101,102,103,104]) % left Go response
                dat.resp(iTrial) = 101; % actually instructed left Go response is 101
                dat.isgo(iTrial) = 1; % any Go
                
        elseif ismember(seq.resp(iTrial),[65,66,67,68,97,98,99,100]) % right Go response
                dat.resp(iTrial) = 97; % actually instructed right Go response is 97
                dat.isgo(iTrial) = 1; % any Go
                
        else % NoGo response
                dat.resp(iTrial) = 0; % no response = NoGo; implicitly also sets NAs to NoGo
                dat.isgo(iTrial) = 0; % any NoGo
        end
        
        %% c) Reconstruct required response (based on stimulus number):
        
        if ismember(dat.stim(iTrial),[1, 3]) % left Go response
                dat.reqresp(iTrial) = 101; % actually instructed left Go response is 101
                dat.reqgo(iTrial)   = 1;
                
        elseif ismember(dat.stim(iTrial),[2, 4]) % right Go response
                dat.reqresp(iTrial) = 97; % actually instructed right Go response is 97
                dat.reqgo(iTrial)   = 1;
                
        else % NoGo response
                dat.reqresp(iTrial) = 0; % no response = NoGo; implicitly also sets NAs to NoGo
                dat.reqgo(iTrial)   = 0;
        end
        
        %% d) Recode RTs:
        
        if seq.RT(iTrial) < 0.2 % if too short: delete
            dat.RT(iTrial)      = NaN;
            
        elseif seq.RT(iTrial) > 1.3035 % if too long: delete, set to NoGo
            
            dat.RT(iTrial)      = NaN;
            dat.resp(iTrial)    = 0; % late response = NoGo
            dat.isgo(iTrial)    = 0; % late response = NoGo 
        else
            
            dat.RT(iTrial) = seq.RT(iTrial); % keep RT
        end
        
        %% e) Recode accuracy based on recorded responses:
        
        if dat.reqresp(iTrial) ~= dat.resp(iTrial) || seq.RT(iTrial) > 1.3035 % if mismatch required and performed action OR too late response: wrong
            dat.accuracy(iTrial) = 0; % 0 if incorrect
            
        else
            dat.accuracy(iTrial) = 1; % 1 if correct
        end
        
    end % end of iTrial
    
    %% Recode outcome types:
    
    fprintf('Recode outcome\n')
    
    dat.fb.abs                  = dat.outcome; % absolute outcome: +1, 0, -1
    dat.fb.rel                  = 1 - dat.outcome; % Avoid cues: 0-->1 becomes good, -1-->2 becomes bad; mind next line for completion
    dat.fb.rel(dat.valence==1)  = dat.fb.rel(dat.valence==1)+1; % Win cues: one up : one up (1-->0-->1 becomes 1, 0-->1-->2 becomes 2)
    dat.fb.all                  = dat.fb.rel+2*(dat.valence==2); % 1-4: reward, no-reward, no-punishment, punishment
%     [dat.valence dat.fb.abs dat.fb.rel dat.fb.all] % check valence, abs, rel, all

    %% Determine win stay/ lose stay behavior (after previous recoding): 
    
    fprintf('Determine stay vs. shift on next cue occurrence \n')
    
    dat.stay        = nan(nTrial,1); % initialize in same format 
    
    for iTrial = 1:nTrial % iTrial = 1;
        
        thisStim    = dat.stimID(iTrial); % retrieve stimulus identifier
        thisRep     = dat.stimRep(iTrial); % retrieve cumulative number of repetitions of this stimulus
        nextTrial   = find(dat.stimID == thisStim & dat.stimRep == (thisRep+1)); % next trial with same stimulus
        
        if thisRep < nRep % if there is still a "next trial" left
            dat.stay(iTrial) = double(dat.resp(iTrial) == dat.resp(nextTrial)); % decide whether exact response (resp) repeated or not
        end
    end
    % End data pre-processing
    
    % ---------------------------------------------------------------------
    %% Store data per subject:
    
    fprintf('Save subject data to big matrix to be stored as .csv \n')

    idx             = ((iSub-1)*nTrial+1):(iSub*nTrial); % trial indices for this subject
    saveData(idx,1) = iSub; % subject
    saveData(idx,2) = 1:nTrial; % trialnr
    saveData(idx,3) = dat.stim; % stimulus ID
    saveData(idx,4) = dat.valence; % valence
    saveData(idx,5) = dat.reqaction; % required action
    saveData(idx,6) = dat.isgo; % action
    saveData(idx,7) = dat.resp; % response
    saveData(idx,8) = dat.outcome; % outcome
    saveData(idx,9) = dat.stay; % stay
    % End trial-by-trial-data, start aggregation
    
    % ---------------------------------------------------------------------
    %% Save stimulus order per subject
    
    fprintf('Save stimulus order in out object \n')
    
    out.stim(iSub,:)                    = seq.stimID;
    out.isgo(iSub,:)                    = dat.isgo;
    out.outcome(iSub,:)                 = dat.outcome;
    
    %% Sort choices/ accuracy/ RTs per CUE:
    
    fprintf('Sort data per cue\n');
    % pGoCond of size nSub x nCond x nRep
    % pCorrectCond of size nSub x nCond x nRep

    for iStim = 1:nStim % iStim = 1;
        
        idx = find(dat.stimID == iStim); % indices of trials with this stimulus:
        if length(idx) ~= nRep; warning('Subject %02d iStim = %d: idx is %d trials',iSub,iStim,length(idx)); end
        
        out.pGoStim(iSub,iStim,:)       = dat.isgo(idx);
        out.pCorrectStim(iSub,iStim,:)  = dat.accuracy(idx);
        out.stayStim(iSub,iStim,:)      = dat.stay(idx);
        out.RTStim(iSub,iStim,:)        = dat.RT(idx); % average over time
        
    end
    
    %% Aggregate choices/ accuracy/ RT per CONDITION:
    
    fprintf('Aggregate choices/ accuracy per condition\n');
    % pGoCond of size nSub x nCond x nRep
    % pCorrectCond of size nSub x nCond x nRep
    
    for iCond = 1:nCond % iCond = 1;

        idx                             = [iCond*2-1 iCond*2 iCond*2-1+8 iCond*2+8]; % indices of cues in this condition
        
        out.pGoCond(iSub,iCond,:)       = mean(out.pGoStim(iSub,idx,:),2);
        out.pCorrectCond(iSub,iCond,:)  = mean(out.pCorrectStim(iSub,idx,:),2);
        out.stayCond(iSub,iCond,:)      = nanmean(out.stayStim(iSub,idx,:),2);
        out.RTCond(iSub,iCond,:)        = nanmean(out.RTStim(iSub,idx,:),2);
        
    end % end iCond
    
    %% RTs: correct Go, incorrect Go (other Go), incorrect Go (NoGo)
    
    fprintf('Sort RTs per correct Go/ incorrect Go (other Go)/ incorrect Go (NoGo) \n')
    
    allReqGo    = [1 1 1 1 0 0];
    allAcc      = [1 1 0 0 0 0];
    allVal      = [1 2 1 2 1 2];
    
    for iCond = 1:6 % iCond = 1;
        
        idx = dat.reqgo==allReqGo(iCond) & dat.accuracy==allAcc(iCond) & dat.valence==allVal(iCond); % trials in this condition
        out.RTAcc(iSub,iCond) = nanmean(dat.RT(idx)); % mean RT in this condition
        
    end
    
    %% Win-stay/ lose stay (without Go):
    
    fprintf('Compute Win-stay/ lose-shift \n')
    
    allVals = [1 0];
    for iVal = 1:length(allVals)
        valIdx                      = allVals(iVal);
%         idx                         = find(dat.outVal==valIdx); % positive/ negative
        idx                         = dat.outVal==valIdx; % positive/ negative
        out.stayOutVal(iSub,iVal)   = nanmean(dat.stay(idx));
    end
    
    %% Win-stay/ lose stay separate for Go/NoGo:
    
    allGos  = [1 0];
    allVals = [1 2 3 4];
    iCond = 0;
    for iGo = 1:length(allGos)
        for iVal = 1:length(allVals)
            iCond = iCond + 1;
            goIdx                           = allGos(iGo);
            valIdx                          = allVals(iVal);
            idx                             = dat.fb.all==valIdx & dat.isgo==goIdx;
            out.stayGoOutVal(iSub,iCond)    = nanmean(dat.stay(idx));
        end
    end
    
end % end iSub

%% Save as csv, but also return as function argument:

fprintf('Save file of all subjects as .csv\n')
csvwrite(fullfile(dirs.root,'Log/Behavior/Matlab_behav_recode.csv'),saveData);

fprintf('Finished loading all subjects\n');

end % end of function