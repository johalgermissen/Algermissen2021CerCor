function taft_runonce_create_onsets()

% taft_runonce_create_onsets()
% 
% Loop over subject, load behavioral data files, extract timing of cue/
% response/ outcome for each trial for each subject, save under subject-/
% block-specific 'AROMA' directory in fMRI directory.
% Mind setting the root directory in taft_set_rootdir().
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Adapted from Tobias Hauser (https://github.com/tuhauser/TAfT).
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/CueLockedAnalyses/TAfT/

%% Set directories:

dirs.root   = taft_set_rootdir(); % /project/3017042.02
dirs.log    = fullfile(dirs.root,'Log'); 
dirs.behav  = fullfile(dirs.log,'Behavior/Data_beh_mat');    
dirs.fMRI   = fullfile(dirs.log,'fMRI');

%% Loop through subjects:

for iSub = 1:36 % iSub=1;

    fprintf('Load subject %03d\n',iSub);
    
    %% Load behavioral data file:
    
    behavfile = fullfile(dirs.behav,sprintf('3017042.02_emmvdij_%03d_001_results.mat',iSub));
    behavior = load(behavfile);

    %% Retrieve start of each block:
    
    % Retrieve absolute timing of start (trial 1) of each of first 3 blocks, repeat for each of the following 110 trials per block:
    t0part1 = reshape(repmat(behavior.tm.learn{1}.stim([1 111 221]),1,110)',[330 1]); % one long vector with 330 elements
    t0part1 = t0part1(1:320); % drop last 10 entries because block 3 was 10 trials shorter
    
    % Retrieve absolute timing of start (trial 1) of each of second 3 blocks, repeat for each of the 110 trials per block:
    t0part2 = reshape(repmat(behavior.tm.learn{2}.stim([1 111 221]),1,110)',[330 1]);
    t0part2 = t0part2(1:320); % drop last 10 entries because block 6 was 10 trials shorter
    t0      = [t0part1; t0part2] - 10; % concatenate; now subtract 10 seconds that were used for getting MRI signal steady-state (dummy scans)

    %% Retrieve timings of events corrected for t0 (onset of block):
    
    % Cue onset:
    tCueAll     = [behavior.tm.learn{1}.stim; behavior.tm.learn{2}.stim] - t0; % take cue onset; subtract start of each block

    % Response onset:
    tRespAll    = [behavior.tm.learn{1}.response; behavior.tm.learn{2}.response] - t0; % note that in case of nogo, the cue offset time was recorded.
    isGo        = [behavior.prep.seq.learn.resp{1}; behavior.prep.seq.learn.resp{2} ~= 0]; % whether Go or not
    tRespAll(isGo==0) = NaN; % delete timing for NoGos
    
    %  Feedback onset:
    tFBPrel        = [behavior.tm.learn{1}.outcome; behavior.tm.learn{2}.outcome]; % take onset of learning, so far uncorrected
    tWarningAll    = [behavior.tm.learn{1}.warning; behavior.tm.learn{2}.warning]; % take timing of warning messages, so far uncorrected
    tFBPrel(isnan(tFBPrel)) = tWarningAll; % replace NaNs with timing of warnings
    tFBAll         = tFBPrel - t0; % subtract start of each block --> now corrected
    if size(tFBAll,1) ~= 640; error('Feedback file has incorrect length'); end
%     [tCueAll tRespAll tFBAll] % check

    %% Create onset file for each block:
    
    trials2use = {1:110,111:220,221:320,321:430,431:540,541:640}; % trials per block
    
    for iBlock = 1:6 % for each block iBlock = 1;
        
        fprintf('Start subject %03d block %d\n',iSub,iBlock);
        
        targetDir = fullfile(dirs.fMRI,sprintf('sub-%03d/FEAT_Block%d.feat/AROMA',iSub,iBlock)); % target directory where to save:
        
        % Cue onset:
        out = tCueAll(trials2use{iBlock});
        outfile = fullfile(targetDir,'cueOnsets.txt');
        save(outfile,'out','-ascii');
        
        % Response onset:
        out = tRespAll(trials2use{iBlock});
        outfile = fullfile(targetDir,'responseOnsets.txt');
        save(outfile,'out','-ascii');
        
        % Outcome onset:
        out = tFBAll(trials2use{iBlock});
        outfile = fullfile(targetDir,'outcomeOnsets.txt');
        save(outfile,'out','-ascii');
        
    end % end iBlock
end % end iSub

end % end of function.
