function out = load_recode_behav(rootdir,iSub,rejectedtrials)

% out = load_recode_behav(rootdir,iSub,rejectedtrials) 
%
% Loads behavioral data and filters out rejected trials.
%
% INPUT:
% rootdir           = string, root directory of project (needed to find
% behavioral data)
% iSub              = scalar integer, number of subject which to retrieve
% data from
% rejectedtrials    = vector of length 1 x number of trials containing
% whether trials rejected (1) or not (0)
%
% OUTPUT:
% out               = pre-processed behavioral data
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

%% Indices of valid trials.

trlidx = find(~rejectedtrials);

%% Load behavioral data:

fprintf('Subject %03d: Load behavioral data\n',iSub);
load(fullfile(rootdir,sprintf('Log/Behavior/Data_beh_mat/3017042.02_emmvdij_%03d_001_results.mat',iSub)))

%% Combine both sessions:

fprintf('Subject %03d: Combine both sessions of behavioral data\n',iSub);

seq.stim        = [prep.seq.learn.stim{1}; prep.seq.learn.stim{2}];
seq.resp        = [results.learn{1}.response; results.learn{2}.response];
seq.RT          = [results.learn{1}.RT; results.learn{2}.RT];
seq.outcome     = [results.learn{1}.outcome; results.learn{2}.outcome];
seq.accuracy    = [results.learn{1}.acc; results.learn{2}.acc];
seq.go          = [results.learn{1}.go; results.learn{2}.go];
seq.splithalf   = [ones(160,1); 2*ones(160,1);ones(160,1); 2*ones(160,1)];

%% Select data only for non-rejected trials:

fprintf('Subject %03d: Remove rejected trials from behavioral data\n',iSub);

out.trlidx          = trlidx;

% Stimuli:
out.stim            = seq.stim(trlidx);
out.goCue           = 1:4;
out.winCue          = [1 2 5 6];
out.action          = ~ismember(out.stim,out.goCue) + 1; % 1 if yes, 2 if no
out.val             = ~ismember(out.stim,out.winCue) + 1; % 1 if yes, 2 if no

% Response:
out.resp            = seq.resp(trlidx);
out.RT              = seq.RT(trlidx);
out.accuracy        = seq.accuracy(trlidx);
out.go              = seq.go(trlidx);

% Outcome:
out.outcome         = seq.outcome(trlidx); % what subject see: -1, 0, 1
out.fb              = out.outcome + 1; % 1 becomes bad, 2 becomes good
out.fb(out.val==2)  = out.fb(out.val==2)+1; % Avoid cues: one up
out.splithalf       = seq.splithalf(trlidx);

end % end function