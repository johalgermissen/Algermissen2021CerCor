function taft_realignment_per_condition()

% Load summary measure of realignment parameters (run taft_relignment_sum
% first), upsample and epoch into trials, split per behavioral conditions, 
% t-test for difference between relevant conditions.
% Automatically loads realignment parameter and behavior.

% Use as 
%   taft_realignment_behav()

% By Johannes Algermissen, 2019, adapted from Tobias Hauser (https://github.com/tuhauser/TAfT)
% works in Matlab 2018b.

%% Delete open objects/ screens/ clean console

clc
clear all
close all

%% Root directory:

dirs.root   = taft_set_rootdir(); % /project/3017042.02

%% Add necessary functions:

addpath(fullfile(rootdir,'EEG_Scripts/CueLockedAnalyses/TAfT'))

%% Directories:

dirs.log    = fullfile(dirs.root,'Log');
dirs.behav  = fullfile(dirs.log,'Behavior/Data_beh_mat');    
dirs.fMRI   = fullfile(dirs.log,'fMRI');    

%% Settings:

nSub        = 36;
nBlock      = 6;
subMotion   = nan(nSub,640);

%% Loop through subjects:

for iSub = 1:nSub % loop through subjects
    
    blockMotion = []; % initialize block motion
    % unclear how many volumes, thus hard to pre-allocate object size

    for iBlock = 1:nBlock
        
        % Upsample and epoch data:
        ROIdef.rawfMRIfile          = {fullfile(dirs.fMRI,sprintf('sub-%03d/FEAT_Block%d.feat/AROMA/mp_Fellner.txt',iSub,iBlock))}; % summary measure of motion
        ROIdef.onsets               = {fullfile(dirs.fMRI,sprintf('sub-%03d/FEAT_Block%d.feat/AROMA/trialOnsets.txt',iSub,iBlock))}; % trial onsets
        [data,~,~,~]                = taft_filter_upsample_epoch(ROIdef,128,10,1,1.4,2); % upsample and epoch in into trials
        data                        = nanmean(data,2); % average over time within trials (just take mean per trial)
        blockMotion                 = [blockMotion data']; % concatenate

    end
    
    subMotion(iSub,:)               = blockMotion; % store in matrix with row for each subject
    
    % Load behavior and store in field:
    
    behav(iSub)                     = taft_load_behavior(iSub); % store in field

end

%% Select valid subjects:

% Select valid subjects:
% invalidSubs = []; % include all subjects
% invalidSubs = [15 25]; % bad co-registrations
invalidSubs = [1 11 15 19 21 25 26]; % invalid co-registrations and outliers
validSubs   = setdiff(1:nSub,invalidSubs);
nSubValid   = length(validSubs);
fprintf('Extract valid subjects, exclude subjects %s\n',num2str(invalidSubs));

%% Mean motion for behavioral conditions:

for iSub = 1:nSubValid %loop over subjects:

    subID               = validSubs(iSub); % respective ID of valid subject
    motionGo(iSub)      = mean(subMotion(subID,behav(subID).isgo'==1));
    motionNoGo(iSub)    = mean(subMotion(subID,behav(subID).isgo'==0));
    motionWin(iSub)     = mean(subMotion(subID,behav(subID).iswin'==1));
    motionAvoid(iSub)   = mean(subMotion(subID,behav(subID).iswin'==0));
    motionCong(iSub)    = mean(subMotion(subID,behav(subID).isgo'==behav(subID).iswin'));
    motionIncong(iSub)  = mean(subMotion(subID,behav(subID).isgo'~=behav(subID).iswin'));

end

%% T-tests:

% Action:
fprintf('Go vs. NoGo: \n');
motionTest              = motionGo - motionNoGo;
[~,p,~,stats]           = ttest(motionTest); % index
d                       = mean(motionTest)/std(motionTest);
fprintf('Result is t(%d) = %.03f, p = %.03f, d = %.03f\n',stats.df,stats.tstat,p,d);

% Valence:
fprintf('Win vs. Avoid: \n');
motionTest              = motionWin - motionAvoid;
[~,p,~,stats]           = ttest(motionTest); % index
d                       = mean(motionTest)/std(motionTest);
fprintf('Result is t(%d) = %.03f, p = %.03f, d = %.03f\n',stats.df,stats.tstat,p,d);

% Congruency:
fprintf('Congruent vs. Incongruent: \n');
motionTest              = motionCong - motionIncong;
[~,p,~,stats]           = ttest(motionTest); % index
d                       = mean(motionTest)/std(motionTest);
fprintf('Result is t(%d) = %.03f, p = %.03f, d = %.03f\n',stats.df,stats.tstat,p,d);

% END
