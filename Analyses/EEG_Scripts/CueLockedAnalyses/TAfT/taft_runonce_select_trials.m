function taft_runonce_select_trials()

% taft_runonce_select_trials()
%
% Load behavior, retrieves indices of trials on which certain relevant
% behaviors/ conditions happen (Go/ NoGo, Win/ Avoid, congruent/
% incongruent), save indices as separate files under TAfT_Betas/selTrials.
% Mind setting the root directory in taft_set_rootdir().
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Adapted from Tobias Hauser (https://github.com/tuhauser/TAfT).
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/CueLockedAnalyses/TAfT/

%% Settings:

nSub = 36; % number subjects

%% Directories:

dirs.root     = taft_set_rootdir(); % /project/3017042.02

dirs.target   = fullfile(dirs.root,'Log/EEG/CueLockedResults/TAfT_Betas/selTrials');
% Create target directory if not existing yet:
if ~exist(dirs.target,'dir'); mkdir(dirs.target); end

%% Loop through subjects:

for iSub = 1:nSub
    
    fprintf('Start subject %03d\n',iSub);
    
    %% Load behavior:    
    
    job.behavFile   = fullfile(dirs.root,'/Log/Behavior/Data_beh_mat',...
    sprintf('3017042.02_emmvdij_%03d_001_results.mat',iSub));    

    out             = taft_preprocess_load_behavior(job);    
    
    %% Select indices of where relevant behavior/ condition occurs, save:
    
    % Go trials:
    fprintf('Go actions\n')
    selTrials   = find(out.isgo);
    save(fullfile(dirs.target,sprintf('selTrials_GoTrials_Sub%03d.mat',iSub)),'selTrials');
    
    % NoGo trials:
    fprintf('NoGo actions\n')
    selTrials   = find(~out.isgo);
    save(fullfile(dirs.target,sprintf('selTrials_NoGoTrials_Sub%03d.mat',iSub)),'selTrials');
    
    % Win trials:
    fprintf('Win cues\n')
    selTrials   = find(out.iswin);
    save(fullfile(dirs.target,sprintf('selTrials_WinTrials_Sub%03d.mat',iSub)),'selTrials');
    
    % Avoid trials:
    fprintf('Avoid cues\n')
    selTrials   = find(~out.iswin);
    save(fullfile(dirs.target,sprintf('selTrials_AvoidTrials_Sub%03d.mat',iSub)),'selTrials');
    
    % Congruent trials:
    fprintf('Congruent actions\n')
    selTrials   = find(out.isgo == out.iswin);
    save(fullfile(dirs.target,sprintf('selTrials_CongTrials_Sub%03d.mat',iSub)),'selTrials');
    
    % Incongruent trials:
    fprintf('Incongruent actions\n')
    selTrials   = find(out.isgo ~= out.iswin);
    save(fullfile(dirs.target,sprintf('selTrials_IncongTrials_Sub%03d.mat',iSub)),'selTrials');

    % Win2Go trials:
    fprintf('Go actions to Win cues\n')
    selTrials   = find(out.isgo == 1 & out.iswin == 1);
    save(fullfile(dirs.target,sprintf('selTrials_Win2GoTrials_Sub%03d.mat',iSub)),'selTrials');

    % Win2NoGo trials:
    fprintf('NoGo actions to Win cues\n')
    selTrials   = find(out.isgo == 0 & out.iswin == 1);
    save(fullfile(dirs.target,sprintf('selTrials_Win2NoGoTrials_Sub%03d.mat',iSub)),'selTrials');

    % Avoid2Go trials:
    fprintf('Go actions to Avoid cues\n')
    selTrials   = find(out.isgo == 1 & out.iswin == 0);
    save(fullfile(dirs.target,sprintf('selTrials_Avoid2GoTrials_Sub%03d.mat',iSub)),'selTrials');

    % Avoid2NoGo trials:
    fprintf('NoGo actions to Avoid cues\n')
    selTrials   = find(out.isgo == 0 & out.iswin == 0);
    save(fullfile(dirs.target,sprintf('selTrials_Avoid2NoGoTrials_Sub%03d.mat',iSub)),'selTrials');
    
end % end iSub.

end % end of function.
