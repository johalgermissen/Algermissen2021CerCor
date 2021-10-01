function EEGfMRIPav_CueLocked_1_preICA(rootdir,subVec)

% EEGfMRIPav_CueLocked_1_preICA(rootdir,subVec)
% 
% Read in complete data of one subject (separate blocks), concatenate
% blocks, re-reference, filter, linear baseline correction 
% (complete pre-processing up to ICA).
% 
% INPUTS:
% rootdir           = string, root directory of project
% subVec            = vector of integers, numbers of subjects to process
%
% OUTPUTS:
% saves pre-processed data to disk under dirs.preICA
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% By J.C.Swart, 2016/ J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Preprocessing/

%% Setup.

% Clean workspace and console:
% clear all; close all; clc
% dbstop if error

%% Complete input settings:

% clear all; close all; clc
% dbstop if error

if ~exist('rootdir','var')
    rootdir = preprocessing_set_rootdir(); % '/project/3017042.02';
    fprintf('rootdir unspecified, assume %s\n',rootdir)
end

if ~exist('subVec','var')
    subVec = 1:36; 
    fprintf('subVec unspecified, assume %s\n',strjoin(string(subVec),', '));
end

%% Set directories:

% Add helper files to path:
fprintf('Add helper files to path\n');
addpath(fullfile(rootdir,'Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Helpers'));

% Set directories and parameters:
dirs    = set_dirs(rootdir);
par     = set_par();

%% Load channels to be excluded:

load(fullfile(dirs.data,'channel2repair.mat'));

%% Detect behavioral data files:

behavFiles  = dir(fullfile(dirs.behavior,'*_results.mat'));
behavFiles  = {resultsfile.name};

%% Detect EEG raw data directories:

% Detect all available input folders (all subjects):
tmp         = dir(fullfile(dirs.raw));
tmp         = {tmp.name};
dirList     = tmp(3:end); % delete first two entries
nInput      = length(dirList);
fprintf('Found raw data folders for %d subjects\n',nInput);

% Extract subject numbers from folders available:
subNames    = str2double(extractBetween(dirList,8,10));

% Extract indices of selected subjects:
selIndices  = find(ismember(subNames,subVec));
nSub        = length(selIndices);

%% Loop over subjects, read in data, epoch:

% Sub 8 exception because only 5 blocks (& together; see trials2use below) 

for iSub = selIndices % iSub = 1;
    
    fprintf('Subject %03d: Start!\n',iSub);
    
    %% Create name of output object:
    
    outputFile  = fullfile(dirs.preICA,sprintf('MRIpav_%03d_preICA.mat',iSub));
    fprintf('EEG output file is %s\n',outputFile);
   
    % Check if output file already exists:
    if exist(outputFile,'file')
        warning('Subject %03d: File %s already exists, skip subject\n', iSub, outputFile)
        continue;
    end
    
    %% Load behavioral data:
    
    fprintf('Subject %03d: Load behavioral data\n',iSub);
    load(behavFiles{iSub})
    
    %% Retrieve subject-specific directory:
    
    dirs.sub    = fullfile(dirs.raw, dirList{iSub}); % go to subject directory
    fprintf('EEG input directory is %s\n',dirs.sub);
    
    %% Locate EEG files of all blocks:
    
    EEGfiles    = dir(fullfile(dirs.sub,'*.eeg')); % read all eeg files for this subject
    EEGfiles    = {EEGfiles.name};
    nBlock      = length(EEGfiles);
    fprintf('Subject %03d: Found %d files (blocks)\n',iSub,nBlock);
    
   
    %% Loop over blocks:

    % Initialize empty objects:
    blockData   = cell(nBlock,1);
    trialCount = 0;                     % count trials to explicitly set trial numbers
    
    for iBlock = 1:nBlock % iBlock = 1;
        
        fprintf('Subject %03d Block %d: Start loading \n',iSub,iBlock);

        % ----------------------------------------------------- %
        %% Determine trials that should be contained in respective block:
        
        if iSub==8 % for subject 8: blocks 1 and 2 together
            trials2use = {1:220,221:320,321:430,431:540,541:640}; 
        else
            trials2use = {1:110,111:220,221:320,321:430,431:540,541:640}; % which trial numbers are in which block
        end
        
        % ----------------------------------------------------- %
        %% Detect events:
        
        cfg                     = [];
        cfg.dataset             = fullfile(dirs.sub,EEGfiles{iBlock});
        cfg.trialdef            = [];
        cfg.trialdef.eventvalue = par.epoch.eventCode;
        cfg.trialdef.eventtype  = 'Stimulus';
        cfg.trialdef.prestim    = abs(par.epoch.epochtime(1));
        cfg.trialdef.poststim   = abs(par.epoch.epochtime(2));
        
        fprintf('Subject %03d Block %d: Read events \n',iSub,iBlock);
        cfg.event               = ft_read_event(dirs.fileArrayBlock);
        
        % ----------------------------------------------------- %
        % Define trials:
        fprintf('Subject %03d Block %d: Define trials \n',iSub,iBlock);
        cfg                     = ft_definetrial(cfg);

        % ----------------------------------------------------- %
        % Read in epochs:
        fprintf('>> Subject %03d Block %d: Define epochs \n',iSub,iBlock);
        blockData{iBlock}       = ft_preprocessing(cfg); 

        % ----------------------------------------------------- %
        %% Remove redundant triggers:
        
        % NOTE: subjects might have had one or two additional triggers. Remove cue
        % trigger that follows previous cue trigger within 1,300 ms (cue 
        % presentation is 1,300 ms) and check if trial sequence from behavioural
        % and EEG file than correspond.

        % Check dependent on block: 1 2 4 and 5 are 110 trials long, 3 and
        % 6 only 100 trials long; if violated:
        
        nEEGTrials  = size(blockData{iBlock}.trial,2);
        
        if  (ismember(iBlock,[1 2 4 5]) && nEEGTrials > 110) || (ismember(iBlock,[3 6]) && nEEGTrials > 100)

            % Print subject ID with number of trials.
            fprintf('Subject %03d Block %d had %d trials in EEG data!\n\n',iSub,iBlock,numel(blockData{iBlock}.trial));

            % Retrieve trial sequence of EEG cues:
            EEGCues     = blockData{iBlock}.trialinfo-110; % order of cues presented based on EEG data

            % Check if any two adjacent triggers are less than 1,300 ms
            % apart (cue presentation is 1300 ms):
            tooShortTriggers    = [0; diff(blockData{iBlock}.sampleinfo(:,1))<1300];
            
            % Remove such triggers:
            EEGCues(tooShortTriggers==1) = []; % delete those triggers
            
            % Determine triggers based on behavioral data:
            behavCues   = [prep.seq.learn.stim{1};prep.seq.learn.stim{2}]; % recover entire cue order from behavioural data
            trlIdx      = trials2use{iBlock}; % trial numbers in this block
            behavCues   = behavCues(trlIdx); % order of cues presented in this block based on behavioral data
            
            % Check if trial sequence from behavioural and EEG file correspond now:
            if any(EEGstim-behavCues)
                warning('Number EEG triggers and behavioral trials still do not correspond: Check triggers carefully!\n')
                warning('Skip subject\n')
                continue;
            end

            % Remove trials with extra triggers:
            cfg                 = [];
            cfg.trials          = ~tooShortTriggers;
            blockData{iBlock}   = ft_selectdata(cfg, blockData{iBlock});
            
        end % end remove redundant trigger loop

        %% Explicitly set trial numbers and timing:
        
        % Overwrite trial numbers (trialinfo):
        blockLength                 = length(blockData{iBlock}.trialinfo);
        blockData{iBlock}.trialinfo = [(trialCount+1):(trialCount+blockLength)]';
        trialCount              = trialCount + blockLength;
        
        % Overwrite exact timing (sampleinfo):
        % Add 1,000,000 to sampleinfo with increasing block number:
        blockData{iBlock}.sampleinfo = blockData{iBlock}.sampleinfo + 1000000*(iBlock-1);
        
        % Or just delete timing manually:
        % blockdata{iBlock} = rmfield(blockdata{iBlock},'sampleinfo'); 
    
    end % end iBlock-loop
        
    %% Concatenate blocks:
    
    fprintf('Subject %03d: Concatenate blocks \n',iSub);
    cfg = [];
    % cfg.keepsampleinfo = 'no'; % delete sampleinfo
    data = ft_appenddata(cfg, blockData{1:par.nBlocks}); % concatenate blocks
    
    %% Re-reference EEG channels:

    % Determine invalid channels:
    invalidChan             = cell2mat(channel2repair(iSub,2)); % previously discarded from visual inspection
    validChan               = setdiff(par.chan.recordChan,[par.chan.ECG par.chan.HR par.chan.RESP invalidChan]); % not ECG, not HR, not RESP, not discarded
    % still excludes channel 67 (to-be-recovered in re-referencing)
    
    % Use grand mean reference, recover old reference as new channel FCz:
    
    fprintf('Subject %03d: Re-reference to grand mean, recover channel %s \n',iSub,par.epoch.implicitRef);

    cfg                     = [];
    cfg.reref               = 'yes';
    cfg.implicitref         = par.epoch.implicitRef; % recover old reference as new channel, called FCz; check with data.label{1,67}
    cfg.refchannel          = validChan;
    cfg.refmethod           = 'avg';
    data                    = ft_preprocessing(cfg, data);
    
    %% Band-pass filter:
    
    fprintf('Subject %03d: Band pass filter %.1f Hz to %d Hz \n',iSub,par.epoch.hpf,par.epoch.lpf);

    cfg                     = [];
    cfg.bpfilter            = 'yes';
    cfg.bpfreq              = [par.epoch.hpf par.epoch.lpf]; % 1-15 Hz
    data                    = ft_preprocessing(cfg,data);

    %% Linear baseline correction.

    fprintf('Subject %03d: Apply linear baseline correction from %.03f - %.03f sec.\n',iSub,par.epoch.base4correction(1),par.epoch.base4correction(2));

    cfg                     = [];
    cfg.demean              = 'yes';
    cfg.baselinewindow      = par.base4correction;
    data                    = ft_preprocessing(cfg,data);
    
    %% Save pre-processed data:
    
    fprintf('Subject %03d: Start saving data\n',iSub);
    
    hist.epoch.dirs       = dirs;
    hist.epoch.par        = par;
    
    save(outputFile,'data','hist','-v7.3');
    clear data hist prep trials2use
    
    fprintf('Subject %03d: Finished saving components\n',iSub);

    %% Save complete list of position labels for once:
    
    if iSub == 1
        
        % Save complete data labels for first subject:
        fprintf('Subject %03d: Saving all channel labels\n',iSub);

        posLabels =  blockData{1}.label;
        posLabels{67} = 'FCz'; % add FCz as label (to be sure)
        save(fullfile(dirs.polhemus,'posLabels.mat'),'posLabels')
        
    end
    
end % end iSub-loop.

fprintf('Done :-)\n')

end % end of function.
