function EEGfMRIPav_CueLocked_4_rejICs(rootdir,subVec)

% EEGfMRIPav_CueLocked_4_rejICs(rootdir,subVec)
% 
% Reject identified components (after ICA and visual inspection of plots),
% save data.
% 
% INPUTS:
% rootdir           = string, root directory of project
% subVec            = vector of integers, numbers of subjects to process
%
% OUTPUTS:
% saves data after component rejection to disk under dirs.postICA.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% By J.C.Swart, 2016/ J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Preprocessing/

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

%% Overview of components to reject:

comp2remove                 = { % sID, [IC to remove]
    1, [1 2 3 4 5 6 8 10 11 14 16 17 55]; % heart beat: 1 2 3 4 5 6 8 10 14 16 17; eye-blinks: 11; single channel: 55.
    2, [1 2 3 4 5 6 8 10 12 14 15 18]; % heart beat: 2 3 4 5 6 8 10 12 14 15; eye-blinks: 1; saccades: 18.
    3, [1 2 3 5 6 7 8 9 10 11 12 39 55 61]; % heart beat: 1 2 3 6 9 10 12; eye-blinks: 8; saccades: 5 7 11; cable bundle: 38 55; single channel: 61.
    4, [1 2 3 4 5 6 7 8 9 10 13 14 16 26 33 43 54]; % heart beat: 1 2 5 6 7 8 9 10 13 14 16; eye-blinks: 3 4; saccades: 11; cable bundle: 39, 55; single channel: 26 43 54.
    5, [1 2 3 5 8 10 11 12 15 49]; % heart beat: 1 2 5 8 10 11 12 15 49; eye-blinks: 3; cable bundle: 49.
    6, [1 2 3 4 8 9 10 11 14 20 5]; % heart beat: 1 2 3 4 8 10 11 14 20; eye-blinks: 9; cable bundle: 50.
    7, [1 2 3 4 5 6 7 8 9 10 11 12 18 19 40 47]; % heart beat: 3 4 5 6 7 8 9 10 11 12 19; eye-blinks: 1 2; saccades: 18; cable bundle: 47; single channel: 40.
    8, [1 2 3 5 6 7 8 9 10 52]; % heart beat: 1 3 8 9 10; eye-blinks: 2 7; saccades: 5 6; cable bundle: 52.
    9, [1 2 3 4 5 7 11 12 15 17 18 37 53]; % heart beat: 2 3 4 5 7 11 12; eye-blinks: 1; saccades: 15 17 18; cable bundle: 53; single channel: 37.
    10, [1 2 3 4 5 6 7 10 13 14 45]; % heart beat: 1 3 4 5 6 13 14; eye-blinks: 2 7 10; cable bundle: 45.
    11, [1 2 3 4 5 6 7 8 9 10 13 15 16 18 19 35 58]; % heart beat: 2 3 4 5 6 7 8 9 10 13 15 16 18 19; eye-blinks: 1; cable bundle: 58; single channel: 35.
    12, [1 3 4 5 6 7 9 11 16 29]; % heart beat: 1 3 4 6 9 11; eye-blinks: 2 16; saccades: 5 7; single channel: 29.
    13, [1 2 3 4 5 6 7 8 9 10 12 13 14 18 20 50 54]; % heart beat: 1 2 3 5 7 9 12 18 20; eye-blinks: 4 10 13; saccades: 6 8 14; cable bundle: 54; single channel: 50.
    14, [1 2 3 4 5 6 7 10 11 13 14 18 55]; % heart beat: 1 2 3 5 7 10 11 13 18; eye-blinks: 4 6 14; single channel: 55.
    15, [1 2 3 4 5 6 8 13 14 15 20 61]; % heart beat: 1 2 3 4 5 6 13 14 15 20; eye-blinks: 8; cable bundle: 61.
    16, [1 2 3 4 5 6 8 10 12 14 17 18 19]; % heart beat: 1 2 3 5 6 10 15 17 18 19; eye-blinks: 4; saccades: 8 12; single channel: 25.
    17, [1 3 4 5 7 8 9 12 16 17 18 37 42]; % heart beat: 3 4 5 7 8 9; eye-blinks: 1 17; saccades: 12; cable bundle: 37; single channel: 16, 18, 42.
    18, [1 2 3 4 5 7 10 11 12 15 19 36 51 52]; % heart beat: 2 3 4 5 7 12 15 19; eye-blinks: 1; saccades: 10; cable bundle: 36 52; single channel: 51.
    19, [1 2 3 4 5 6 7 10 12 17 18 19 52]; % heart beat: 1 4 5 6 7 12 17 18 19; eye-blinks: 2 3; saccades: 8 10; cable bundle: 52.
    20, [1 2 3 4 5 6 7 8 9 11 18 52]; % heart beat: 2 4 6 18; eye-blinks: 1 3 7 8 ; saccades: 5 9 11; cable bundle: 52.
    21, [1 3 4 6 7 8 11 12 15 18 35 49]; % heart beat: 1 4 6 7 8 15 18; eye-blinks: 3; saccades: 11 12; single channel: 35 49.
    22, [1 2 3 4 5 12 13 18 44]; % heart beat: 2 3 5 12 13 18; eye-blinks: 1 4; single channel: 44.
    23, [1 2 3 4 5 6 8 9 10 11 15 16 17 19 60 61]; % heart beat: 1 2 4 5 6 8 9 10 11 15 16 17 19; eye-blinks: 3; cable bundle: 60; single channel: 61.
    24, [1 2 3 4 6 8 10 57]; % heart beat: 2 3 4 6 8; eye-blinks: 1; saccades: 10; cable bundle: 57.
    25, [1 2 3 4 5 6 7 8 9 10 11 12 13 14 16 17 18 20 49]; % heart beat: 1 2 3 4 5 8 9 10 11 12 13 20; eye-blinks: 6 14 18; saccades: 7 16; single channel: 17 49.
    26, [1 2 3 4 5 7 12 15]; % heart beat: 1 3 4 5 7 12 15; eye-blinks: 2.
    27, [1 2 3 4 5 7 8 9 23 47 49 53]; % heart beat: 1 2 3 4 5 9 13; eye-blinks: 7; saccades: 8 23; single channel: 47 49 53.
    28, [1 2 3 4 5 6 8 9 12 14 15 20 46]; % heart beat: 2 3 4 5 8 9 14 15 20; eye-blinks: 1; saccades: 6 12; single channel: 46.
    29, [1 2 3 4 6 7 8 9 14 15 16 17 18 20 21 22]; % heart beat: 2 3 6 7 8 9 14 15 16 17 18 20 21 22; eye-blinks: 1 4.
    30, [2 3 4 5 6 7 8 9 10 11 14 16 17 19 20]; % heart beat: 2 3 4 5 6 7 8 9 11 14 16 17 19 20; eye-blinks: 10.
    31, [1 2 3 4 5 7 14 22 51 52]; % heart beat: 2 3 5 7; eye-blinks: 1; saccades: 4 14 22.
    32, [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16]; % heart beat: 1 2 3 7 10 11 12 13; eye-blinks: 14; saccades: 4 5 6 8 9 15 16.
    33, [1 2 3 4 5 6 7 8 10 11 19 20 44 53]; % heart beat: 1 3 4 5 6 7 8 10 19 20; eye-blinks: 2; cable bundle: 11; single channel: 44 53.
    34, [1 2 3 4 5 6 7 8 9 10 13 34]; % heart beat: 2 3 5 6 7 8 9 10; eye-blinks: 1 4; saccades: 13; cable bundle: 34.
    35, [1 2 3 4 5 6 7 8 25 48 54]; % heart beat: 2 3 5 7; eye-blinks: 1 4 6; saccades: 8; cable bundle: 25; single channel: 48 54.
    36, [1 2 3 4 5 7 8 9 12 13 15 19 44 47]; % heart beat: 1 2 3 5 7 8 9 12 15 19; eye-blinks: 4; saccades: 13; cable bundle: 44; single channel: 57.
    };

%% Check how many got removed:

nSub    = 36;
nComp   = nan(nSub,1); % initialize
for iSub = 1:nSub   
    nComp(iSub) = length(comp2remove{iSub,2});    
end

fprintf('Number of rejected components: M = %.03f, SD = %.03f, range %d - %d \n',mean(nComp),std(nComp),min(nComp),max(nComp));
% Number of rejected components: M = 12.944, SD = 2.661, range 8 - 19 

fprintf('Number of rejected components: %s\n',num2str(sort(nComp,'descend')',' %d')); % pretty continuous at high end, makes not much sense to set arbitrary cut-off...
% Number of rejected components: 19 17 17 17 16 16 16 16 15 14 14 14 14 13 13 13 13 13 13 13 12 12 12 12 12 12 11 11 11 10 10 10 10  9  8  8

%% Detect input files:

% Detect all available input folders (all subjects):
tmp         = dir(fullfile(dirs.ICA,'*.mat'));
fileList    = {tmp.name};
nInput      = length(fileList);
fprintf('Found data files from %d subjects\n',nInput);

% Extract subject numbers from folders available:
subNames    = str2double(extractBetween(fileList,8,10));

% Extract indices of selected subjects:
selIndices  = find(ismember(subNames,subVec));

%% Loop over subjects.

for iSub = selIndices % iSub = 1;
    
    %% Create output file name:
    
    outputFile      = fullfile(dirs.postICA, sprintf('MRIpav_%03d_postIC.mat',iSub));
    fprintf('EEG output file is %s\n',outputFile);
    
    % Check if output file already exists:
    if exist(outputFile,'file')
        warning('Subject %03d: File %s already exists, skip subject\n', iSub, outputFile)
        continue;
    end

    %% Create input file name:
    
    inputFile       = fullfile(dirs.ICA, fileList{iSub}); % go to subject directory
    fprintf('EEG input file is %s\n',inputFile);
    
    %% Load data:
    
    fprintf('Subject %03d: Loading data ... (can take a few seconds)\n',iSub);
    load(inputFile) % 
    fprintf('Subject %03d: Loading completed\n',iSub);
    
    
    %% Remove selected ICs:
    
    rejComps = comp2remove{[comp2remove{:,1}]==iSub,2};
    fprintf('Subject %03d: Reject components %s \n',iSub,strjoin(string(rejComps),', '));

    cfg                     = [];
    cfg.component           = rejComps;
    data                    = ft_rejectcomponent(cfg, comp);
   
    %% Save data: 

    fprintf('Subject %03d: Start saving data\n',iSub);
    
    hist.postICA.dirs       = dirs;
    hist.postICA.par        = par;
    
    save(outputFile,'data','hist','-v7.3');
    clear comp data hist

    fprintf('Subject %03d: Finished saving data\n',iSub);
    
end % end iSub-loop.

fprintf('Done :-)\n');

end % end of function.
