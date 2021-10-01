function EEGfMRIPav_extract_rawdata()

% EEGfMRIPav_extract_rawdata()
% 
% Reads the .mat raw data files, extracts relevant fields, saves as .csv
% per subject per session.
%
% Mind setting root directory and target directory for saving.
% 
% INPUT:
% none.
%
% Output:
% Saves .csv files.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

%% Set directories:

fprintf('Set directories\n')

rootdir     = '/project/3017042.02'; % root directory--needs to be adapted to users' folder structure

behavdir    = fullfile(rootdir,'Log/Behavior');
matdir      = fullfile(behavdir,'Data_beh_mat'); % where the data comes from
csvdir      = fullfile(behavdir,'Data_beh_csv'); % where the data goes to

%% Extract subject numbers:

fprintf('Detect files and extract subject numbers\n');
subtmp  = dir(fullfile(matdir,'*_results*')); % select proper files
subtmp  = {subtmp.name}; % keep only column "name"
nSub    = length(subtmp);
subList = nan(nSub,1);

fprintf('Found %d subjects\n',nSub)

for iSub = 1:nSub % 1 till number of elements in subtmp
    subList(iSub) = str2double(string(extractBetween(subtmp{iSub},'dij_','_001_results')));  % extracts the actual "sub" number of each subject --> check which ones found
end
subjects = sort(subList); % create a final subjectlist

fprintf('Found subject numbers are %s\n',num2str(subjects',' %d'));

%% To-be-extracted data:

% load(fullfile(dataPath_mat,'3017042.02_emmvdij_001_001_results.mat'))
% prep.seq.learn
%         stim: {[320x1 double]  [320x1 double]} % number of stimulus presented in each trial
%     feedback: {2x8x3 cell} % each 2 (round) x 8 (stimulus type) x 3 (response type) matrices with 40 entries of 1 or 0: indicate for each time (40) a stimulus (8) is encountered whether given response (3) will result in valid or invalid feedback (?)
%          SFI: {2x1 cell} % outcome presentation onset (?)
%         resp: {[320x1 double]  [320x1 double]} % correct response key (0 = NoGo, 37 = Left Go, 39 = Right Go)
%          ITI: {[320x1 double]  [320x1 double]} % Inter-trial-interval (1.25, 1.5, 1.75, 2)
% results.learn{1} % 320 trials
%           RT: [320x1 double] % RT per trial, NaN for NoGo
%          acc: [320x1 double] % 1 for correct, 0 for incorrect response
%     response: [320x1 double] % pressed response key (0 for NoGo, 37 = Left Go, 39 = Right Go)
%      outcome: [320x1 double] % 1 for reward, 0 for neutral, -1 for punishment
%           go: [320x1 double] % 1 for (any) Go, 0 for NoGo

%% Load data:

fprintf('Loop over subjects and sessions to extract data\n');

for iSub = 1:nSub % for each subject
    
    data = nan(320,10); % initialize output object
    
    fprintf('Subject %03d: Load data\n',iSub);
    load(fullfile(matdir,subtmp{iSub})); % load input data
    
    for iSes = 1:2 % for each session
        
        data(:,1) = repelem(subjects(iSub),320); % subject number
        data(:,2) = repelem(iSes,320); % session number
        data(:,3) = 1:320; % trial number
        data(:,4) = prep.seq.learn.stim{iSes}; % stimulus presented
        data(:,5) = prep.seq.learn.resp{iSes}; % correct response
        data(:,6) = results.learn{iSes}.RT; % RT per trial, NaN for NoGo
        data(:,7) = results.learn{iSes}.acc; % 1 for correct, 0 for incorrect
        data(:,8) = results.learn{iSes}.response; % pressed response key (0 for NoGo, 37 = Left Go, 39 == Right Go)
        data(:,9) = results.learn{iSes}.outcome; % 1 for reward, 0 for neutral, -1 for punishment
        data(:,10) = results.learn{iSes}.go; % 1 for (any) Go, 0 for NoGo
        
        % Target file name:
        targetFile = sprintf('EEGfMRIPav_%d_%d.csv',subjects(iSub),iSes); % name of output file

        % Save data:
        fprintf('Subject %03d: Save data under \n', iSub, outputfile);
        csvwrite(fullfile(csvdir,targetFile),data); % write file as .csv

    end % end iSes
end % end iSub

end % end of function.