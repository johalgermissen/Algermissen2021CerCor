function EEGfMRIPav_0_Polhemus_Template(rootdir)

% EEGfMRIPav_0_Polhemus_Template(rootdir)
% 
% Interactive script (visually inspect channel locations).
% Load polhemus data of all (valid) subjects, create mean as template for
% invalid subjects.
% 
% INPUTS:
% rootdir           = string, root directory of project
%
% OUTPUTS:
% saves template to disk under dirs.polhemus
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
    rootdir = '/project/3017042.02';
    fprintf('rootdir unspecified, assume %s\n',rootdir)
end

%% Set directories:

% Add helper files to path:
fprintf('Add helper files to path\n');
addpath(fullfile(rootdir,'Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Helpers'));

% Set directories and parameters:
dirs    = set_dirs(rootdir);
par     = set_par();

%% Read position labels:
% created when doing initial pre-processing in EEGfMRIPav_CueLocked_1_preICA.m

fprintf('Load overall position labels\n');

load(fullfile(dirs.polhemus,'posLabels.mat'));
% object posLabels contains all 67 channels collected

%% Check for existing Polhemus files:

fileList             = dir(fullfile(dirs.polhemus,'*.pos'));
fileList             = {fileList.name};
nSub                 = length(fileList);

fprintf('Found polhemus files from %d subjects \n',nSub);

%% Extract existing subject numbers:

subNames = nan(nSub,1);
for iSub = 1:nSub
    subNames(iSub) = str2double(fileList{iSub}((end-6):(end-3)));
end

fprintf('Subject numbers are %s\n',strjoin(string(subNames),', '));

%% Load all existing Polhemus files into one object:

subIdx = 0; % initialize index
allElec = struct([]); % initialize object for all data files

for iSub = 1:nSub
    subID = subNames(iSub);
    fprintf('Subject %03d: Start reading polhemus file\n', subID)
    subIdx = subIdx + 1; % increment
    allElec{subIdx} = ft_read_sens(fullfile(dirs.polhemus,sprintf('MRIpav_001_%03d.pos',subID)));
    fprintf('Found %d channels\n',size(allElec{subIdx}.label,1));
end
fprintf('Finished reading Polhemus data\n');

% All subjects should have 70 channels

% sub 04 has 71 channels
% sub 05 has 69 channels

%% Plot locations for all subjects:

for iSub = 1:nSub
    
    sens = allElec{iSub}; % select subject
    sens.label(1:par.chan.nChan) = cell(posLabels); % change numbers into channel names
    
    % Plot:
    figure('name',sprintf('Subject %03d',subNames(iSub)),'Position',[50 100 800 600])
    sensplot = ft_plot_sens(sens,'label','on','elec',true);
    w = waitforbuttonpress; % wait for click
    close gcf; % close plot again; loop to next subject

end

%% Subject with concerning pattern:

fprintf('Inspect problematic subjects\n');

% ----------------------------------------------- %
% Subject 4:
fprintf('Inspect subject 4\n');
sens = allElec{4};
sens.label(1:67) = cell(posLabels);
ft_plot_sens(sens,'label','on','elec',true);
% some channels (FT9, AF3, 02, F5, Oz, POz, P2, PO8, PO4) extremely off
% --> not clear how to resolve the issue...? Use overall template

% ----------------------------------------------- %
% Subject 5:
fprintf('Inspect subject 5\n');
sens = allElec{5};
sens.label(1:67) = cell(posLabels);
ft_plot_sens(sens,'label','on','elec',true);
% Recorded labels are normally 1:67, nasion, left, right
% 64 = CPz; 65 = HR; 66 = RESP; 67 = FCz
% Subject 5 has only 1:66, so one missing
% --> CPz (64) correct, but FCz (67) incorrect --> missed HR (65) or RESP (66)?
% --> HR visible at very left (next to C3), so probably forgot RESP
% --> FCz actually on position 66 

% ----------------------------------------------- %
% Correct sens.label for subject 5, plot again:
sens = allElec{5};
sens.label(1:66) = cell(posLabels([1:65,67]));
ft_plot_sens(sens,'label','on','elec',true);
% --> corrected, looks better now

%% Average all channel positions:
 
fprintf('Correct sens.label for all subjects\n');

% Overwrite channel labels: 

for iSub = 1:length(allElec)
    if iSub==5
        allElec{iSub}.label(1:66) = cell(posLabels([1:65,67]));
    else
        allElec{iSub}.label(1:67) = cell(posLabels);
    end
end

% Select valid subjects:
subIdx = [1:3,6:33]; % only subjects with all channels (without 4 and 5)

% Average electrodes:
fprintf('Average channel positions of all valid subjects (%s)\n',strjoin(string(subIdx),', '));
avgGrand = ft_average_sens([allElec{subIdx}]);

%% Compute mean absolute deviation from template for each subject:

fprintf('Compute deviation of each subject from template\n');

devElec         = cell(length(subNames),1);
meanDevElec     = nan(nSub,3); % 

nValidChannels  = 66;

for iSub = 1:length(subNames)
    % Compute average minus each subjects' individual positions:
    devElec{iSub}       = avgGrand.chanpos(1:nValidChannels,:) - allElec{iSub}.chanpos(1:nValidChannels,:);
    meanDevElec(iSub,:) = mean(abs(devElec{iSub})); % mean absolute deviation
end

% Check (expect subject 4 to be way off the mean):
threshold   = 2;
fprintf('Deviation > %d for each subject (rows) for each coordinate (columns\n'.threshold);
(meanDevElec > threshold) % > 1 only for subject 4
% meanDevElec % see strong average deviation for subject 4

%% Save template:

fprintf('Save template (average of valid subjects)\n');

save(fullfile(dirs.polhemus,'PolhemusAvg.mat'),'avgGrand');

fprintf('Done :-)\n')

end % end of function.