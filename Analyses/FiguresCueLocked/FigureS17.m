% Figure02.m

% Plots for Figure S17 in supplementary material.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/FiguresCueLocked

% Set root directory:
rootdir = figures_set_rootdir(); % '/project/3017042.02';

%% Figure S17 A-B. SLICE DISPLAY MAPS:

addpath(fullfile(rootdir,'Analyses/fMRI_Scripts/Display_Results/'));

dirs.root           = rootdir;
dirs.save           = fullfile(dirs.root,'Log/CueLockedPaperPlots');
job.lockSettings    = 'stimlocked';
job.type            = 'standard';
job.cLim            = 15;
job.zLim            = [0 3]; % 
job.GLMID           = '3';

% Sagittal:
job.iCope = 10; job.iView = 'sagittal'; job.iSlice = -38; % 
EEGfMRIPav_sliceDisplay_cope(dirs,job)
job.iCope = 10; job.iView = 'sagittal'; job.iSlice = 56; % 
EEGfMRIPav_sliceDisplay_cope(dirs,job)
job.iCope = 11; job.iView = 'sagittal'; job.iSlice = 4; % 
EEGfMRIPav_sliceDisplay_cope(dirs,job)

% Coronal:
job.iCope = 11; job.iView = 'coronal'; job.iSlice = -30; % 
EEGfMRIPav_sliceDisplay_cope(dirs,job)

% Axial:
job.iCope = 10; job.iView = 'axial'; job.iSlice = 30; % 
EEGfMRIPav_sliceDisplay_cope(dirs,job)
job.iCope = 10; job.iView = 'axial'; job.iSlice = 40; % 
EEGfMRIPav_sliceDisplay_cope(dirs,job)
job.iCope = 11; job.iView = 'axial'; job.iSlice = 10; % 
EEGfMRIPav_sliceDisplay_cope(dirs,job)

% END.
