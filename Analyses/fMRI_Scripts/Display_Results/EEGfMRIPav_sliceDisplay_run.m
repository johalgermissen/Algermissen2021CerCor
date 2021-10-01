% EEGfMRIPav_sliceDisplay_run.m

% Before using: mind to set directories/paths in slice_display_set_dirs to respective toolboxes needed!

% Directories:
dirs        = []; 
dirs.root   = '/project/3017042.02'; % root directory--needs to be adapted to users' folder structure

addpath(fullfile(dirs.root,'Analyses/fMRI_Scripts/Display_Results'));

%% Stimulus-locked: GLM1 (displayed in main result)

dirs.save           = fullfile(dirs.root,'Log/CueLockedPaperPlots');
job.lockSettings    = 'stimlocked';
job.type            = 'standard';
job.zLim            = [0 3]; % 
job.GLMID           = '1';
% job.GLMID           = '1without6';

% Sagittal:
job.iCope = 1; job.iView = 'sagittal'; job.iSlice = 2; job.cLim = 100; % 
EEGfMRIPav_sliceDisplay_cope(dirs,job)

job.iCope = 2; job.iView = 'sagittal'; job.iSlice = 2; job.cLim = 100; % 
EEGfMRIPav_sliceDisplay_cope(dirs,job)

job.iCope = 3; job.iView = 'sagittal'; job.iSlice = 2; job.cLim = 30; % 
EEGfMRIPav_sliceDisplay_cope(dirs,job)

% Coronal:
job.postFix = '_striatum';
job.iCope = 1; job.iView = 'coronal'; job.iSlice = -10; job.cLim = 100; % 
EEGfMRIPav_sliceDisplay_cope(dirs,job)
job = rmfield(job,'postFix');

job.postFix = '_striatum';
job.iCope = 1; job.iView = 'coronal'; job.iSlice = 4; job.cLim = 100; % 
EEGfMRIPav_sliceDisplay_cope(dirs,job)
job = rmfield(job,'postFix');

job.postFix = '_ACC_JBL';
job.iCope = 3; job.iView = 'coronal'; job.iSlice = 4; job.cLim = 30; % 
EEGfMRIPav_sliceDisplay_cope(dirs,job)
job = rmfield(job,'postFix');

%% Stimulus-locked: GLM1without6 (displayed in supplementary material S01)

dirs.save           = fullfile(dirs.root,'Log/CueLockedPaperPlots');
job.lockSettings    = 'stimlocked';
job.type            = 'standard';
job.zLim            = [0 3]; % 
job.GLMID           = '1without6';

% Sagittal:
job.iCope = 1; job.iView = 'sagittal'; job.iSlice = 2; job.cLim = 100; % 
EEGfMRIPav_sliceDisplay_cope(dirs,job)

job.iCope = 2; job.iView = 'sagittal'; job.iSlice = 2; job.cLim = 100; % 
EEGfMRIPav_sliceDisplay_cope(dirs,job)

job.iCope = 3; job.iView = 'sagittal'; job.iSlice = 2; job.cLim = 30; % 
EEGfMRIPav_sliceDisplay_cope(dirs,job)

% Coronal:
job.postFix = '_striatum';
job.iCope = 1; job.iView = 'coronal'; job.iSlice = -10; job.cLim = 100; % 
EEGfMRIPav_sliceDisplay_cope(dirs,job)
job = rmfield(job,'postFix');

job.postFix = '_striatum';
job.iCope = 1; job.iView = 'coronal'; job.iSlice = 4; job.cLim = 100; % 
EEGfMRIPav_sliceDisplay_cope(dirs,job)
job = rmfield(job,'postFix');

job.postFix = '_ACC_JBL';
job.iCope = 3; job.iView = 'coronal'; job.iSlice = 4; job.cLim = 30; % 
EEGfMRIPav_sliceDisplay_cope(dirs,job)
job = rmfield(job,'postFix');

%% Stimulus-locked: GLM2 (split per correct/incorrect; presented in Supplementary Material S06)

dirs.save           = fullfile(dirs.root,'Log/CueLockedPaperPlots');
job.lockSettings    = 'stimlocked';
job.type            = 'standard';
job.zLim            = [0 3]; % 
job.GLMID           = '2';

% Sagittal:
job.iCope = 4; job.iView = 'sagittal'; job.iSlice = 2; job.cLim = 100; % 
EEGfMRIPav_sliceDisplay_cope(dirs,job)
job.iCope = 5; job.iView = 'sagittal'; job.iSlice = 2; job.cLim = 100; % 
EEGfMRIPav_sliceDisplay_cope(dirs,job)
job.iCope = 6; job.iView = 'sagittal'; job.iSlice = 2; job.cLim = 30; % 
EEGfMRIPav_sliceDisplay_cope(dirs,job)

% Coronal:
job.iCope = 4; job.iView = 'coronal'; job.iSlice = -10; job.cLim = 100; % 
EEGfMRIPav_sliceDisplay_cope(dirs,job)
job.iCope = 4; job.iView = 'coronal'; job.iSlice = 4; job.cLim = 100; % 
EEGfMRIPav_sliceDisplay_cope(dirs,job)
job.iCope = 6; job.iView = 'coronal'; job.iSlice = -10; job.cLim = 100; % 
EEGfMRIPav_sliceDisplay_cope(dirs,job)

%% Stimulus-locked: GLM3 (alpha power stimulus-locked, theta-power response-locked; presented in Supplementary Material S17)

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

% Axial:
job.iCope = 11; job.iView = 'coronal'; job.iSlice = -30; % 
EEGfMRIPav_sliceDisplay_cope(dirs,job)

% Axial:
job.iCope = 10; job.iView = 'axial'; job.iSlice = 30; % 
EEGfMRIPav_sliceDisplay_cope(dirs,job)
job.iCope = 10; job.iView = 'axial'; job.iSlice = 40; % 
EEGfMRIPav_sliceDisplay_cope(dirs,job)
job.iCope = 11; job.iView = 'axial'; job.iSlice = 10; % 
EEGfMRIPav_sliceDisplay_cope(dirs,job)

% END
