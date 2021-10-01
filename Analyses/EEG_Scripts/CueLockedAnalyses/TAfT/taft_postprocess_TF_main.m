function taft_postprocess_TF_main()

% taft_postprocess_TF_main()
% 
% Interactive script.
% Main script to load a previously created TAfT object, align data
% (channels) across subjects, do one-sample-t-test across subjects and plot
% as TF plot and as topoplot. 
% Use this script to select fMRI and behavioral regressors as well as
% subset of trials; all other settings are to be set in
% taft_postprocess_load_job.m.
% 
% INPUTS (interactive script):
% EEGdomain 	= string, set to 'TF' (data in TF domain) in this script.
% ROIs2use      = cell, one or more string(s) specifying the file names
%               of fMRI ROIs data files to include in design matrix. For any volume-by-volume regressor (also nuisance regressor).
% behav2use       cell, one or more string(s) specifying the behavioral variables 
%               to include in design matrix (optional). For any trial-by-trial regressor.
% selTrials     = string, sub-selection of trials to use (optional).
%
% OUTPUTS:
% Will print results to console and create plots.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Adapted from Tobias Hauser (https://github.com/tuhauser/TAfT).
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/CueLockedAnalyses/TAfT/

%% Load EEG-fMRI regression weights, specify directories:

EEGdomain   = 'TF';

% ROIs2use    = {}; % empty
ROIs2use    = {'GLM1StriatumAction','GLM1CingulateAnteriorAction','GLM1LeftMotorHand','GLM1RightMotorHand','GLM1vmPFCValenceMan'}; % final paper
% With ACC Valence contrast (instead of Action contrast) ROI:
% ROIs2use    = {'GLM1StriatumAction','GLM1CingulateAnteriorValence','GLM1LeftMotorHand','GLM1RightMotorHand','GLM1vmPFCValenceMan'};
% With striatum split into regions coding positive/negative valence effect:
% ROIs2use    = {'GLM1LeftPutamenValence','GLM1MedialCaudateValence','GLM1CingulateAnteriorAction','GLM1LeftMotorHand','GLM1RightMotorHand','GLM1vmPFCValenceMan'};
% With relative displacement (realignment parameters summary metric based on Fellner et al., 2016):
% ROIs2use    = {'GLM1StriatumAction','GLM1CingulateAnteriorAction','GLM1LeftMotorHand','GLM1RightMotorHand','GLM1vmPFCValenceMan','mp_Fellner'};

% 2) Behavioral regressors:
% behav2use   = {}; % none
behav2use   = {'isgo','iswin','isW2G'}; % main effects and interaction --> final paper

% 3) Perform only on selected trials:
selTrials   = 'all';

[job, dirs, betas] = taft_postprocess_load_job(EEGdomain,ROIs2use,behav2use,selTrials);

%% TF plot: One single ROI:

rng(20190822) % set random number generator for constant p-values
job.invalidSubs = [1 11 15 19 21 25]; % invalid co-registrations and spike outliers (> 5)

iROI            = 1; % iROI to be tested/ plotted
[sortBetas,~]   = taft_postprocess_TF_selectData(job,betas,iROI); % select data, align channels across subjects
sortBetas       = taft_postprocess_FisherZ(sortBetas); % Perform Fisher z-transform:
selChans        = {'FCz','Cz'}; % significant Go/NoGo contrast
[~,~]           = taft_postprocess_TF_TFplot(job,dirs,sortBetas,iROI,selChans,2,1000,2,false); % job,dirs,sortBetas,iROI,selChans,thresh,nP,zlim,isNormalize,isSave
pause(1)
close gcf

%% TF plot: Loop over all ROIs:

rng(20190823) % set random number generator for constant p-values

for iROI = 1:length(job.regNames)
    [sortBetas,~]   = taft_postprocess_TF_selectData(job,betas,iROI);
    sortBetas       = taft_postprocess_FisherZ(sortBetas); % Perform Fisher z-transform:
    selChans        = {'FCz','Cz'}; % significant Go/NoGo contrast
    [~,~]           = taft_postprocess_TF_TFplot(job,dirs,sortBetas,iROI,selChans,2,1000,2,false); % iROI, selChans, thresh,nP
    pause(3)
    close gcf
end

%% Topoplot: one single ROI, for all frequency bands:

iROI                = 1; % select ROI
[~,Tvalues]         = taft_postprocess_TF_selectData(job,betas,iROI); % extract data for ROI for selected subjects
taft_postprocess_TF_topoplot(job,dirs,Tvalues,iROI,[1 4 8 13],[4 8 13 15]) % ROI, startFreq, endFreq
pause(3)
% close all

%% Topoplot: one single ROI, separate command for each frequency band:

% Settings:
startTime = 0; endTime = 1.1; step = 0.1; % stimulus-locked
% startTime = -1.0; endTime = 0.4; step = 0.1; % response-locked

% Selected ROI:
iROI = 1;

% Data selection:
[~,Tvalues] = taft_postprocess_TF_selectData(job,betas,iROI); % extract data for ROI for selected subjects

% Each band:
taft_postprocess_TF_topoplot(job,dirs,Tvalues,iROI,1,4,startTime,endTime,step,3,false) % iROI, startFreq, endFreq
taft_postprocess_TF_topoplot(job,dirs,Tvalues,iROI,4,8,startTime,endTime,step,3,false) % iROI, startFreq, endFreq
taft_postprocess_TF_topoplot(job,dirs,Tvalues,iROI,8,13,startTime,endTime,step,3,true) % iROI, startFreq, endFreq
taft_postprocess_TF_topoplot(job,dirs,Tvalues,iROI,[1 4 8 13],[4 8 13 15],0,1.1,0.1,3,false) % ROI, startFreq, endFreq
pause(5)
close all

%% Topoplot: Loop through ROIs, all frequency bands:

for iROI = 1:5 % iROI = 1;
    
    [~,Tvalues] = taft_postprocess_TF_selectData(job,betas,iROI); % extract data for ROI for selected subjects
    taft_postprocess_TF_topoplot(job,dirs,Tvalues,iROI,[1 4 8 13],[4 8 13 15],0,1.1,0.1,3,false) % ROI, startFreq, endFreq, startTime, endTime, steps, zlim, isSave
    close all

end

end % end of function.
