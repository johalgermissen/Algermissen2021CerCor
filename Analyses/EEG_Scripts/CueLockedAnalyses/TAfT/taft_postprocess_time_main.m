function taft_postprocess_time_main()

% taft_postprocess_time_main()
% 
% Interactive script.
% Main script to load a previously created TAfT object, align data
% (channels) across subjects, do one-sample-t-test across subjects and plot
% as ER (line-) plot and as topoplot. 
% Use this script to select fMRI and behavioral regressors as well as
% subset of trials; all other settings are to be set in
% taft_postprocess_load_job.m.
% 
% INPUTS (interactive script):
% EEGdomain 	= string, set to 'time' (data in time domain) in this script.
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

EEGdomain   = 'time';

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

% 3) Perform only on selected trials
selTrials   = 'all';

[job, dirs, betas] = taft_postprocess_load_job(EEGdomain,ROIs2use,behav2use,selTrials);

%% ERplot: One single ROI:

rng(20190823) % set random number generator for constant p-values

iROI            = 1;
[sortBetas,~]   = taft_postprocess_time_selectData(job,betas,iROI); % extract data for ROI for selected subjects
sortBetas       = taft_postprocess_FisherZ(sortBetas); % Perform Fisher z-transform

% Select channels:
selChans = {'FCz','Cz'}; % Go/NoGo contrast
% selChans = {'F1','F3','FCz','FC1','FC3','Cz','C1','C3'}; % identified as P2 modulation

% Without saving:
[~,~] = taft_postprocess_time_ERplot(job,dirs,sortBetas,iROI,selChans,2.0,1000,false);
% pause(3)
% close gcf

%% ERplot: Loop over all ROIs:

rng(20190823) % set random number generator for constant p-values

for iROI = 1:length(job.regNames)
    [sortBetas,~]   = taft_postprocess_TF_selectData(job,betas,iROI);
    sortBetas       = taft_postprocess_FisherZ(sortBetas); % Perform Fisher z-transform:
    selChans        = {'FCz','Cz'}; % significant Go/NoGo contrast
    [~,~]           = taft_postprocess_time_ERplot(job,dirs,sortBetas,iROI,selChans,2.0,1000,false);
    pause(3)
    close gcf
end

%% Topoplot: one single ROI:

iROI        = 1;
[~,Tvalues] = taft_postprocess_time_selectData(job,betas,iROI);

taft_postprocess_time_topoplot(job,dirs,Tvalues,iROI,0,0.75,0.10,3,false) % in steps of 100 ms
% taft_postprocess_time_topoplot(job,dirs,Tvalues,iROI,0,0.75,0.05,3,false) % in steps of 50 ms
pause(3)
close gcf

end % end of function.
