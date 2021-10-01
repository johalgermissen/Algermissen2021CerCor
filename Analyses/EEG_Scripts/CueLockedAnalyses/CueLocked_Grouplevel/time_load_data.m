function [job, data] = time_load_data(job)

% [job, data] = time_load_data(job) 
%
% Loads time-domain EEG data based on job settings.
%
% INPUTS:
% job               = cell, needs at least fields:
%   inputFileName   = string, name of file to be loaded, created by 
%   TF_update_job.m.
%   dirs.timegroup  = string, directory where data aggregated per condition
%   per subject are stored.
%
% OUTPUTS:
% job               = cell, field 'nSub' gets added
% data              = cell, loaded data in Fieldtrip format.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel/

%% Load EEG data:

fprintf('Start loading %s\n',job.inputFileName);
data = load(fullfile(job.dirs.timegroup,job.inputFileName));
fprintf('Finished loading %s\n',job.inputFileName);   

%% Determine number of subjects:

job.nSub = size(data.ERPdata,1); % update subject number

end % end of function.
