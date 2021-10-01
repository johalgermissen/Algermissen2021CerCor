function [job, data] = TF_prepare_generic_data(job, data)

% [job, data] = TF_prepare_generic_data(job, data)
%
% Add aggregated TF-domain data sets per subject across conditions, per 
% condition across subjects, and both.
% Creates data sets that are only dependent on the conditions used for
% aggregation; NOT on contrast/ band/ channels etc.
% Those data sets are then used in TF_prepare_contrast_data.m
% 
% INPUTS:
% job               = cell, created via update_job.m, needs at least fields 
%   .nSub           = integer, number of subjects
%   .nCond          = integer, number of conditions
%   .validSubs      = numeric vector, subject numbers of valid subjects to
%   be included in analyses
% data              = cell, need at least the following fields:
% 	TFall{iSub}.ValxAct{iCond} = aggregated data per subject per condition,
% 	created via EEGfMRIPav_Cuelocked_8_TF_grouplevel_create.m
%
% OUTPUTS:
% job               = cell, no fields added
% data              = cell, with the following fields added:
%   .tmpSubCond     = average data per subject per condition in 2D cell
%   format.
%   .Sub            = grand average per subject across conditions.
%   .Cond           = grand average per condition across subjects.
%   .mu             = grand average across both conditions and subjects.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel/

fprintf('Prepare generic data sets\n')

%% Average within subjects over conditions (grand average per subject):

data.tmpSubCond = cell(job.nSub,job.nCond);
data.Sub        = cell(job.nSub,1);

for iSub = job.validSubs % 1:job.nSub

    % First recode subjects and conditions in 2D cell:
    for iCondi = 1:job.nCond
        fprintf('Subject %d, Condition %d \n', iSub, iCondi);
        data.tmpSubCond{iSub,iCondi} = data.TFall{iSub}.ValxAct{iCondi}; % for this condition extract this subject
        data.tmpSubCond{iSub,iCondi} = rmfield(data.tmpSubCond{iSub,iCondi},'cumtapcnt');
    end
    
    % Average over conditions (grand average per subject):
    cfg              = [];
    cfg.parameter    = 'powspctrm';
    data.Sub{iSub} = ft_freqgrandaverage(cfg,data.tmpSubCond{iSub,:}); % average over all conditions
    
end

%% Average within conditions over subjects (grand average per condition):
% Use only valid subjects:

data.Cond   = cell(job.nCond,1);

for iCondi = 1:job.nCond % iCondi = 1;
    fprintf('Condition %d \n', iCondi);
    cfg                 = [];
    cfg.parameter       = 'powspctrm';
    data.Cond{iCondi} = ft_freqgrandaverage(cfg,data.tmpSubCond{job.validSubs,iCondi}); % average over subjects
end

%% Overall average of subjects and conditions (reference for time, job.channels, etc.):
% Use output of previous step, thus use only valid subjects:

cfg                   = [];
cfg.parameter         = 'powspctrm';
data.mu               = ft_freqgrandaverage(cfg,data.Cond{:}); % average over conditions

end % end of function.
