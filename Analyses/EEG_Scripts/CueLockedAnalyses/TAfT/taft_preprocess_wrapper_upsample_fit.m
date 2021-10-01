function out = taft_preprocess_wrapper_upsample_fit(job)

% out = taft_preprocess_wrapper_upsample_fit(job)
% 
% Loop over ROIs, loop over blocks, call function for fitting HRF
% amplitude, reshape into design matrix.
%
% INPUTS:
% job               = structure with settings for creating TAfT object, specifically:
% .HRFtype          = string, whether to perform estimation of HRF amplitude for each trial separately ('trial') or in a GLM featuring all trials of a given block ('block').
% .ROIs{iROI}       = further settings/file names to use for particular ROI for particular block of particular subject, specifically:
%   .ROIname        = string, name of particular ROI. 
%   .ROIdef{iBlock}.betafMRIfile= string, file to which to save trial-by-trial HRF amplitude estimates (optional).
%   .ROIdef{iBlock}.fitType     = string, type of turning volume-by-volume data into trial-by-trial data, either 'HRF' (default) or 'avg' (for mp_Fellner regressor).
% .save             = Boolean, whether to save estimates of trial-by-trial HRF amplitude to disk.
% .subID            = numeric scalar, subject ID.
% .TR               = numeric scalar, repetition time of MRI sequence (used for upsampling volumes) in seconds.
% .ups              = numeric scalar, factor by which to upsample MRI sequence.
%
% Output:
% out.ROIs{iROI}    = cell with following settings per ROI:
% .ROIname          = string, ROI name.
% .ROIdef{iBlock}   = numeric vector, trial-by-trial HRF amplitude in given ROI for given block.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Adapted from Tobias Hauser (https://github.com/tuhauser/TAfT).
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/CueLockedAnalyses/TAfT/

%% Check settings in job:

if isempty(job.ups) || job.ups < 0 || ~isnumeric(job.ups)
    error('Upsampling factor must be a positive number'); end
if isempty(job.TR) || job.TR < 0 || ~isnumeric(job.TR)
    error('Repetition Time must be a positive number'); end

%% For each ROI: upsample volumes, fit trial-by-trial betas:

out = []; % initialize empty output

% Loop through ROIs:
for iROI = 1:length(job.ROIs)

    % Carry forward ROI name to out object:
    out.ROIs{iROI}.ROIname  = job.ROIs(iROI).ROIname;
    
    % Loop through blocks:
    for iBlock = 1:length(job.ROIs(iROI).ROIdef) % for each block % iBlock = 6
        
        fprintf('Subject %03d Block %d: Start processing fMRI data\n',job.subID,iBlock);

        % Upsample ROI-data:
        [job,data,data_ups,i_ons] = taft_preprocess_filter_upsample_epoch(job,iROI,iBlock);
        
        % Fit HRF per trial:
        if strcmp(job.ROIs(iROI).ROIdef(iBlock).fitType,'HRF')
            
            fprintf('Fit HRF per trial, get trial-by-trial betas\n');
            
            % Separate GLM per trial:
            if strcmp(job.HRFtype,'trial')
                out.ROIs{iROI}.ROIdef{iBlock}   = taft_preprocess_fit_HRF_trial(job,data);
            
            % One single GLM across all trials of entire block:
            elseif strcmp(job.HRFtype,'block')
                out.ROIs{iROI}.ROIdef{iBlock}   = taft_preprocess_fit_HRF_block(job,data_ups,i_ons);
            
            else
                 error('Unknown GLM type')
            end
            
        % Or just take average of epoch (for volume-by-volume nuisance regressors):
        elseif strcmp(job.ROIs(iROI).ROIdef(iBlock).fitType,'avg')

            out.ROIs{iROI}.ROIdef{iBlock} = nanmean(data,2)';
            
        else
            error('Unknown fitting type')
        end
        
    end % end iBlock
    
    % Save trial-by-trial HRF amplitudes as separate file:
    if job.save
        outfile = job.ROIs(iROI).ROIdef(iBlock).betafMRIfile{:}; %
        save(outfile,'out');
    end
    
end % end iROI

end % end of function.
