function betas = taft_preprocess_combine_EEG_fMRI(job,X,Y)

% betas = taft_preprocess_combine_EEG_fMRI(job,X,Y)
%
% Take a design matrix X of trial-by-trial BOLD HRF estimates and a 
% matrix Y including trial-by-trial EEG data, 
% add behavioral (nuisance) regressors to X, compute the pseudo-inverse,
% perform a multiple linear regression across trials for each
% channel-time(-frequency) bin, output regression weights for each
% regressor. Separately per subject.
%
% INPUTS:
% job           = a structure with settings to be specified via taft_preprocess_initialize_job, specifically:
% .add_intercept= Boolean, whether to include intercept to design matrix or not. 
% .behavReg 	= cell of strings, behavioral regressors to include as loaded and specified in taft_preprocess_load_behavior (optional).
% .dBconversion = Boolean, whether to convert outcome variable (EEG data) to decibel (10*log10(Y)) or not.
% .EEGdomain 	= string, either 'TF' (data in TF domain) or 'time' (data in time domain).
% .goodTrlIdx 	= numeric vector, indices of trials to include (after EEG trial rejection).
% .demeanX      = Boolean, whether to demean design matrix X or not.
% .demeanY      = Boolean, whether to demean outcome variable (EEG data) or not.
% .regType      = string, type of regression over trials performed, either 'pinv', 'fitlm', 'robustfit'
% .subID        = numeric scalar, subject ID
% X             = a n x m (number trials x number fMRI seed regions) matrix
% Y             = a 3D (trials x channels x time) or 4D (trials x channels x frequency x time) matrix
%
% OUTPUTS:
% betas         = 3D (regressor/ channel/ time) or 4D (regressor/ channel/
%               frequency/ time) matrix with regression weights of given
%               of given fMRI/ behavioral regressor on EEG data.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Adapted from Tobias Hauser (https://github.com/tuhauser/TAfT).
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/CueLockedAnalyses/TAfT/

fprintf('Subject %03d: Fit GLM across trials\n', job.subID);

%% Retrieve number of BOLD regressors (to-be-retrieved from betas later):

nReg    = size(X,2); % count number of variables in current design matrix

%% Add behavioral regressors:

if isfield(job,'behavReg')
    
    nBehav  = length(job.behavReg); % number behavioral variables
    nReg    = nReg + nBehav; % update number of to-be-extracted regressors
    behav   = taft_preprocess_load_behavior(job); % load behavior
    
    XBehav  = nan(size(X,1),nBehav); % initialize behavioral design matrix
    
    for iBehav = 1:nBehav
        tmp = behav.(char(job.behavReg(iBehav))); % extract relevant behavior variables
        XBehav(:,iBehav) = tmp(job.goodTrlIdx,:); % select good trials after EEG trial rejection
    end
    
    % Demean:
    if job.demeanX
        XBehav = XBehav - nanmean(XBehav,1);
    end
    
    % Add to design matrix:
    X(:,(end+1):(end+nBehav)) = XBehav; % add to design matrix

end
%% Add intercept:

if job.add_intercept
    X(:,end+1) = 1;
    fprintf('Add intercept\n');
end

% Note: regressors might have been demeaned above, but intercept isn't (would be zero if demeaned)

%% Compute pseudo-inverse:

pX = taft_pinv(X); % Compute pseudoinverse of design matrix (fMRI + behavior).

%% Apply to each channel/timing/ frequency:

if strcmp(job.EEGdomain, 'time') % data in time domain

    imax = size(Y,1)*size(Y,2); iIter = 0; prog = 0; % initialize progress bar
    fprintf(1,'Perform %s regression for each channel/time bin: %3d%%', job.regType, prog);
    
    betas   = zeros(nReg,size(Y,1),size(Y,2)); % Initialize empty beta weight matrix (time sample, channel,:)
    
    for iChan = 1:size(Y,1) % loop through channels
    
        for iTime = 1:size(Y,2) % loop through timings
        
            y = squeeze(Y(iChan,iTime,:)); % extract EEG data for this channel/ time bin across trials
                        
            if job.demeanY % demean EEG data across trials (so no intercept necessary)
                y = y - nanmean(y);
            end
            
            if strcmp(job.regType,'pinv')
                binBetas = pX * y; % Y is in format ntrials, nchannels, nsamples
            
            elseif strcmp(job.regType,'fitlm')
                warnings off
                m = fitlm(X,y);
                warnings on
                binBetas = m.Coefficients.Estimate;
            
            elseif strcmp(job.regType,'robustfit')
                warning off
                [binBetas, ~] = robustfit(X,y);
                warning on
                
            else
                error('Unknown regression type')
            end
            
            % Retrieve coefficient per regressor (drop intercept at the end):
            for iReg = 1:nReg
                betas(iReg,iChan,iTime) = binBetas(iReg);
            end
            
            % Update progress bar:
            iIter = iIter + 1; % increment iteration
            prog = 100*(iIter/imax); % update percent completion
            if ismember(ceil(prog),10:10:100)
                fprintf(1,'\b\b\b\b%3.0f%%',prog); % update progress bar
            end
            
        end % end of iFreq
    end % end of iChan
    fprintf('\n');
    
elseif strcmp(job.EEGdomain, 'TF') % data in time-frequency domain
    
    imax = size(Y,1)*size(Y,2)*size(Y,3); iIter = 0; prog = 0; % initialize progress bar
    fprintf(1,'Perform %s regression for each channel/frequency/time bin: %3d%%', job.regType, prog);
    
    betas = zeros(nReg,size(Y,1),size(Y,2),size(Y,3)); % Initialize empty beta weight matrix (time sample, channel,:)
    
    for iChan = 1:size(Y,1) % loop through channels
        
        for iFreq = 1:size(Y,2) % loop through frequencies
        
            for iTime = 1:size(Y,3) % loop through timings
            
                y = squeeze(Y(iChan,iFreq,iTime,:)); % extract EEG data for this channel/ frequency/ time bin across trials
            
                if job.dBconversion % decibel conversion
                    y = 10*log10(y);
                end
            
                if job.demeanY % demean EEG data across trials (so no intercept necessary)
                    y = y - nanmean(y);
                end
                if strcmp(job.regType,'pinv')
                    binBetas = pX * y; % Y is in format ntrials, nchannels, nsamples

                elseif strcmp(job.regType,'fitlm')
                    warnings off
                    m = fitlm(X,y);
                    binBetas = m.Coefficients.Estimate;
                    warnings on
                
                elseif strcmp(job.regType,'robustfit')
                    warning off
                    [binBetas, ~] = robustfit(X,y);
                    warning on
                
                else
                    error('Unknown regression type')
                end
                
                % Retrieve coefficient per regressor (drop intercept at the end):
                for iReg = 1:nReg
                    betas(iReg,iChan,iFreq,iTime) = binBetas(iReg);
                end
                
                % Update progress bar:
                iIter = iIter + 1; % increment iteration
                prog = 100*(iIter/imax); % update percent completion
                if ismember(ceil(prog),10:10:100)
                    fprintf(1,'\b\b\b\b%3.0f%%',prog); % update progress bar
                end
                
            end % end of iFreq
        end % end of iTime
    end % end of iChan

    fprintf('\n'); % end of progress bar
    
else
    error('unknown data format')
end

end % end of function.