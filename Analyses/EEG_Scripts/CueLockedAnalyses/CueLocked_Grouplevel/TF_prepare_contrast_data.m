function [job, data] = TF_prepare_contrast_data(job, data)

% [job, data] = TF_prepare_generic_data(job, data)
%
% Add aggregated TF-domain data sets for selected contrast, channels, 
% frequency bands to job object previously set up via TF_update_job.m
% Preliminary data sets created via TF_prepare_contrast_data.m
% 
% INPUTS:
% job               = cell, created via TF_update_job.m, needs at least
% fields:
%   .nSub           = integer, number of subjects.
%   .nCond          = integer, number of conditions.
%   .validSubs      = numeric vector, subject numbers of valid subjects to
%   be included in analyses.
%   .channels       = vector of strings, selected channels.
%   .freq           = numeric vector, frequencies bins to be included.
%   .contrastType   = string, contrast to be used: 'Congruency', 'Go',
%   'Valence', 'GoValence', 'GoAccuracy', 'GoLeftRight', 'Accuracy',
%   'CongruentAccuracy', or 'IncongruentAccuracy'.
%   .responseSettings = string, response setting 'Go' or 'Hand' to be used.
%   .accSettings    = string, accuracy setting 'correct', 'incorrect', or
%   'bothAcc' to be used.
%   .bin.Num        = numeric, number of bins applied (2 or 3; optional).
% data              = cell, need at least the following fields:
% 	TFall{iSub}.ValxAct{iCond} = aggregated data per subject per condition,
% 	created via EEGfMRIPav_Cuelocked_8_TF_grouplevel_create.m
%   .Cond           = grand average per condition across subjects.
%   .mu             = grand average across both conditions and subjects.
%
% OUTPUTS:
% job               = cell, with the following fields added:
%   .chanIdx        = numeric vector, indices of selected channels.
%   .freqIdx        = numeric vector, indices of selected frequencies.
% data              = cell, with the following fields added:
%   .SubCondTime    = per subject per condition, over time, averaged over
%   frequencies and channels.
%   .TF1 and TF2    = averaged over conditions, for permutation tests.
%   .mat1 and mat2  = averaged over channels/ frequencies/ conditions, for
%   two-line plots.
%   .matMean        = mean of mat1 and mat2 (over all conditions).
%   .matGrandMean   = matMean averaged over subjects.
%   .SubTime        = SubCondTime averaged over conditions.
%   .GrandTime      = SubTime averaged over subjects.
%   .topo2plot      = averaged over conditions/ subjects, for topoplots.
%   .TF2plot        = averaged over channels/ conditions/ subjects, for
%   time-frequency plots.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel/

fprintf('Prepare data for contrast %s, channels %s, frequencies %d-%d\n',job.contrastType,strjoin(job.channels,'/'),job.freq(1),job.freq(2))

%% Retrieve indices of selected channels and frequencies:

% Channel indices:
job.chanIdx = find(ismember(data.mu.label, job.channels)); % determine indices of channels per subject

% Frequency indices:
job.freqIdx = dsearchn(data.mu.freq',job.freq'); % indices pf lowest/highest frequency in that range
job.freqIdx = job.freqIdx(1):job.freqIdx(2); % include all indices in between

%% Keep subjects and conditions, average over selected frequencies and channels:

% a) Average within subject and within subject/condition over frequencies and channels:
data.SubCondTime = zeros(job.nSub,job.nCond,length(data.mu.time));

for iSub = job.validSubs % 1:parAll.nSub
    for iCondi = 1:job.nCond

        fprintf('Subject %d, Condition %d \n', iSub, iCondi);
        subChanIdx = ismember(data.TFall{iSub}.ValxAct{iCondi}.label, job.channels); % determine indices of channels per subject
        data.SubCondTime(iSub,iCondi,:) = nanmean(nanmean(data.TFall{iSub}.ValxAct{iCondi}.powspctrm(subChanIdx,job.freqIdx,:),2)); % average first over frequencies, then channels

    end
end

%% Set and initialize data objects 

% Initialize objects:

% a) Permutation test:
data.TF1 = cell(job.nSub,1);
data.TF2 = cell(job.nSub,1);

% b) 2-line plot:
data.mat1 = nan(1,length(data.mu.time)); 
data.mat2 = nan(1,length(data.mu.time)); 

% c) Topoplot:
data.topo2plot = data.Cond{1}; 

%% Compute objects that require looping over subjects (i.e. permutation test):

for iSub = 1:job.nSub
    
    if isfield(job,'bin')

        if strcmp(job.bin.Type,'RT')
            if job.bin.Num == 2
                data.TF1{iSub} = data.TFall{iSub}.ValxAct{1}; % initialize low.
                data.TF1{iSub}.powspctrm = real(data.TFall{iSub}.ValxAct{1}.powspctrm + data.TFall{iSub}.ValxAct{3}.powspctrm)./2; % mean of low
                data.TF2{iSub} = data.TFall{iSub}.ValxAct{2}; % initialize high.
                data.TF2{iSub}.powspctrm = real(data.TFall{iSub}.ValxAct{2}.powspctrm + data.TFall{iSub}.ValxAct{4}.powspctrm)./2; % mean of high

            elseif job.bin.Num == 3
                data.TF1{iSub} = data.TFall{iSub}.ValxAct{1}; % initialize low.
                data.TF1{iSub}.powspctrm = real(data.TFall{iSub}.ValxAct{1}.powspctrm + data.TFall{iSub}.ValxAct{4}.powspctrm)./2; % mean of low
                data.TF2{iSub} = data.TFall{iSub}.ValxAct{2}; % initialize medium.
                data.TF2{iSub}.powspctrm = real(data.TFall{iSub}.ValxAct{2}.powspctrm + data.TFall{iSub}.ValxAct{5}.powspctrm)./2; % mean of medium
                data.TF3{iSub} = data.TFall{iSub}.ValxAct{3}; % initialize high.
                data.TF3{iSub}.powspctrm = real(data.TFall{iSub}.ValxAct{3}.powspctrm + data.TFall{iSub}.ValxAct{6}.powspctrm)./2; % mean of high

            else
                error('Bin number %d not implemented yet',job.bin.Num);
            end
        
        else % if BOLD
            
            data.TF1{iSub} = data.TFall{iSub}.ValxAct{1}; % initialize low.
            data.TF2{iSub} = data.TFall{iSub}.ValxAct{2}; % initialize medium.
            
        end
                        
    elseif strcmp(job.contrastType,'Congruency')

        if strcmp(job.responseSettings,'Go')
        
            if strcmp(job.accSettings,'correct') || strcmp(job.accSettings,'incorrect') 
                data.TF1{iSub} = data.TFall{iSub}.ValxAct{1}; % initialize congruent.
                data.TF1{iSub}.powspctrm = real(data.TFall{iSub}.ValxAct{1}.powspctrm + data.TFall{iSub}.ValxAct{4}.powspctrm)./2; % mean of Go2Win and NoGo2Avoid
                data.TF2{iSub} = data.TFall{iSub}.ValxAct{2}; % initialize incongruent.
                data.TF2{iSub}.powspctrm = real(data.TFall{iSub}.ValxAct{2}.powspctrm + data.TFall{iSub}.ValxAct{3}.powspctrm)./2; % mean of Go2Avoid and NoGo2Win
            
            elseif strcmp(job.accSettings,'bothAcc') 
                data.TF1{iSub} = data.TFall{iSub}.ValxAct{1}; % initialize congruent.
                data.TF1{iSub}.powspctrm = real(data.TFall{iSub}.ValxAct{1}.powspctrm + data.TFall{iSub}.ValxAct{4}.powspctrm  + data.TFall{iSub}.ValxAct{5}.powspctrm + data.TFall{iSub}.ValxAct{8}.powspctrm)./4; % mean of Go2Win and NoGo2Avoid
                data.TF2{iSub} = data.TFall{iSub}.ValxAct{2}; % initialize incongruent.
                data.TF2{iSub}.powspctrm = real(data.TFall{iSub}.ValxAct{2}.powspctrm + data.TFall{iSub}.ValxAct{3}.powspctrm + data.TFall{iSub}.ValxAct{6}.powspctrm + data.TFall{iSub}.ValxAct{7}.powspctrm)./4; % mean of Go2Avoid and NoGo2Win
            
            else
                error('Unknown accuracy setting')
            end
            
        else
            error('Unknown response setting')
        end
        
    elseif strcmp(job.contrastType,'Go')
    
        if strcmp(job.responseSettings,'Go')
            data.TF1{iSub} = data.TFall{iSub}.ValxAct{1}; % initialize Go.
            data.TF1{iSub}.powspctrm = real(data.TFall{iSub}.ValxAct{1}.powspctrm + data.TFall{iSub}.ValxAct{2}.powspctrm)./2; % mean of Go2Win and NoGo2Avoid
            data.TF2{iSub} = data.TFall{iSub}.ValxAct{2}; % initialize NoGo.
            data.TF2{iSub}.powspctrm = real(data.TFall{iSub}.ValxAct{3}.powspctrm + data.TFall{iSub}.ValxAct{4}.powspctrm)./2; % mean of Go2Avoid and NoGo2Win

        elseif strcmp(job.responseSettings,'Hand')
            data.TF1{iSub} = data.TFall{iSub}.ValxAct{1}; % initialize Go.
            data.TF1{iSub}.powspctrm = real(data.TFall{iSub}.ValxAct{1}.powspctrm + data.TFall{iSub}.ValxAct{2}.powspctrm + data.TFall{iSub}.ValxAct{3}.powspctrm + data.TFall{iSub}.ValxAct{4}.powspctrm)./4;
            data.TF2{iSub} = data.TFall{iSub}.ValxAct{2}; % initialize NoGo.
            data.TF2{iSub}.powspctrm = real(data.TFall{iSub}.ValxAct{5}.powspctrm + data.TFall{iSub}.ValxAct{6}.powspctrm)./2;
        
        else
            error('Unknown response setting')
        end
        
    elseif strcmp(job.contrastType,'Valence')
        
        if strcmp(job.responseSettings,'Go')
            data.TF1{iSub} = data.TFall{iSub}.ValxAct{1}; % initialize Win.
            data.TF1{iSub}.powspctrm = real(data.TFall{iSub}.ValxAct{1}.powspctrm + data.TFall{iSub}.ValxAct{3}.powspctrm)./2; % mean of Go2Win and NoGo2Avoid
            data.TF2{iSub} = data.TFall{iSub}.ValxAct{2}; % initialize Avoid.
            data.TF2{iSub}.powspctrm = real(data.TFall{iSub}.ValxAct{2}.powspctrm + data.TFall{iSub}.ValxAct{4}.powspctrm)./2; % mean of Go2Avoid and NoGo2Win
        
        elseif strcmp(job.responseSettings,'Hand')
            data.TF1{iSub} = data.TFall{iSub}.ValxAct{1}; % initialize Win.
            data.TF1{iSub}.powspctrm = real((data.TFall{iSub}.ValxAct{1}.powspctrm + data.TFall{iSub}.ValxAct{3}.powspctrm)./2 + data.TFall{iSub}.ValxAct{5}.powspctrm)./2;
            data.TF2{iSub} = data.TFall{iSub}.ValxAct{3}; % initialize Avoid.
            data.TF2{iSub}.powspctrm = real((data.TFall{iSub}.ValxAct{2}.powspctrm + data.TFall{iSub}.ValxAct{4}.powspctrm)./2 + data.TFall{iSub}.ValxAct{6}.powspctrm)./2;
        
        else
            error('Unknown response setting')
        end
        
    elseif strcmp(job.contrastType,'GoValence')
        
        if strcmp(job.responseSettings,'Go')
            data.TF1{iSub} = data.TFall{iSub}.ValxAct{1}; % initialize Go2Win.
            data.TF1{iSub}.powspctrm = real(data.TFall{iSub}.ValxAct{1}.powspctrm); % only Go2Win
            data.TF2{iSub} = data.TFall{iSub}.ValxAct{2}; % initialize Go2Avoid.
            data.TF2{iSub}.powspctrm = real(data.TFall{iSub}.ValxAct{2}.powspctrm); % only Go2Avoid
        
        elseif strcmp(job.responseSettings,'Hand')
            data.TF1{iSub} = data.TFall{iSub}.ValxAct{1}; % initialize Go2Win.
            data.TF1{iSub}.powspctrm = real(data.TFall{iSub}.ValxAct{1}.powspctrm + data.TFall{iSub}.ValxAct{3}.powspctrm)./2; % only Go2Win
            data.TF2{iSub} = data.TFall{iSub}.ValxAct{3}; % initialize Go2Avoid.
            data.TF2{iSub}.powspctrm = real(data.TFall{iSub}.ValxAct{2}.powspctrm + data.TFall{iSub}.ValxAct{4}.powspctrm)./2; % only Go2Avoid
        
        else
            error('Unknown response setting')
        end
        
    elseif strcmp(job.contrastType,'GoAccuracy')
        
        if strcmp(job.responseSettings,'Go')
            data.TF1{iSub} = data.TFall{iSub}.ValxAct{1}; % initialize GoCorrect.
            data.TF1{iSub}.powspctrm = real(data.TFall{iSub}.ValxAct{1}.powspctrm + data.TFall{iSub}.ValxAct{2}.powspctrm)./2; % GoCorrect
            data.TF2{iSub} = data.TFall{iSub}.ValxAct{5}; % initialize GoIncorrect.
            data.TF2{iSub}.powspctrm = real(data.TFall{iSub}.ValxAct{5}.powspctrm + data.TFall{iSub}.ValxAct{6}.powspctrm)./2; % only GoIncorrect
        
        elseif strcmp(job.responseSettings,'Hand')
            data.TF1{iSub} = data.TFall{iSub}.ValxAct{1}; % initialize GoCorrect.
            data.TF1{iSub}.powspctrm = real(data.TFall{iSub}.ValxAct{1}.powspctrm + data.TFall{iSub}.ValxAct{2}.powspctrm + data.TFall{iSub}.ValxAct{3}.powspctrm + data.TFall{iSub}.ValxAct{4}.powspctrm)./4; % mean of LeftGo2Win and LeftGo2Avoid
            data.TF2{iSub} = data.TFall{iSub}.ValxAct{1}; % initialize GoIncorrect.
            data.TF2{iSub}.powspctrm = real(data.TFall{iSub}.ValxAct{7}.powspctrm + data.TFall{iSub}.ValxAct{8}.powspctrm + data.TFall{iSub}.ValxAct{9}.powspctrm + data.TFall{iSub}.ValxAct{10}.powspctrm)./4; % mean of RightGo2Win and RightGo2Avoid
        
        else
            error('Unknown response setting')
        end
        
    elseif strcmp(job.contrastType,'GoLeftRight') && strcmp(job.responseSettings,'Hand')
    
        if strcmp(job.accSettings,'correct') || strcmp(job.accSettings,'incorrect')
            data.TF1{iSub} = data.TFall{iSub}.ValxAct{1}; % initialize Left.
            data.TF1{iSub}.powspctrm = real(data.TFall{iSub}.ValxAct{1}.powspctrm + data.TFall{iSub}.ValxAct{2}.powspctrm)./2; % mean of LeftGo2Win and LeftGo2Avoid
            data.TF2{iSub} = data.TFall{iSub}.ValxAct{3}; % initialize Right.
            data.TF2{iSub}.powspctrm = real(data.TFall{iSub}.ValxAct{3}.powspctrm + data.TFall{iSub}.ValxAct{4}.powspctrm)./2; % mean of RightGo2Win and RightGo2Avoid
        
        elseif strcmp(job.accSettings,'bothAcc')
            data.TF1{iSub} = data.TFall{iSub}.ValxAct{1}; % initialize Left.
            data.TF1{iSub}.powspctrm = real(data.TFall{iSub}.ValxAct{1}.powspctrm + data.TFall{iSub}.ValxAct{2}.powspctrm + data.TFall{iSub}.ValxAct{7}.powspctrm + data.TFall{iSub}.ValxAct{8}.powspctrm)./4; % mean of LeftGo2Win and LeftGo2Avoid
            data.TF2{iSub} = data.TFall{iSub}.ValxAct{1}; % initialize Right.
            data.TF2{iSub}.powspctrm = real(data.TFall{iSub}.ValxAct{3}.powspctrm + data.TFall{iSub}.ValxAct{4}.powspctrm + data.TFall{iSub}.ValxAct{9}.powspctrm + data.TFall{iSub}.ValxAct{10}.powspctrm)./4; % mean of RightGo2Win and RightGo2Avoid
        
        else
            error('Unknown accuracy setting')
        end
               
    elseif strcmp(job.contrastType,'Accuracy') && strcmp(job.accSettings,'bothAcc')
        data.TF1{iSub} = data.TFall{iSub}.ValxAct{1}; % initialize incongruent.
        data.TF1{iSub}.powspctrm = real(data.TFall{iSub}.ValxAct{1}.powspctrm + data.TFall{iSub}.ValxAct{2}.powspctrm + data.TFall{iSub}.ValxAct{3}.powspctrm + data.TFall{iSub}.ValxAct{4}.powspctrm)./2; % mean of correct Go2Avoid and NoGo2Win
        data.TF2{iSub} = data.TFall{iSub}.ValxAct{1}; % initialize incongruent.
        data.TF2{iSub}.powspctrm = real(data.TFall{iSub}.ValxAct{5}.powspctrm + data.TFall{iSub}.ValxAct{6}.powspctrm + data.TFall{iSub}.ValxAct{7}.powspctrm + data.TFall{iSub}.ValxAct{8}.powspctrm)./2; % mean of correct Go2Avoid and NoGo2Win

    elseif strcmp(job.contrastType,'CongruentAccuracy') && strcmp(job.accSettings,'bothAcc')
        data.TF1{iSub} = data.TFall{iSub}.ValxAct{1}; % initialize incongruent.
        data.TF1{iSub}.powspctrm = real(data.TFall{iSub}.ValxAct{1}.powspctrm + data.TFall{iSub}.ValxAct{4}.powspctrm)./2; % mean of correct Go2Avoid and NoGo2Win
        data.TF2{iSub} = data.TFall{iSub}.ValxAct{1}; % initialize incongruent.
        data.TF2{iSub}.powspctrm = real(data.TFall{iSub}.ValxAct{5}.powspctrm + data.TFall{iSub}.ValxAct{8}.powspctrm)./2; % mean of incorrect Go2Avoid and NoGo2Win
    
    elseif strcmp(job.contrastType,'IncongruentAccuracy') && strcmp(job.accSettings,'bothAcc') 
        data.TF1{iSub} = data.TFall{iSub}.ValxAct{1}; % initialize incongruent.
        data.TF1{iSub}.powspctrm = real(data.TFall{iSub}.ValxAct{2}.powspctrm + data.TFall{iSub}.ValxAct{3}.powspctrm)./2; % mean of correct Go2Avoid and NoGo2Win
        data.TF2{iSub} = data.TFall{iSub}.ValxAct{1}; % initialize incongruent.
        data.TF2{iSub}.powspctrm = real(data.TFall{iSub}.ValxAct{6}.powspctrm + data.TFall{iSub}.ValxAct{7}.powspctrm)./2; % mean of incorrect Go2Avoid and NoGo2Win
       
    else
        error('Unknown contrast setting')
    end
end

%% COMPUTE OBJECTS WITHOUT LOOPING over subjects:

if isfield(job,'bin')

    if strcmp(job.bin.Type,'RT')
        
        if job.bin.Num == 2
            % a) 2-line plot (averaged over subjects):
            data.mat1 = squeeze((data.SubCondTime(:,1,:) + data.SubCondTime(:,3,:))/2); % low
            data.mat2 = squeeze((data.SubCondTime(:,2,:) + data.SubCondTime(:,4,:))/2); % high
            % b) Topoplot:
            data.topo2plot.powspctrm = (data.Cond{2}.powspctrm + data.Cond{4}.powspctrm)./2 ...
            - (data.Cond{1}.powspctrm + data.Cond{3}.powspctrm)./2; % high minus low
            % c) TF plot:
            if length(job.chanIdx) == 1 % if only 1 data.channel: don't average over channels
                data.TF2plot = squeeze(data.Cond{2}.powspctrm(job.chanIdx,:,:) + data.Cond{4}.powspctrm(job.chanIdx,:,:) ...
                    - data.Cond{1}.powspctrm(job.chanIdx,:,:) - data.Cond{3}.powspctrm(job.chanIdx,:,:))./2;
            else
                data.TF2plot = squeeze(nanmean(data.Cond{2}.powspctrm(job.chanIdx,:,:) + data.Cond{4}.powspctrm(job.chanIdx,:,:) ...
                    - data.Cond{1}.powspctrm(job.chanIdx,:,:) - data.Cond{3}.powspctrm(job.chanIdx,:,:)))./2;
            end
    
        elseif job.bin.Num == 3
            % a) 2-line plot (averaged over subjects):
            data.mat1 = squeeze((data.SubCondTime(:,1,:) + data.SubCondTime(:,4,:))/2); % low
%            data.mat2 = squeeze((data.SubCondTime(:,2,:) + data.SubCondTime(:,5,:))/2); % medium
            data.mat2 = squeeze((data.SubCondTime(:,3,:) + data.SubCondTime(:,6,:))/2); % high
            % b) Topoplot:
            data.topo2plot.powspctrm = (data.Cond{3}.powspctrm + data.Cond{6}.powspctrm)./2 ...
            - (data.Cond{1}.powspctrm + data.Cond{4}.powspctrm)./2; % high minus low
            % c) TF plot:
            if length(job.chanIdx) == 1 % if only 1 data.channel: don't average over channels
                data.TF2plot = squeeze(data.Cond{3}.powspctrm(job.chanIdx,:,:) + data.Cond{6}.powspctrm(job.chanIdx,:,:) ...
                    - data.Cond{1}.powspctrm(job.chanIdx,:,:) - data.Cond{3}.powspctrm(job.chanIdx,:,:))./2;
            else
                data.TF2plot = squeeze(nanmean(data.Cond{3}.powspctrm(job.chanIdx,:,:) + data.Cond{6}.powspctrm(job.chanIdx,:,:) ...
                    - data.Cond{1}.powspctrm(job.chanIdx,:,:) - data.Cond{3}.powspctrm(job.chanIdx,:,:)))./2;
            end
    
        else
            error('Bin number %d not implemented yet',job.bin.Num);
        end
        
    else % if BOLD
        
        % a) 2-line plot (averaged over subjects):        
        data.mat1 = squeeze(data.SubCondTime(:,1,:)); % high
        data.mat2 = squeeze(data.SubCondTime(:,2,:)); % low
        % b) Topoplot:
        data.topo2plot.powspctrm = data.Cond{1}.powspctrm...
            - data.Cond{2}.powspctrm; % high minus low
        % c) TF plot:
        if length(job.chanIdx) == 1 % if only 1 data.channel: don't average over channels
            data.TF2plot = squeeze(data.Cond{2}.powspctrm(job.chanIdx,:,:) ...
                - data.Cond{1}.powspctrm(job.chanIdx,:,:))./2;
        else
            data.TF2plot = squeeze(nanmean(data.Cond{2}.powspctrm(job.chanIdx,:,:) ...
                - data.Cond{1}.powspctrm(job.chanIdx,:,:)))./2;
        end
        
    end
    
elseif strcmp(job.contrastType,'Congruency')
    
    if strcmp(job.responseSettings,'Go') 

        if strcmp(job.accSettings,'correct') || strcmp(job.accSettings,'incorrect')
            % a) 2-line plot (averaged over subjects):
            data.mat1 = squeeze((data.SubCondTime(:,2,:) + data.SubCondTime(:,3,:))/2); % incongruent
            data.mat2 = squeeze((data.SubCondTime(:,1,:) + data.SubCondTime(:,4,:))/2); % congruent
            % b) Topoplot:
            data.topo2plot.powspctrm = (data.Cond{2}.powspctrm+data.Cond{3}.powspctrm)./2 ...
            - (data.Cond{1}.powspctrm + data.Cond{4}.powspctrm)./2; % incongruent - congruent
            % c) TF plot:
            if length(job.chanIdx) == 1 % if only 1 data.channel: don't average over channels
                data.TF2plot = squeeze(data.Cond{2}.powspctrm(job.chanIdx,:,:) + data.Cond{3}.powspctrm(job.chanIdx,:,:) ...
                    - data.Cond{1}.powspctrm(job.chanIdx,:,:) - data.Cond{4}.powspctrm(job.chanIdx,:,:))./2;
            else
                data.TF2plot = squeeze(nanmean(data.Cond{2}.powspctrm(job.chanIdx,:,:) + data.Cond{3}.powspctrm(job.chanIdx,:,:) ...
                    - data.Cond{1}.powspctrm(job.chanIdx,:,:) - data.Cond{4}.powspctrm(job.chanIdx,:,:)))./2;
            end
            
        elseif strcmp(job.accSettings,'bothAcc')
            % a) 2-line plot (averaged over subjects):
            data.mat1 = squeeze((data.SubCondTime(:,2,:) + data.SubCondTime(:,3,:) + data.SubCondTime(:,6,:) + data.SubCondTime(:,7,:))/4); % incongruent
            data.mat2 = squeeze((data.SubCondTime(:,1,:) + data.SubCondTime(:,4,:) + data.SubCondTime(:,5,:) + data.SubCondTime(:,8,:))/4); % congruent
            % b) Topoplot:
            data.topo2plot.powspctrm = (data.Cond{2}.powspctrm + data.Cond{3}.powspctrm + data.Cond{6}.powspctrm + data.Cond{7}.powspctrm)./4 ...
            - (data.Cond{1}.powspctrm + data.Cond{4}.powspctrm + data.Cond{5}.powspctrm + data.Cond{8}.powspctrm)/4;
            % c) TFplot:             
            if length(job.chanIdx) == 1 % if only 1 job.channel: don't average over channels
                data.TF2plot = squeeze(data.Cond{2}.powspctrm(job.chanIdx,:,:) + data.Cond{3}.powspctrm(job.chanIdx,:,:) + data.Cond{6}.powspctrm(job.chanIdx,:,:) + data.Cond{7}.powspctrm(job.chanIdx,:,:) ...
                    - data.Cond{1}.powspctrm(job.chanIdx,:,:) - data.Cond{4}.powspctrm(job.chanIdx,:,:) - data.Cond{5}.powspctrm(job.chanIdx,:,:) - data.Cond{8}.powspctrm(job.chanIdx,:,:))./2;
            else
                data.TF2plot = squeeze(nanmean(data.Cond{2}.powspctrm(job.chanIdx,:,:) + data.Cond{3}.powspctrm(job.chanIdx,:,:) + data.Cond{6}.powspctrm(job.chanIdx,:,:) + data.Cond{7}.powspctrm(job.chanIdx,:,:) ...
                    - data.Cond{1}.powspctrm(job.chanIdx,:,:) - data.Cond{4}.powspctrm(job.chanIdx,:,:) - data.Cond{5}.powspctrm(job.chanIdx,:,:) - data.Cond{8}.powspctrm(job.chanIdx,:,:)))./2;
            end
            
        else
            error('Invalid accuracy setting')
        end
        
    elseif strcmp(job.responseSettings,'Hand')
    
        if strcmp(job.accSettings,'correct') || strcmp(job.accSettings,'incorrect') 
            % a) 2- line plot:
            data.mat1 = squeeze(((data.SubCondTime(:,2,:) + data.SubCondTime(:,4,:))/2 + data.SubCondTime(:,5,:))/2); % incongruent
            data.mat2 = squeeze(((data.SubCondTime(:,1,:) + data.SubCondTime(:,3,:))/2 + data.SubCondTime(:,6,:))/2); % congruent 
            % b) Topoplot:
            data.topo2plot.powspctrm = ((data.Cond{2}.powspctrm + data.Cond{4}.powspctrm)./2 + data.Cond{5}.powspctrm)./2 ...
            - ((data.Cond{1}.powspctrm + data.Cond{3}.powspctrm)./2 + data.Cond{6}.powspctrm)./2;
            % c) TF2plot:
            if length(job.chanIdx) == 1 % if only 1 data.channel: don't average over channels
                data.TF2plot = squeeze((data.Cond{2}.powspctrm(job.chanIdx,:,:) + data.Cond{4}.powspctrm(job.chanIdx,:,:))./2 + data.Cond{5}.powspctrm(job.chanIdx,:,:) ...
                    - (data.Cond{1}.powspctrm(job.chanIdx,:,:) + data.Cond{3}.powspctrm(job.chanIdx,:,:))./2 - data.Cond{6}.powspctrm(job.chanIdx,:,:))./2;
            else
                data.TF2plot = squeeze(nanmean((data.Cond{2}.powspctrm(job.chanIdx,:,:) + data.Cond{4}.powspctrm(job.chanIdx,:,:))./2 + data.Cond{5}.powspctrm(job.chanIdx,:,:) ...
                    - (data.Cond{1}.powspctrm(job.chanIdx,:,:) + data.Cond{3}.powspctrm(job.chanIdx,:,:))./2 - data.Cond{6}.powspctrm(job.chanIdx,:,:)))./2;        
            end
            
        elseif strcmp(job.accSettings,'bothAcc')
            % a) 2- line plot:
            data.mat1 = squeeze(((data.SubCondTime(:,2,:) + data.SubCondTime(:,4,:))/2 + data.SubCondTime(:,5,:) + (data.SubCondTime(:,8,:) + data.SubCondTime(:,10,:))/2 + data.SubCondTime(:,11,:))/4); % incongruent
            data.mat2 = squeeze(((data.SubCondTime(:,1,:) + data.SubCondTime(:,3,:))/2 + data.SubCondTime(:,6,:) + (data.SubCondTime(:,7,:) + data.SubCondTime(:,9,:))/2 + data.SubCondTime(:,12,:))/4); % congruent 
            % b) Topoplot:
            data.topo2plot.powspctrm = ((data.Cond{2}.powspctrm + data.Cond{4}.powspctrm)./2 + data.Cond{5}.powspctrm + (data.Cond{8}.powspctrm + data.Cond{10}.powspctrm)./2 + data.Cond{11}.powspctrm)./4 ...
                - ((data.Cond{1}.powspctrm + data.Cond{3}.powspctrm)./2 + data.Cond{6}.powspctrm + (data.Cond{7}.powspctrm + data.Cond{9}.powspctrm)./2 + data.Cond{12}.powspctrm)./4;
            % c) TFplot:
            if length(job.chanIdx) == 1 % if only 1 data.channel: don't average over channels
                data.TF2plot = squeeze((data.Cond{2}.powspctrm(job.chanIdx,:,:) + data.Cond{4}.powspctrm(job.chanIdx,:,:))./2 + data.Cond{5}.powspctrm(job.chanIdx,:,:) + (data.Cond{8}.powspctrm(job.chanIdx,:,:) + data.Cond{10}.powspctrm(job.chanIdx,:,:))./2 + data.Cond{11}.powspctrm(job.chanIdx,:,:) ...
                    - (data.Cond{1}.powspctrm(job.chanIdx,:,:) + data.Cond{3}.powspctrm(job.chanIdx,:,:))./2 - data.Cond{6}.powspctrm(job.chanIdx,:,:) - (data.Cond{7}.powspctrm(job.chanIdx,:,:) + data.Cond{9}.powspctrm(job.chanIdx,:,:))./2 - data.Cond{12}.powspctrm(job.chanIdx,:,:))./2;
            else
                data.TF2plot = squeeze(nanmean((data.Cond{2}.powspctrm(job.chanIdx,:,:) + data.Cond{4}.powspctrm(job.chanIdx,:,:))./2 + data.Cond{5}.powspctrm(job.chanIdx,:,:) + (data.Cond{8}.powspctrm(job.chanIdx,:,:) + data.Cond{10}.powspctrm(job.chanIdx,:,:))./2 + data.Cond{11}.powspctrm(job.chanIdx,:,:) ...
                    - (data.Cond{1}.powspctrm(job.chanIdx,:,:) + data.Cond{3}.powspctrm(job.chanIdx,:,:))./2 - data.Cond{6}.powspctrm(job.chanIdx,:,:) - (data.Cond{7}.powspctrm(job.chanIdx,:,:) + data.Cond{9}.powspctrm(job.chanIdx,:,:))./2 - data.Cond{12}.powspctrm(job.chanIdx,:,:)))./2;
            end
            
        else
            error('Unknown accuracy setting')
        end
        
    else
        error('Inknown response setting')
    end
    
elseif strcmp(job.contrastType,'Go')

    if strcmp(job.responseSettings,'Go')
        % a) 2-line plot:
        data.mat1 = squeeze((data.SubCondTime(:,1,:) + data.SubCondTime(:,2,:))/2);
        data.mat2 = squeeze((data.SubCondTime(:,3,:) + data.SubCondTime(:,4,:))/2);
        % b) Topoplot:
        data.topo2plot.powspctrm = (data.Cond{1}.powspctrm + data.Cond{2}.powspctrm)./2 ...
        - (data.Cond{3}.powspctrm + data.Cond{4}.powspctrm)/2;
        % c) TFplot:
        if length(job.chanIdx) == 1 % if only 1 data.channel: don't average over channels
            data.TF2plot = squeeze(data.Cond{1}.powspctrm(job.chanIdx,:,:) + data.Cond{2}.powspctrm(job.chanIdx,:,:) ...
                - data.Cond{3}.powspctrm(job.chanIdx,:,:) - data.Cond{4}.powspctrm(job.chanIdx,:,:))./2;
        else
            data.TF2plot = squeeze(nanmean(data.Cond{1}.powspctrm(job.chanIdx,:,:) + data.Cond{2}.powspctrm(job.chanIdx,:,:) ...
                - data.Cond{3}.powspctrm(job.chanIdx,:,:) - data.Cond{4}.powspctrm(job.chanIdx,:,:)))./2;
        end
        
    elseif strcmp(job.responseSettings,'Hand')
        % a) 2-D line:
        data.mat1 = squeeze((data.SubCondTime(:,1,:) + data.SubCondTime(:,2,:) + data.SubCondTime(:,3,:) + data.SubCondTime(:,4,:))/4); % Go 
        data.mat2 = squeeze((data.SubCondTime(:,5,:) + data.SubCondTime(:,6,:))/2); % NoGo
        % b) Topoplot:
        data.topo2plot.powspctrm = (data.Cond{1}.powspctrm + data.Cond{2}.powspctrm + data.Cond{3}.powspctrm + data.Cond{4}.powspctrm)./4 ...
        - (data.Cond{5}.powspctrm + data.Cond{6}.powspctrm)./2;
        % c) TFplot:
        if length(job.chanIdx) == 1 % if only 1 data.channel: don't average over channels
            data.TF2plot = squeeze((data.Cond{1}.powspctrm(job.chanIdx,:,:) + data.Cond{2}.powspctrm(job.chanIdx,:,:) + data.Cond{3}.powspctrm(job.chanIdx,:,:) + data.Cond{4}.powspctrm(job.chanIdx,:,:))/4 ...
                - (data.Cond{5}.powspctrm(job.chanIdx,:,:) + data.Cond{6}.powspctrm(job.chanIdx,:,:))./2)./2;
        else
            data.TF2plot = squeeze(nanmean((data.Cond{1}.powspctrm(job.chanIdx,:,:) + data.Cond{2}.powspctrm(job.chanIdx,:,:) + data.Cond{3}.powspctrm(job.chanIdx,:,:) + data.Cond{4}.powspctrm(job.chanIdx,:,:))/4 ...
                - (data.Cond{5}.powspctrm(job.chanIdx,:,:) + data.Cond{6}.powspctrm(job.chanIdx,:,:))./2))./2;
        end
        
    else
        error('Unknown response setting')
    end
    
elseif strcmp(job.contrastType,'Valence')
    
    if strcmp(job.responseSettings,'Go')
        % a) 2-line plot
        data.mat1 = squeeze((data.SubCondTime(:,1,:) + data.SubCondTime(:,3,:))/2);
        data.mat2 = squeeze((data.SubCondTime(:,2,:) + data.SubCondTime(:,4,:))/2);
        % b) Topoplot:
        data.topo2plot.powspctrm = (data.Cond{1}.powspctrm+data.Cond{3}.powspctrm)./2 ...
        - (data.Cond{2}.powspctrm + data.Cond{4}.powspctrm)./2;
        % c) TFplot:
        if length(job.chanIdx) == 1 % if only 1 data.channel: don't average over channels
            data.TF2plot = squeeze(data.Cond{1}.powspctrm(job.chanIdx,:,:) + data.Cond{3}.powspctrm(job.chanIdx,:,:) ...
                - data.Cond{2}.powspctrm(job.chanIdx,:,:) - data.Cond{4}.powspctrm(job.chanIdx,:,:))./2;
        else
            data.TF2plot = squeeze(nanmean(data.Cond{1}.powspctrm(job.chanIdx,:,:) + data.Cond{3}.powspctrm(job.chanIdx,:,:) ...
                - data.Cond{2}.powspctrm(job.chanIdx,:,:) - data.Cond{4}.powspctrm(job.chanIdx,:,:)))./2;
        end
        
    elseif strcmp(job.responseSettings,'Hand')
        % a) 2-line plot:
        data.mat1 = squeeze(((data.SubCondTime(:,1,:) + data.SubCondTime(:,3,:))/2 + data.SubCondTime(:,5,:))/2); 
        data.mat2 = squeeze(((data.SubCondTime(:,2,:) + data.SubCondTime(:,4,:))/2 + data.SubCondTime(:,6,:))/2); 
        % b) Topoplot:
        data.topo2plot.powspctrm = ((data.Cond{1}.powspctrm + data.Cond{3}.powspctrm)./2 + data.Cond{5}.powspctrm)./2 ...
        - ((data.Cond{2}.powspctrm + data.Cond{4}.powspctr)./2 + data.Cond{6}.powspctrm)./2;
        % c) TFplot:
        if length(job.chanIdx) == 1 % if only 1 data.channel: don't average over channels
            data.TF2plot = squeeze((data.Cond{1}.powspctrm(job.chanIdx,:,:) + data.Cond{3}.powspctrm(job.chanIdx,:,:))./2 + data.Cond{5}.powspctrm(job.chanIdx,:,:) ...
                - (data.Cond{2}.powspctrm(job.chanIdx,:,:) + data.Cond{4}.powspctrm(job.chanIdx,:,:))./2 - data.Cond{6}.powspctrm(job.chanIdx,:,:))./2;
        else
            data.TF2plot = squeeze(nanmean((data.Cond{1}.powspctrm(job.chanIdx,:,:) + data.Cond{3}.powspctrm(job.chanIdx,:,:))./2 + data.Cond{5}.powspctrm(job.chanIdx,:,:) ...
                - (data.Cond{2}.powspctrm(job.chanIdx,:,:) + data.Cond{4}.powspctrm(job.chanIdx,:,:))./2 - data.Cond{6}.powspctrm(job.chanIdx,:,:)))./2;
        end
        
    else
        error('Unknown response setting')
    end
    
elseif strcmp(job.contrastType,'GoValence')
    
    if strcmp(job.responseSettings,'Go')
        
        if strcmp(job.accSettings,'correct') || strcmp(job.accSettings,'incorrect')
            % b) 2-line plot:
            data.mat1 = squeeze(data.SubCondTime(:,1,:)); % Go2Win
            data.mat2 = squeeze(data.SubCondTime(:,2,:)); % Go2Avoid
            % c) Topoplot:
            data.topo2plot.powspctrm = data.Cond{1}.powspctrm - data.Cond{2}.powspctrm;
            % d) TFplot:
            if length(job.chanIdx) == 1 % if only 1 data.channel: don't average over channels
                data.TF2plot = squeeze(data.Cond{1}.powspctrm(job.chanIdx,:,:) - data.Cond{2}.powspctrm(job.chanIdx,:,:));
            else
                data.TF2plot = squeeze(nanmean(data.Cond{1}.powspctrm(job.chanIdx,:,:) - data.Cond{2}.powspctrm(job.chanIdx,:,:)))/2;            
            end
            
        elseif strcmp(job.accSettings,'bothAcc')
            % b) 2-line plot:
            data.mat1 = squeeze((data.SubCondTime(:,1,:) + data.SubCondTime(:,5,:))/2); % Go2Win
            data.mat2 = squeeze((data.SubCondTime(:,2,:) + data.SubCondTime(:,6,:))/2); % Go2Avoid
            % c) Topoplot:
            data.topo2plot.powspctrm = data.Cond{1}.powspctrm - data.Cond{2}.powspctrm;
            % d) TFplot:
            if length(job.chanIdx) == 1 % if only 1 data.channel: don't average over channels
                data.TF2plot = squeeze(data.Cond{1}.powspctrm(job.chanIdx,:,:) + data.Cond{5}.powspctrm(job.chanIdx,:,:) ...
                - data.Cond{2}.powspctrm(job.chanIdx,:,:) - data.Cond{6}.powspctrm(job.chanIdx,:,:))./2;
            else
                data.TF2plot = squeeze(nanmean(data.Cond{1}.powspctrm(job.chanIdx,:,:) + data.Cond{5}.powspctrm(job.chanIdx,:,:) ...
                - data.Cond{2}.powspctrm(job.chanIdx,:,:) - data.Cond{6}.powspctrm(job.chanIdx,:,:)))./2;
            end
            
        else
            error('unknown accuracy setting')
        end
        
    elseif strcmp(job.responseSettings,'Hand')
        % a) 2-line plot:
        data.mat1 = squeeze((data.SubCondTime(:,1,:) + data.SubCondTime(:,3,:))/2); % Go2Win
        data.mat2 = squeeze((data.SubCondTime(:,2,:) + data.SubCondTime(:,4,:))/2); % Go2Avoid
        % b) Topoplot:
        data.topo2plot.powspctrm = (data.Cond{1}.powspctrm + data.Cond{3}.powspctrm)./2 ...
            - (data.Cond{2}.powspctrm + data.Cond{4}.powspctrm)./2;
        % c) TFplot:
        if length(job.chanIdx) == 1 % if only 1 data.channel: don't average over channels
            data.TF2plot = squeeze(data.Cond{1}.powspctrm(job.chanIdx,:,:)) + nanmean(data.Cond{3}.powspctrm(job.chanIdx,:,:)) ...
            - nanmean(data.Cond{2}.powspctrm(job.chanIdx,:,:)) - nanmean(data.Cond{4}.powspctrm(job.chanIdx,:,:))./2;
        else
            data.TF2plot = squeeze(nanmean(data.Cond{1}.powspctrm(job.chanIdx,:,:)) + nanmean(data.Cond{3}.powspctrm(job.chanIdx,:,:)) ...
            - nanmean(data.Cond{2}.powspctrm(job.chanIdx,:,:)) - nanmean(data.Cond{4}.powspctrm(job.chanIdx,:,:)))./2;
        end
        
    else
        error('Unknown response setting')
    end
    
elseif strcmp(job.contrastType,'GoAccuracy') && strcmp(job.accSettings,'bothAcc')
    
    if strcmp(job.responseSettings,'Go')
        % a) 2-line plot:
        data.mat1 = squeeze((data.SubCondTime(:,1,:) + data.SubCondTime(:,2,:))/2); % Go correct 
        data.mat2 = squeeze((data.SubCondTime(:,5,:) + data.SubCondTime(:,6,:))/2); % Go incorrect
        % b) Topoplot:
        data.topo2plot.powspctrm = (data.Cond{1}.powspctrm + data.Cond{2}.powspctrm)./2 ...
        - (data.Cond{5}.powspctrm + data.Cond{6}.powspctrm)./2;
        % d) TFplot:
        if length(job.chanIdx) == 1 % if only 1 job.channel: don't average over channels
            data.TF2plot = squeeze(data.Cond{1}.powspctrm(job.chanIdx,:,:) + data.Cond{2}.powspctrm(job.chanIdx,:,:) ...
            - data.Cond{5}.powspctrm(job.chanIdx,:,:) - data.Cond{6}.powspctrm(job.chanIdx,:,:))./2;
        else
            data.TF2plot = squeeze(nanmean(data.Cond{1}.powspctrm(job.chanIdx,:,:) + data.Cond{2}.powspctrm(job.chanIdx,:,:) ...
            - data.Cond{5}.powspctrm(job.chanIdx,:,:) - data.Cond{6}.powspctrm(job.chanIdx,:,:)))./2;
        end
        
    elseif strcmp(job.responseSettings,'Hand')
        % Decide whether to collapse left/right before averaging or not:
        % a) 2-line plot:
        data.mat1 = squeeze((data.SubCondTime(:,1,:) + data.SubCondTime(:,2,:) + data.SubCondTime(:,3,:) + data.SubCondTime(:,4,:))/4);
        data.mat2 = squeeze((data.SubCondTime(:,7,:) + data.SubCondTime(:,8,:) + data.SubCondTime(:,9,:) + data.SubCondTime(:,10,:))/4);
        % b) Topoplot:
        data.topo2plot.powspctrm = (data.Cond{1}.powspctrm + data.Cond{2}.powspctrm + data.Cond{3}.powspctrm + data.Cond{4}.powspctrm)./4 ...
        - (data.Cond{7}.powspctrm + data.Cond{8}.powspctrm + data.Cond{9}.powspctrm + data.Cond{10}.powspctrm)./4;
        % d) TFplot:
        if length(job.chanIdx) == 1 % if only 1 job.channel: don't average over channels
            data.TF2plot = squeeze(data.Cond{1}.powspctrm(job.chanIdx,:,:) + data.Cond{2}.powspctrm(job.chanIdx,:,:) + data.Cond{3}.powspctrm(job.chanIdx,:,:) + data.Cond{4}.powspctrm(job.chanIdx,:,:)  ...
            - data.Cond{7}.powspctrm(job.chanIdx,:,:) - data.Cond{8}.powspctrm(job.chanIdx,:,:) - data.Cond{9}.powspctrm(job.chanIdx,:,:) - data.Cond{10}.powspctrm(job.chanIdx,:,:))./2;
        else
            data.TF2plot = squeeze(nanmean(data.Cond{1}.powspctrm(job.chanIdx,:,:) + data.Cond{2}.powspctrm(job.chanIdx,:,:) + data.Cond{3}.powspctrm(job.chanIdx,:,:) + data.Cond{4}.powspctrm(job.chanIdx,:,:) ...
            - data.Cond{7}.powspctrm(job.chanIdx,:,:) - data.Cond{8}.powspctrm(job.chanIdx,:,:) - data.Cond{9}.powspctrm(job.chanIdx,:,:) - data.Cond{10}.powspctrm(job.chanIdx,:,:)))./2;
        end
        
    else
        error('invalid response settings')
    end
    
elseif strcmp(job.contrastType,'GoLeftRight') && strcmp(job.responseSettings,'Hand')
    
    if strcmp(job.accSettings,'correct') || strcmp(job.accSettings,'incorrect')
        % a) 2-line plot:
        data.mat1 = squeeze((data.SubCondTime(:,1,:) + data.SubCondTime(:,2,:))/2); % LeftGo
        data.mat2 = squeeze((data.SubCondTime(:,3,:) + data.SubCondTime(:,4,:))/2); % Right Go
        % b) Topoplot:
        data.topo2plot.powspctrm = (data.Cond{1}.powspctrm + data.Cond{2}.powspctrm)./2 ...
        - (data.Cond{3}.powspctrm + data.Cond{4}.powspctrm)./2; % Left minus Right
        % c) TFplot:
        if length(job.chanIdx) == 1 % if only 1 data.channel: don't average over channels
            data.TF2plot = squeeze(data.Cond{1}.powspctrm(job.chanIdx,:,:) + data.Cond{2}.powspctrm(job.chanIdx,:,:) ...
                - data.Cond{3}.powspctrm(job.chanIdx,:,:) - data.Cond{4}.powspctrm(job.chanIdx,:,:))./2;
        else
            data.TF2plot = squeeze(nanmean(data.Cond{1}.powspctrm(job.chanIdx,:,:) + data.Cond{2}.powspctrm(job.chanIdx,:,:) ...
                - data.Cond{3}.powspctrm(job.chanIdx,:,:) - data.Cond{4}.powspctrm(job.chanIdx,:,:)))./2;
        end
        
    elseif strcmp(job.accSettings,'bothAcc')
        % a) 2-line plot:
        data.mat1 = squeeze((data.SubCondTime(:,1,:) + data.SubCondTime(:,2,:) + data.SubCondTime(:,7,:) + data.SubCondTime(:,8,:))/4); % Goleft
        data.mat2 = squeeze((data.SubCondTime(:,3,:) + data.SubCondTime(:,4,:) + data.SubCondTime(:,9,:) + data.SubCondTime(:,10,:))/4); % GoRight
        % b) Topoplot:
        data.topo2plot.powspctrm = (data.Cond{1}.powspctrm + data.Cond{2}.powspctrm + data.Cond{7}.powspctrm + data.Cond{8}.powspctrm)./4 ...
        - (data.Cond{3}.powspctrm + data.Cond{4}.powspctrm + data.Cond{9}.powspctrm + data.Cond{10}.powspctrm)./4;
        % a) TFplot:
        if length(job.chanIdx) == 1 % if only 1 job.channel: don't average over channels
            data.TF2plot = squeeze(data.Cond{1}.powspctrm(job.chanIdx,:,:) + data.Cond{2}.powspctrm(job.chanIdx,:,:) + data.Cond{7}.powspctrm(job.chanIdx,:,:) + data.Cond{8}.powspctrm(job.chanIdx,:,:) ...
                - data.Cond{3}.powspctrm(job.chanIdx,:,:) - data.Cond{4}.powspctrm(job.chanIdx,:,:) - data.Cond{9}.powspctrm(job.chanIdx,:,:) - data.Cond{10}.powspctrm(job.chanIdx,:,:))./2;
        else
            data.TF2plot = squeeze(nanmean(data.Cond{1}.powspctrm(job.chanIdx,:,:) + data.Cond{2}.powspctrm(job.chanIdx,:,:) + data.Cond{7}.powspctrm(job.chanIdx,:,:) + data.Cond{8}.powspctrm(job.chanIdx,:,:) ...
                - data.Cond{3}.powspctrm(job.chanIdx,:,:) - data.Cond{4}.powspctrm(job.chanIdx,:,:)  - data.Cond{9}.powspctrm(job.chanIdx,:,:) - data.Cond{10}.powspctrm(job.chanIdx,:,:)))./2;
        end
        
    else
        error('Unknown accuracy setting')
    end
    
elseif strcmp(job.contrastType,'Accuracy') && strcmp(job.accSettings,'bothAcc')
    
    if strcmp(job.responseSettings,'Go')
        % a) 2-line plot:
        data.mat1 = squeeze((data.SubCondTime(:,1,:) + data.SubCondTime(:,2,:) + data.SubCondTime(:,3,:) + data.SubCondTime(:,4,:))/4); % correct 
        data.mat2 = squeeze((data.SubCondTime(:,5,:) + data.SubCondTime(:,6,:) + data.SubCondTime(:,7,:) + data.SubCondTime(:,8,:))/4); % incorrect
        % b) Topoplot:
        data.topo2plot.powspctrm = (data.Cond{1}.powspctrm + data.Cond{2}.powspctrm + data.Cond{3}.powspctrm + data.Cond{4}.powspctrm)./4 ...
        - (data.Cond{5}.powspctrm + data.Cond{6}.powspctrm + data.Cond{7}.powspctrm + data.Cond{8}.powspctrm)./4;
        % d) TFplot:
        if length(job.chanIdx) == 1 % if only 1 job.channel: don't average over channels
            data.TF2plot = squeeze(data.Cond{1}.powspctrm(job.chanIdx,:,:) + data.Cond{2}.powspctrm(job.chanIdx,:,:) + data.Cond{3}.powspctrm(job.chanIdx,:,:) + data.Cond{4}.powspctrm(job.chanIdx,:,:) ...
            - data.Cond{5}.powspctrm(job.chanIdx,:,:) - data.Cond{6}.powspctrm(job.chanIdx,:,:) - data.Cond{7}.powspctrm(job.chanIdx,:,:) - data.Cond{8}.powspctrm(job.chanIdx,:,:))./2;
        else
            data.TF2plot = squeeze(nanmean(data.Cond{1}.powspctrm(job.chanIdx,:,:) + data.Cond{2}.powspctrm(job.chanIdx,:,:) + data.Cond{3}.powspctrm(job.chanIdx,:,:) + data.Cond{4}.powspctrm(job.chanIdx,:,:) ...
            - data.Cond{5}.powspctrm(job.chanIdx,:,:) - data.Cond{6}.powspctrm(job.chanIdx,:,:) - data.Cond{7}.powspctrm(job.chanIdx,:,:) - data.Cond{8}.powspctrm(job.chanIdx,:,:)))./2;
        end
        
    elseif strcmp(job.responseSettings,'Hand')
        % Decide whether to collapse left/right before averaging or not:
        % a) 2-line plot:
        data.mat1 = squeeze((data.SubCondTime(:,1,:) + data.SubCondTime(:,2,:) + data.SubCondTime(:,3,:) + data.SubCondTime(:,4,:) + data.SubCondTime(:,5,:) + data.SubCondTime(:,6,:))/6);
        data.mat2 = squeeze((data.SubCondTime(:,7,:) + data.SubCondTime(:,8,:) + data.SubCondTime(:,9,:) + data.SubCondTime(:,10,:) + data.SubCondTime(:,11,:) + data.SubCondTime(:,12,:))/6);
        % b) Topoplot:
        data.topo2plot.powspctrm = (data.Cond{1}.powspctrm + data.Cond{2}.powspctrm + data.Cond{3}.powspctrm + data.Cond{4}.powspctrm + data.Cond{5}.powspctrm + data.Cond{6}.powspctrm)./6 ...
        - (data.Cond{7}.powspctrm + data.Cond{8}.powspctrm + data.Cond{9}.powspctrm + data.Cond{10}.powspctrm + data.Cond{11}.powspctrm + data.Cond{12}.powspctrm)./6;
        % d) TFplot:
        if length(job.chanIdx) == 1 % if only 1 job.channel: don't average over channels
            data.TF2plot = squeeze(data.Cond{1}.powspctrm(job.chanIdx,:,:) + data.Cond{2}.powspctrm(job.chanIdx,:,:) + data.Cond{3}.powspctrm(job.chanIdx,:,:) + data.Cond{4}.powspctrm(job.chanIdx,:,:)  + data.Cond{5}.powspctrm(job.chanIdx,:,:)  + data.Cond{6}.powspctrm(job.chanIdx,:,:) ...
            - data.Cond{7}.powspctrm(job.chanIdx,:,:) - data.Cond{8}.powspctrm(job.chanIdx,:,:) - data.Cond{9}.powspctrm(job.chanIdx,:,:) - data.Cond{10}.powspctrm(job.chanIdx,:,:) - data.Cond{11}.powspctrm(job.chanIdx,:,:) - data.Cond{12}.powspctrm(job.chanIdx,:,:))./2;
        else
            data.TF2plot = squeeze(nanmean(data.Cond{1}.powspctrm(job.chanIdx,:,:) + data.Cond{2}.powspctrm(job.chanIdx,:,:) + data.Cond{3}.powspctrm(job.chanIdx,:,:) + data.Cond{4}.powspctrm(job.chanIdx,:,:)  + data.Cond{5}.powspctrm(job.chanIdx,:,:)  + data.Cond{6}.powspctrm(job.chanIdx,:,:) ...
            - data.Cond{7}.powspctrm(job.chanIdx,:,:) - data.Cond{8}.powspctrm(job.chanIdx,:,:) - data.Cond{9}.powspctrm(job.chanIdx,:,:) - data.Cond{10}.powspctrm(job.chanIdx,:,:) - data.Cond{11}.powspctrm(job.chanIdx,:,:) - data.Cond{12}.powspctrm(job.chanIdx,:,:)))./2;
        end
        
    else
        error('invalid response settings')
    end
    
elseif strcmp(job.contrastType,'CongruentAccuracy') && strcmp(job.accSettings,'bothAcc')
    
    if strcmp(job.responseSettings,'Go')
        % a) 2-line plot:
        data.mat1 = squeeze((data.SubCondTime(:,1,:) + data.SubCondTime(:,4,:))/2); % congruent correct
        data.mat2 = squeeze((data.SubCondTime(:,5,:) + data.SubCondTime(:,8,:))/2); % congruent incorrect
        % b) Topoplot:
        data.topo2plot.powspctrm = (data.Cond{1}.powspctrm+data.Cond{4}.powspctrm)./2 ...
        - (data.Cond{5}.powspctrm + data.Cond{8}.powspctrm)./2;
        % d) TFplot:
        if length(job.chanIdx) == 1 % if only 1 job.channel: don't average over channels
            data.TF2plot = squeeze(data.Cond{1}.powspctrm(job.chanIdx,:,:) + data.Cond{4}.powspctrm(job.chanIdx,:,:) ...
            - data.Cond{5}.powspctrm(job.chanIdx,:,:) - data.Cond{8}.powspctrm(job.chanIdx,:,:))./2;
        else
            data.TF2plot = squeeze(nanmean(data.Cond{1}.powspctrm(job.chanIdx,:,:) + data.Cond{4}.powspctrm(job.chanIdx,:,:) ...
            - data.Cond{5}.powspctrm(job.chanIdx,:,:) - data.Cond{8}.powspctrm(job.chanIdx,:,:)))./2;
        end
        
    elseif strcmp(job.responseSettings,'Hand')
        % a) 2-line plot:
        data.mat1 = squeeze(((data.SubCondTime(:,1,:) + data.SubCondTime(:,3,:))/2 + data.SubCondTime(:,6,:))/2); % congruent correct
        data.mat2 = squeeze(((data.SubCondTime(:,7,:) + data.SubCondTime(:,9,:))/2 + data.SubCondTime(:,12,:))/2); % congruent incorrect
        % b) Topoplot:
        data.topo2plot.powspctrm = ((data.Cond{1}.powspctrm + data.Cond{3}.powspctrm)./2 + data.Cond{6}.powspctrm)./2 ...
            - ((data.Cond{7}.powspctrm + data.Cond{9}.powspctrm)./2 + data.Cond{12}.powspctrm)./2;
        % d) TFplot:
        if length(job.chanIdx) == 1 % if only 1 data.channel: don't average over channels
            data.TF2plot = squeeze((data.Cond{1}.powspctrm(job.chanIdx,:,:)+data.Cond{3}.powspctrm(job.chanIdx,:,:))./2 + data.Cond{6}.powspctrm(job.chanIdx,:,:) ...
                - (data.Cond{7}.powspctrm(job.chanIdx,:,:)+data.Cond{9}.powspctrm(job.chanIdx,:,:))./2 - data.Cond{12}.powspctrm(job.chanIdx,:,:))./2;
        else
            data.TF2plot = squeeze(nanmean((data.Cond{1}.powspctrm(job.chanIdx,:,:) + data.Cond{3}.powspctrm(job.chanIdx,:,:))./2 + data.Cond{6}.powspctrm(job.chanIdx,:,:) ...
                - (data.Cond{7}.powspctrm(job.chanIdx,:,:) + data.Cond{9}.powspctrm(job.chanIdx,:,:))./2 - data.Cond{12}.powspctrm(job.chanIdx,:,:)))./2;
        end
        
    else
        error('unknown response settings')
    end
    
elseif strcmp(job.contrastType,'IncongruentAccuracy') && strcmp(job.accSettings,'bothAcc')
    
    if strcmp(job.responseSettings,'Go')
        % a) 2-line plot:
        data.mat1 = squeeze((data.SubCondTime(:,2,:) + data.SubCondTime(:,3,:))/2); % incongruent correct
        data.mat2 = squeeze((data.SubCondTime(:,6,:) + data.SubCondTime(:,7,:))/2); % incongruent incorrect
        % b) Topoplot:
        data.topo2plot.powspctrm = (data.Cond{2}.powspctrm+data.Cond{3}.powspctrm)./2 ...
        - (data.Cond{6}.powspctrm + data.Cond{7}.powspctrm)./2;
        % c) TFplot:             
        if length(job.chanIdx) == 1 % if only 1 job.channel: don't average over channels
            data.TF2plot = squeeze(data.Cond{2}.powspctrm(job.chanIdx,:,:) + data.Cond{3}.powspctrm(job.chanIdx,:,:) ...
                - data.Cond{6}.powspctrm(job.chanIdx,:,:) - data.Cond{7}.powspctrm(job.chanIdx,:,:))./2;
        else
            data.TF2plot = squeeze(nanmean(data.Cond{2}.powspctrm(job.chanIdx,:,:) + data.Cond{3}.powspctrm(job.chanIdx,:,:) ...
                - data.Cond{6}.powspctrm(job.chanIdx,:,:) - data.Cond{7}.powspctrm(job.chanIdx,:,:)))./2;
        end 
        
    elseif strcmp(job.responseSettings,'Hand')
        % a) 2-line plot:
        data.mat1 = squeeze(((data.SubCondTime(:,2,:) + data.SubCondTime(:,4,:))/2 + data.SubCondTime(:,5,:))/2); % incongruent correct
        data.mat2 = squeeze(((data.SubCondTime(:,8,:) + data.SubCondTime(:,10,:))/2 + data.SubCondTime(:,11,:))/2); % incongruent correct
        % b) Topoplot:
        data.topo2plot.powspctrm = ((data.Cond{2}.powspctrm + data.Cond{4}.powspctrm)./2 + data.Cond{5}.powspctrm)./2 ...
            - ((data.Cond{8}.powspctrm + data.Cond{10}.powspctrm)./2 + data.Cond{11}.powspctrm)./2;
        % d) TFplot:
        if length(job.chanIdx) == 1 % if only 1 data.channel: don't average over channels
            data.TF2plot = squeeze((data.Cond{2}.powspctrm(job.chanIdx,:,:) + data.Cond{4}.powspctrm(job.chanIdx,:,:))/2 + data.Cond{5}.powspctrm(job.chanIdx,:,:) ...
                - (data.Cond{8}.powspctrm(job.chanIdx,:,:) + data.Cond{10}.powspctrm(job.chanIdx,:,:))./2 - data.Cond{11}.powspctrm(job.chanIdx,:,:))./2;
        else
            data.TF2plot = squeeze(nanmean((data.Cond{2}.powspctrm(job.chanIdx,:,:) + data.Cond{4}.powspctrm(job.chanIdx,:,:))/2 + data.Cond{5}.powspctrm(job.chanIdx,:,:) ...
                - (data.Cond{8}.powspctrm(job.chanIdx,:,:) + data.Cond{10}.powspctrm(job.chanIdx,:,:))./2 - data.Cond{11}.powspctrm(job.chanIdx,:,:)))./2;
        end
        
    else
        error('unknown response settings')
    end
    
else 
    error('Unknown contrast type')
end

%% Error bounds for 2-line plots:
% Error bounds based on Cousinea-Morey method:

data.matDiff        = data.mat1 - data.mat2; % difference between conditions for each subject
data.matMean        = (data.mat1 + data.mat2)./2; % mean across conditions for each subject
data.matGrandMean   = nanmean(data.matMean); % grand average across subjects
data.matGrandMean   = nanmean(data.matMean); % for Cousineau method
data.SubTime        = squeeze(nanmean(data.SubCondTime,2)); % average over conditions for Cousineau method
data.GrandTime      = squeeze(nanmean(data.SubTime(job.validSubs,:),1)); % average over subjects for Cousineau method

fprintf('Finished preparing data for contrast %s, channels %s, frequencies %d-%d\n',job.contrastType,strjoin(job.channels,'/'),job.freq(1),job.freq(2))
% END
