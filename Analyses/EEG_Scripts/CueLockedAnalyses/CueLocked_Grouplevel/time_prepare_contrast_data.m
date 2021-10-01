function [job, data] = time_prepare_contrast_data(job, data)

% [job, data] = time_prepare_generic_data(job, data)
%
% Add aggregated time-domain data sets for selected contrast, channels 
% to job object previously set up via TF_update_job.m
% Preliminary data sets created via time_prepare_contrast_data.m
% 
% INPUTS:
% job               = cell, created via time_update_job.m, needs at least
% fields:
%   .nSub           = integer, number of subjects.
%   .nCond          = integer, number of conditions.
%   .validSubs      = numeric vector, subject numbers of valid subjects to
%   .channels       = vector of strings, selected channels.
%   .contrastType   = string, contrast to be used: 'Congruency', 'Go',
%   'Valence', 'GoValence', 'GoAccuracy', 'GoLeftRight', 'Accuracy',
%   'CongruentAccuracy', or 'IncongruentAccuracy'.
%   .responseSettings = string, response setting 'Go' or 'Hand' to be used.
%   .accSettings    = string, accuracy setting 'correct', 'incorrect', or
%   'bothAcc' to be used.
%   .bin.Num        = numeric, number of bins applied (2 or 3; optional).
% data              = cell, need at least the following fields:
%   be included in analyses.
% 	ERPdata{iSub,iCond} = aggregated data per subject per condition,
% 	created via EEGfMRIPav_Cuelocked_8_time_grouplevel_create.m
%   .Cond           = grand average per condition across subjects.
%   .mu             = grand average across both conditions and subjects.
%
% OUTPUTS:
% job               = cell, with the following fields added:
%   .chanIdx        = numeric vector, indices of selected channels.
% data              = cell, with the following fields added:
%   .SubCondTime    = per subject per condition, over time, averaged over
%   channels.
%   .time1 and time2= averaged over conditions, for permutation tests.
%   .mat1 and mat2  = averaged over channels/ conditions, for two-line 
%   plots.
%   .matMean        = mean of mat1 and mat2 (over all conditions).
%   .matGrandMean   = matMean averaged over subjects.
%   .SubTime        = SubCondTime averaged over conditions.
%   .GrandTime      = SubTime averaged over subjects.
%   .topo2plot      = averaged over conditions/ subjects, for topoplots.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel/

fprintf('Prepare data for contrast %s, channels %s\n',job.contrastType,strjoin(job.channels,'/'))

%% Retrieve indices of selected channels:

% Channel indices:
job.chanIdx = find(ismember(data.mu.label, job.channels)); % determine indices of channels per subject

%% Keep subjects and conditions and time, average over selected channels:

% a) Average within subject and within subject/condition over channels:
data.SubCondTime = zeros(job.nSub,job.nCond,length(data.mu.time));

for iSub = job.validSubs % 1:parAll.nSub
    for iCondi = 1:job.nCond
        fprintf('Subject %d, Condition %d: Average over channels %s \n', iSub, iCondi, strjoin(job.channels,'/'));
        subChanIdx = ismember(data.ERPdata{iSub,iCondi}.label, job.channels); % determine indices of channels per subject
        data.SubCondTime(iSub,iCondi,:) = nanmean(data.ERPdata{iSub,iCondi}.avg(subChanIdx,:),1); % average over channels
    end
end

%% Set and initialize data objects 

clear data.time1 data.time2 data.mat1 data.mat2 data.topoplot

% Initialize objects:

% a) Permutation test:
data.time1 = cell(job.nSub,1);
data.time2 = cell(job.nSub,1);

% b) 2-line plot:
data.mat1 = nan(1,length(data.mu.time)); 
data.mat2 = nan(1,length(data.mu.time)); 

% c) Topoplot:
data.topo2plot = data.Cond{2}; 

%% Compute objects WITH LOOPING over subjects (i.e. permutation test):

for iSub = 1:job.nSub
    
    if isfield(job,'bin')

        if strcmp(job.bin.Type,'RT')
            if job.bin.Num == 2
                data.time1{iSub} = data.ERPdata{iSub,1}; % initialize low.
                data.time1{iSub}.avg = real(data.ERPdata{iSub,1}.avg + data.ERPdata{iSub,3}.avg)./2; % mean of low
                data.time2{iSub} = data.ERPdata{iSub,2}; % initialize high.
                data.time2{iSub}.avg = real(data.ERPdata{iSub,2}.avg + data.ERPdata{iSub,4}.avg)./2; % mean of high

            elseif job.bin.Num == 3
                data.time1{iSub} = data.ERPdata{iSub,1}; % initialize low.
                data.time1{iSub}.avg = real(data.ERPdata{iSub,1}.avg + data.ERPdata{iSub,4}.avg)./2; % mean of low
                data.time2{iSub} = data.ERPdata{iSub,2}; % initialize medium.
                data.time2{iSub}.avg = real(data.ERPdata{iSub,2}.avg + data.ERPdata{iSub,5}.avg)./2; % mean of medium
                data.time3{iSub} = data.ERPdata{iSub,3}; % initialize high.
                data.time3{iSub}.avg = real(data.ERPdata{iSub,3}.avg + data.ERPdata{iSub,6}.avg)./2; % mean of high

            else
                error('Bin number %d not implemented yet',job.bin.Num);
            end
        
        else % if BOLD
            
            data.time1{iSub} = data.ERPdata{iSub,1}; % initialize high.
            data.time2{iSub} = data.ERPdata{iSub,2}; % initialize low.
            
        end
            
    elseif strcmp(job.contrastType,'Congruency')

        if strcmp(job.responseSettings,'Go')
        
            if strcmp(job.accSettings,'correct') || strcmp(job.accSettings,'incorrect') 
                data.time1{iSub} = data.ERPdata{iSub,1}; % initialize congruent.
                data.time1{iSub}.avg = real(data.ERPdata{iSub,1}.avg + data.ERPdata{iSub,4}.avg)./2; % mean of Go2Win and NoGo2Avoid
                data.time2{iSub} = data.ERPdata{iSub,2}; % initialize incongruent.
                data.time2{iSub}.avg = real(data.ERPdata{iSub,2}.avg + data.ERPdata{iSub,3}.avg)./2; % mean of Go2Avoid and NoGo2Win
            
            elseif strcmp(job.accSettings,'bothAcc') 
                data.time1{iSub} = data.ERPdata{iSub,1}; % initialize congruent.
                data.time1{iSub}.avg = real(data.ERPdata{iSub,1}.avg + data.ERPdata{iSub,4}.avg  + data.ERPdata{iSub,5}.avg + data.ERPdata{iSub,8}.avg)./4; % mean of Go2Win and NoGo2Avoid
                data.time2{iSub} = data.ERPdata{iSub,2}; % initialize incongruent.
                data.time2{iSub}.avg = real(data.ERPdata{iSub,2}.avg + data.ERPdata{iSub,3}.avg + data.ERPdata{iSub,6}.avg + data.ERPdata{iSub,7}.avg)./4; % mean of Go2Avoid and NoGo2Win
            
            else
                error('Unknown accuracy setting')
            end
            
        else
            error('Unknown response setting')
        end
        
    elseif strcmp(job.contrastType,'Go')
    
        if strcmp(job.responseSettings,'Go')
            data.time1{iSub} = data.ERPdata{iSub,1}; % initialize Go.
            data.time1{iSub}.avg = real(data.ERPdata{iSub,1}.avg + data.ERPdata{iSub,2}.avg)./2; % mean of Go2Win and NoGo2Avoid
            data.time2{iSub} = data.ERPdata{iSub,2}; % initialize NoGo.
            data.time2{iSub}.avg = real(data.ERPdata{iSub,3}.avg + data.ERPdata{iSub,4}.avg)./2; % mean of Go2Avoid and NoGo2Win

        elseif strcmp(job.responseSettings,'Hand')
            data.time1{iSub} = data.ERPdata{iSub,1}; % initialize Go.
            data.time1{iSub}.avg = real(data.ERPdata{iSub,1}.avg + data.ERPdata{iSub,2}.avg + data.ERPdata{iSub,3}.avg + data.ERPdata{iSub,4}.avg)./4;
            data.time2{iSub} = data.ERPdata{iSub,2}; % initialize NoGo.
            data.time2{iSub}.avg = real(data.ERPdata{iSub,5}.avg + data.ERPdata{iSub,6}.avg)./2;
        
        else
            error('Unknown response setting')
        end
        
    elseif strcmp(job.contrastType,'Valence')
        
        if strcmp(job.responseSettings,'Go')
            data.time1{iSub} = data.ERPdata{iSub,1}; % initialize Win.
            data.time1{iSub}.avg = real(data.ERPdata{iSub,1}.avg + data.ERPdata{iSub,3}.avg)./2; % mean of Go2Win and NoGo2Avoid
            data.time2{iSub} = data.ERPdata{iSub,2}; % initialize Avoid.
            data.time2{iSub}.avg = real(data.ERPdata{iSub,2}.avg + data.ERPdata{iSub,4}.avg)./2; % mean of Go2Avoid and NoGo2Win
        
        elseif strcmp(job.responseSettings,'Hand')
            data.time1{iSub} = data.ERPdata{iSub,1}; % initialize Win.
            data.time1{iSub}.avg = real((data.ERPdata{iSub,1}.avg + data.ERPdata{iSub,3}.avg)./2 + data.ERPdata{iSub,5}.avg)./2;
            data.time2{iSub} = data.ERPdata{iSub,3}; % initialize Avoid.
            data.time2{iSub}.avg = real((data.ERPdata{iSub,2}.avg + data.ERPdata{iSub,4}.avg)./2 + data.ERPdata{iSub,6}.avg)./2;
        
        else
            error('Unknown response setting')
        end
        
    elseif strcmp(job.contrastType,'GoValence')
        
        if strcmp(job.responseSettings,'Go')
            data.time1{iSub} = data.ERPdata{iSub,1}; % initialize Go2Win.
            data.time1{iSub}.avg = real(data.ERPdata{iSub,1}.avg); % only Go2Win
            data.time2{iSub} = data.ERPdata{iSub,2}; % initialize Go2Avoid.
            data.time2{iSub}.avg = real(data.ERPdata{iSub,2}.avg); % only Go2Avoid
        
        elseif strcmp(job.responseSettings,'Hand')
            data.time1{iSub} = data.ERPdata{iSub,1}; % initialize Go2Win.
            data.time1{iSub}.avg = real(data.ERPdata{iSub,1}.avg + data.ERPdata{iSub,3}.avg)./2; % only Go2Win
            data.time2{iSub} = data.ERPdata{iSub,3}; % initialize Go2Avoid.
            data.time2{iSub}.avg = real(data.ERPdata{iSub,2}.avg + data.ERPdata{iSub,4}.avg)./2; % only Go2Avoid
        
        else
            error('Unknown response setting')
        end
        
    elseif strcmp(job.contrastType,'GoAccuracy')
        
        if strcmp(job.responseSettings,'Go')
            data.time1{iSub} = data.ERPdata{iSub,1}; % initialize GoCorrect.
            data.time1{iSub}.avg = real(data.ERPdata{iSub,1}.avg + data.ERPdata{iSub,2}.avg)./2; % GoCorrect
            data.time2{iSub} = data.ERPdata{iSub,5}; % initialize GoIncorrect.
            data.time2{iSub}.avg = real(data.ERPdata{iSub,5}.avg + data.ERPdata{iSub,6}.avg)./2; % only GoIncorrect
        
        elseif strcmp(job.responseSettings,'Hand')
            data.time1{iSub} = data.ERPdata{iSub,1}; % initialize GoCorrect.
            data.time1{iSub}.avg = real(data.ERPdata{iSub,1}.avg + data.ERPdata{iSub,2}.avg + data.ERPdata{iSub,3}.avg + data.ERPdata{iSub,4}.avg)./4; % mean of LeftGo2Win and LeftGo2Avoid
            data.time2{iSub} = data.ERPdata{iSub,1}; % initialize GoIncorrect.
            data.time2{iSub}.avg = real(data.ERPdata{iSub,7}.avg + data.ERPdata{iSub,8}.avg + data.ERPdata{iSub,9}.avg + data.ERPdata{iSub,10}.avg)./4; % mean of RightGo2Win and RightGo2Avoid
        
        else
            error('Unknown response setting')
        end
        
    elseif strcmp(job.contrastType,'GoLeftRight') && strcmp(job.responseSettings,'Hand')
    
        if strcmp(job.accSettings,'correct') || strcmp(job.accSettings,'incorrect')
            data.time1{iSub} = data.ERPdata{iSub,1}; % initialize Left.
            data.time1{iSub}.avg = real(data.ERPdata{iSub,1}.avg + data.ERPdata{iSub,2}.avg)./2; % mean of LeftGo2Win and LeftGo2Avoid
            data.time2{iSub} = data.ERPdata{iSub,3}; % initialize Right.
            data.time2{iSub}.avg = real(data.ERPdata{iSub,3}.avg + data.ERPdata{iSub,4}.avg)./2; % mean of RightGo2Win and RightGo2Avoid
        
        elseif strcmp(job.accSettings,'bothAcc')
            data.time1{iSub} = data.ERPdata{iSub,1}; % initialize Left.
            data.time1{iSub}.avg = real(data.ERPdata{iSub,1}.avg + data.ERPdata{iSub,2}.avg + data.ERPdata{iSub,7}.avg + data.ERPdata{iSub,8}.avg)./4; % mean of LeftGo2Win and LeftGo2Avoid
            data.time2{iSub} = data.ERPdata{iSub,1}; % initialize Right.
            data.time2{iSub}.avg = real(data.ERPdata{iSub,3}.avg + data.ERPdata{iSub,4}.avg + data.ERPdata{iSub,9}.avg + data.ERPdata{iSub,10}.avg)./4; % mean of RightGo2Win and RightGo2Avoid
        
        else
            error('Unknown accuracy setting')
        end
        
    elseif strcmp(job.contrastType,'Accuracy') && strcmp(job.accSettings,'bothAcc')
        data.time1{iSub} = data.ERPdata{iSub,1}; % initialize incongruent.
        data.time1{iSub}.avg = real(data.ERPdata{iSub,1}.avg + data.ERPdata{iSub,2}.avg + data.ERPdata{iSub,3}.avg + data.ERPdata{iSub,4}.avg)./2; % mean of correct Go2Avoid and NoGo2Win
        data.time2{iSub} = data.ERPdata{iSub,1}; % initialize incongruent.
        data.time2{iSub}.avg = real(data.ERPdata{iSub,5}.avg + data.ERPdata{iSub,6}.avg + data.ERPdata{iSub,7}.avg + data.ERPdata{iSub,8}.avg)./2; % mean of correct Go2Avoid and NoGo2Win

    elseif strcmp(job.contrastType,'CongruentAccuracy') && strcmp(job.accSettings,'bothAcc')
        data.time1{iSub} = data.ERPdata{iSub,1}; % initialize incongruent.
        data.time1{iSub}.avg = real(data.ERPdata{iSub,1}.avg + data.ERPdata{iSub,4}.avg)./2; % mean of correct Go2Avoid and NoGo2Win
        data.time2{iSub} = data.ERPdata{iSub,1}; % initialize incongruent.
        data.time2{iSub}.avg = real(data.ERPdata{iSub,5}.avg + data.ERPdata{iSub,8}.avg)./2; % mean of incorrect Go2Avoid and NoGo2Win
    
    elseif strcmp(job.contrastType,'IncongruentAccuracy') && strcmp(job.accSettings,'bothAcc') 
        data.time1{iSub} = data.ERPdata{iSub,1}; % initialize incongruent.
        data.time1{iSub}.avg = real(data.ERPdata{iSub,2}.avg + data.ERPdata{iSub,3}.avg)./2; % mean of correct Go2Avoid and NoGo2Win
        data.time2{iSub} = data.ERPdata{iSub,1}; % initialize incongruent.
        data.time2{iSub}.avg = real(data.ERPdata{iSub,6}.avg + data.ERPdata{iSub,7}.avg)./2; % mean of incorrect Go2Avoid and NoGo2Win
    
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
            data.topo2plot.avg = (data.Cond{2}.avg+data.Cond{4}.avg)./2 ...
            - (data.Cond{1}.avg + data.Cond{3}.avg)./2; % high minus low

        elseif job.bin.Num == 3
            % a) 2-line plot (averaged over subjects):
            data.mat1 = squeeze((data.SubCondTime(:,1,:) + data.SubCondTime(:,4,:))/2); % low
%             data.mat2 = squeeze((data.SubCondTime(:,2,:) + data.SubCondTime(:,5,:))/2); % medium
            data.mat2 = squeeze((data.SubCondTime(:,3,:) + data.SubCondTime(:,6,:))/2); % high
            % b) Topoplot:
            data.topo2plot.avg = (data.Cond{3}.avg + data.Cond{6}.avg)./2 ...
            - (data.Cond{1}.avg + data.Cond{4}.avg)./2; % high minus low

        else
            error('Bin number %d not implemented yet',job.bin.Num);
        end
        
    else % if BOLD
        
        % a) 2-line plot (averaged over subjects):        
        data.mat1 = squeeze(data.SubCondTime(:,1,:)); % high
        data.mat2 = squeeze(data.SubCondTime(:,2,:)); % low
        % b) Topoplot:
        data.topo2plot.avg = data.Cond{1}.avg...
            - data.Cond{2}.avg; % high minus low

    end
    
elseif strcmp(job.contrastType,'Congruency')
        
    if strcmp(job.responseSettings,'Go') 

        if strcmp(job.accSettings,'correct') || strcmp(job.accSettings,'incorrect')
            % a) 2-line plot (averaged over subjects):
            data.mat1 = squeeze((data.SubCondTime(:,2,:) + data.SubCondTime(:,3,:))/2); % incongruent
            data.mat2 = squeeze((data.SubCondTime(:,1,:) + data.SubCondTime(:,4,:))/2); % congruent
            % b) Topoplot:
            data.topo2plot.avg = (data.Cond{2}.avg+data.Cond{3}.avg)./2 ...
            - (data.Cond{1}.avg + data.Cond{4}.avg)./2; % incongruent - congruent
            
        elseif strcmp(job.accSettings,'bothAcc')
            % a) 2-line plot (averaged over subjects):
            data.mat1 = squeeze((data.SubCondTime(:,2,:) + data.SubCondTime(:,3,:) + data.SubCondTime(:,6,:) + data.SubCondTime(:,7,:))/4); % incongruent
            data.mat2 = squeeze((data.SubCondTime(:,1,:) + data.SubCondTime(:,4,:) + data.SubCondTime(:,5,:) + data.SubCondTime(:,8,:))/4); % congruent
            % b) Topoplot:
            data.topo2plot.avg = (data.Cond{2}.avg + data.Cond{3}.avg + data.Cond{6}.avg + data.Cond{7}.avg)./4 ...
            - (data.Cond{1}.avg + data.Cond{4}.avg + data.Cond{5}.avg + data.Cond{8}.avg)/4;
            
        else
            error('Invalid accuracy setting')
        end
        
    elseif strcmp(job.responseSettings,'Hand')
    
        if strcmp(job.accSettings,'correct') || strcmp(job.accSettings,'incorrect') 
            % a) 2- line plot:
            data.mat1 = squeeze(((data.SubCondTime(:,2,:) + data.SubCondTime(:,4,:))/2 + data.SubCondTime(:,5,:))/2); % incongruent
            data.mat2 = squeeze(((data.SubCondTime(:,1,:) + data.SubCondTime(:,3,:))/2 + data.SubCondTime(:,6,:))/2); % congruent 
            % b) Topoplot:
            data.topo2plot.avg = ((data.Cond{2}.avg + data.Cond{4}.avg)./2 + data.Cond{5}.avg)./2 ...
            - ((data.Cond{1}.avg + data.Cond{3}.avg)./2 + data.Cond{6}.avg)./2;
            
        elseif strcmp(job.accSettings,'bothAcc')
            % a) 2- line plot:
            data.mat1 = squeeze(((data.SubCondTime(:,2,:) + data.SubCondTime(:,4,:))/2 + data.SubCondTime(:,5,:) + (data.SubCondTime(:,8,:) + data.SubCondTime(:,10,:))/2 + data.SubCondTime(:,11,:))/4); % incongruent
            data.mat2 = squeeze(((data.SubCondTime(:,1,:) + data.SubCondTime(:,3,:))/2 + data.SubCondTime(:,6,:) + (data.SubCondTime(:,7,:) + data.SubCondTime(:,9,:))/2 + data.SubCondTime(:,12,:))/4); % congruent 
            % b) Topoplot:
            data.topo2plot.avg = ((data.Cond{2}.avg + data.Cond{4}.avg)./2 + data.Cond{5}.avg + (data.Cond{8}.avg + data.Cond{10}.avg)./2 + data.Cond{11}.avg)./4 ...
                - ((data.Cond{1}.avg + data.Cond{3}.avg)./2 + data.Cond{6}.avg + (data.Cond{7}.avg + data.Cond{9}.avg)./2 + data.Cond{12}.avg)./4;
            
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
        data.topo2plot.avg = (data.Cond{1}.avg + data.Cond{2}.avg)./2 ...
        - (data.Cond{3}.avg + data.Cond{4}.avg)/2;
        
    elseif strcmp(job.responseSettings,'Hand')
        % a) 2-D line:
        data.mat1 = squeeze((data.SubCondTime(:,1,:) + data.SubCondTime(:,2,:) + data.SubCondTime(:,3,:) + data.SubCondTime(:,4,:))/4); % Go 
        data.mat2 = squeeze((data.SubCondTime(:,5,:) + data.SubCondTime(:,6,:))/2); % NoGo
        % b) Topoplot:
        data.topo2plot.avg = (data.Cond{1}.avg + data.Cond{2}.avg + data.Cond{3}.avg + data.Cond{4}.avg)./4 ...
        - (data.Cond{5}.avg + data.Cond{6}.avg)./2;
        
    else
        error('Unknown response setting')
    end
    
elseif strcmp(job.contrastType,'Valence')
        
    if strcmp(job.responseSettings,'Go')
        % a) 2-line plot
        data.mat1 = squeeze((data.SubCondTime(:,1,:) + data.SubCondTime(:,3,:))/2);
        data.mat2 = squeeze((data.SubCondTime(:,2,:) + data.SubCondTime(:,4,:))/2);
        % b) Topoplot:
        data.topo2plot.avg = (data.Cond{1}.avg+data.Cond{3}.avg)./2 ...
        - (data.Cond{2}.avg + data.Cond{4}.avg)./2;
        
    elseif strcmp(job.responseSettings,'Hand')
        % a) 2-line plot:
        data.mat1 = squeeze(((data.SubCondTime(:,1,:) + data.SubCondTime(:,3,:))/2 + data.SubCondTime(:,5,:))/2); 
        data.mat2 = squeeze(((data.SubCondTime(:,2,:) + data.SubCondTime(:,4,:))/2 + data.SubCondTime(:,6,:))/2); 
        % b) Topoplot:
        data.topo2plot.avg = ((data.Cond{1}.avg + data.Cond{3}.avg)./2 + data.Cond{5}.avg)./2 ...
        - ((data.Cond{2}.avg + data.Cond{4}.powspctr)./2 + data.Cond{6}.avg)./2;
        
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
            data.topo2plot.avg = data.Cond{1}.avg - data.Cond{2}.avg;
            
        elseif strcmp(job.accSettings,'bothAcc')
            % b) 2-line plot:
            data.mat1 = squeeze((data.SubCondTime(:,1,:) + data.SubCondTime(:,5,:))/2); % Go2Win
            data.mat2 = squeeze((data.SubCondTime(:,2,:) + data.SubCondTime(:,6,:))/2); % Go2Avoid
            % c) Topoplot:
            data.topo2plot.avg = data.Cond{1}.avg - data.Cond{2}.avg;
            
        else
            error('unknown accuracy setting')
        end
        
    elseif strcmp(job.responseSettings,'Hand')
        % a) 2-line plot:
        data.mat1 = squeeze((data.SubCondTime(:,1,:) + data.SubCondTime(:,3,:))/2); % Go2Win
        data.mat2 = squeeze((data.SubCondTime(:,2,:) + data.SubCondTime(:,4,:))/2); % Go2Avoid
        % b) Topoplot:
        data.topo2plot.avg = (data.Cond{1}.avg + data.Cond{3}.avg)./2 ...
            - (data.Cond{2}.avg + data.Cond{4}.avg)./2;
        
    else
        error('Unknown response setting')
    end
    
elseif strcmp(job.contrastType,'GoAccuracy') && strcmp(job.accSettings,'bothAcc')
        
    if strcmp(job.responseSettings,'Go')
        % a) 2-line plot:
        data.mat1 = squeeze((data.SubCondTime(:,1,:) + data.SubCondTime(:,2,:))/2); % Go correct 
        data.mat2 = squeeze((data.SubCondTime(:,5,:) + data.SubCondTime(:,6,:))/2); % Go incorrect
        % b) Topoplot:
        data.topo2plot.avg = (data.Cond{1}.avg + data.Cond{2}.avg)./2 ...
        - (data.Cond{5}.avg + data.Cond{6}.avg)./2;
        
    elseif strcmp(job.responseSettings,'Hand')
        % Decide whether to collapse left/right before averaging or not:
        % a) 2-line plot:
        data.mat1 = squeeze((data.SubCondTime(:,1,:) + data.SubCondTime(:,2,:) + data.SubCondTime(:,3,:) + data.SubCondTime(:,4,:))/4);
        data.mat2 = squeeze((data.SubCondTime(:,7,:) + data.SubCondTime(:,8,:) + data.SubCondTime(:,9,:) + data.SubCondTime(:,10,:))/4);
        % b) Topoplot:
        data.topo2plot.avg = (data.Cond{1}.avg + data.Cond{2}.avg + data.Cond{3}.avg + data.Cond{4}.avg)./4 ...
        - (data.Cond{7}.avg + data.Cond{8}.avg + data.Cond{9}.avg + data.Cond{10}.avg)./4;
        
    else
        error('invalid response settings')
    end
    
elseif strcmp(job.contrastType,'GoLeftRight') && strcmp(job.responseSettings,'Hand')
        
    if strcmp(job.accSettings,'correct') || strcmp(job.accSettings,'incorrect')
        % a) 2-line plot:
        data.mat1 = squeeze((data.SubCondTime(:,1,:) + data.SubCondTime(:,2,:))/2); % LeftGo
        data.mat2 = squeeze((data.SubCondTime(:,3,:) + data.SubCondTime(:,4,:))/2); % Right Go
        % b) Topoplot:
        data.topo2plot.avg = (data.Cond{1}.avg + data.Cond{2}.avg)./2 ...
        - (data.Cond{3}.avg + data.Cond{4}.avg)./2; % Left minus Right
        
    elseif strcmp(job.accSettings,'bothAcc')
        % a) 2-line plot:
        data.mat1 = squeeze((data.SubCondTime(:,1,:) + data.SubCondTime(:,2,:) + data.SubCondTime(:,7,:) + data.SubCondTime(:,8,:))/4); % Goleft
        data.mat2 = squeeze((data.SubCondTime(:,3,:) + data.SubCondTime(:,4,:) + data.SubCondTime(:,9,:) + data.SubCondTime(:,10,:))/4); % GoRight
        % b) Topoplot:
        data.topo2plot.avg = (data.Cond{1}.avg + data.Cond{2}.avg + data.Cond{7}.avg + data.Cond{8}.avg)./4 ...
        - (data.Cond{3}.avg + data.Cond{4}.avg + data.Cond{9}.avg + data.Cond{10}.avg)./4;
       
    else
        error('Unknown accuracy setting')
    end
    
elseif strcmp(job.contrastType,'Accuracy') && strcmp(job.accSettings,'bothAcc')
        
    if strcmp(job.responseSettings,'Go')
        % a) 2-line plot:
        data.mat1 = squeeze((data.SubCondTime(:,1,:) + data.SubCondTime(:,2,:) + data.SubCondTime(:,3,:) + data.SubCondTime(:,4,:))/4); % correct 
        data.mat2 = squeeze((data.SubCondTime(:,5,:) + data.SubCondTime(:,6,:) + data.SubCondTime(:,7,:) + data.SubCondTime(:,8,:))/4); % incorrect
        % b) Topoplot:
        data.topo2plot.avg = (data.Cond{1}.avg + data.Cond{2}.avg + data.Cond{3}.avg + data.Cond{4}.avg)./4 ...
        - (data.Cond{5}.avg + data.Cond{6}.avg + data.Cond{7}.avg + data.Cond{8}.avg)./4;
        
    elseif strcmp(job.responseSettings,'Hand')
        % Decide whether to collapse left/right before averaging or not:
        % a) 2-line plot:
        data.mat1 = squeeze((data.SubCondTime(:,1,:) + data.SubCondTime(:,2,:) + data.SubCondTime(:,3,:) + data.SubCondTime(:,4,:) + data.SubCondTime(:,5,:) + data.SubCondTime(:,6,:))/6);
        data.mat2 = squeeze((data.SubCondTime(:,7,:) + data.SubCondTime(:,8,:) + data.SubCondTime(:,9,:) + data.SubCondTime(:,10,:) + data.SubCondTime(:,11,:) + data.SubCondTime(:,12,:))/6);
        % b) Topoplot:
        data.topo2plot.avg = (data.Cond{1}.avg + data.Cond{2}.avg + data.Cond{3}.avg + data.Cond{4}.avg + data.Cond{5}.avg + data.Cond{6}.avg)./6 ...
        - (data.Cond{7}.avg + data.Cond{8}.avg + data.Cond{9}.avg + data.Cond{10}.avg + data.Cond{11}.avg + data.Cond{12}.avg)./6;
        
    else
        error('invalid response settings')
    end
    
elseif strcmp(job.contrastType,'CongruentAccuracy') && strcmp(job.accSettings,'bothAcc')
        
    if strcmp(job.responseSettings,'Go')
        % a) 2-line plot:
        data.mat1 = squeeze((data.SubCondTime(:,1,:) + data.SubCondTime(:,4,:))/2); % congruent correct
        data.mat2 = squeeze((data.SubCondTime(:,5,:) + data.SubCondTime(:,8,:))/2); % congruent incorrect
        % b) Topoplot:
        data.topo2plot.avg = (data.Cond{1}.avg+data.Cond{4}.avg)./2 ...
        - (data.Cond{5}.avg + data.Cond{8}.avg)./2;
        
    elseif strcmp(job.responseSettings,'Hand')
        % a) 2-line plot:
        data.mat1 = squeeze(((data.SubCondTime(:,1,:) + data.SubCondTime(:,3,:))/2 + data.SubCondTime(:,6,:))/2); % congruent correct
        data.mat2 = squeeze(((data.SubCondTime(:,7,:) + data.SubCondTime(:,9,:))/2 + data.SubCondTime(:,12,:))/2); % congruent incorrect
        % b) Topoplot:
        data.topo2plot.avg = ((data.Cond{1}.avg + data.Cond{3}.avg)./2 + data.Cond{6}.avg)./2 ...
            - ((data.Cond{7}.avg + data.Cond{9}.avg)./2 + data.Cond{12}.avg)./2;        
    else
        error('unknown response settings')
    end
    
elseif strcmp(job.contrastType,'IncongruentAccuracy') && strcmp(job.accSettings,'bothAcc')
    
    job.twoLineLabels = {'Incongruent correct','Incongruent incorrect'};
    
    if strcmp(job.responseSettings,'Go')
        % a) 2-line plot:
        data.mat1 = squeeze((data.SubCondTime(:,2,:) + data.SubCondTime(:,3,:))/2); % incongruent correct
        data.mat2 = squeeze((data.SubCondTime(:,6,:) + data.SubCondTime(:,7,:))/2); % incongruent incorrect
        % b) Topoplot:
        data.topo2plot.avg = (data.Cond{2}.avg+data.Cond{3}.avg)./2 ...
        - (data.Cond{6}.avg + data.Cond{7}.avg)./2;
        
    elseif strcmp(job.responseSettings,'Hand')
        % a) 2-line plot:
        data.mat1 = squeeze(((data.SubCondTime(:,2,:) + data.SubCondTime(:,4,:))/2 + data.SubCondTime(:,5,:))/2); % incongruent correct
        data.mat2 = squeeze(((data.SubCondTime(:,8,:) + data.SubCondTime(:,10,:))/2 + data.SubCondTime(:,11,:))/2); % incongruent correct
        % b) Topoplot:
        data.topo2plot.avg = ((data.Cond{2}.avg + data.Cond{4}.avg)./2 + data.Cond{5}.avg)./2 ...
            - ((data.Cond{8}.avg + data.Cond{10}.avg)./2 + data.Cond{11}.avg)./2;
        
    else
        error('unknown response settings')
    end
    
elseif strcmp(job.contrastType,'CorrectIncongruentIncorrectCongruent') && strcmp(job.accSettings,'bothAcc')
    
    job.twoLineLabels = {'Incongruent correct','Congruent incorrect'};
    
    if strcmp(job.responseSettings,'Go')
        % a) 2-line plot:
        data.mat1 = squeeze((data.SubCondTime(:,2,:) + data.SubCondTime(:,3,:))/2); % incongruent correct
        data.mat2 = squeeze((data.SubCondTime(:,5,:) + data.SubCondTime(:,8,:))/2); % incongruent incorrect
        % b) Topoplot:
        data.topo2plot.avg = (data.Cond{2}.avg+data.Cond{3}.avg)./2 ...
        - (data.Cond{5}.avg + data.Cond{8}.avg)./2;
        
    elseif strcmp(job.responseSettings,'Hand')
        error('Not yet implemented')
    end
else 
    error('Unknown contrast type')
end

%% Error bounds for 2-line plots:

data.matDiff        = data.mat1 - data.mat2; % for difference-CIs
data.matMean        = (data.mat1 + data.mat2)./2; % for Cousineau method
data.matGrandMean   = nanmean(data.matMean); % for Cousineau method
data.SubTime        = squeeze(nanmean(data.SubCondTime,2)); % average over conditions for Cousineau method
data.GrandTime      = squeeze(nanmean(data.SubTime(job.validSubs,:),1)); % average over subjects for Cousineau method

fprintf('Finished preparing data for contrast %s, channels %s\n',job.contrastType,strjoin(job.channels,'/'))
% END
