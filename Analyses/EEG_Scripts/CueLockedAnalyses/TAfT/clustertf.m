function [corrp, tg] = clustertf(c,thresh,nP,teststat)

% [corrp, tg] = clustertf(c,thresh,nP,teststat)
% 
% Take 4D input object (subject x regressors x frequency x time), loop over
% regressors, perform one-sample t-test across subjects fore each frequency
% x time grid point, do cluster-based permutation test for grid points
% above threshold, output t-value and corrected p-value for each grid
% point.
%
% INPUTS:
% c             = 4D matrix of TF data with subject in first dimension,
% unclustered 'regressors' (e.g. channels) in second dimension, frequency
% in third and time in fourth dimension.
% thresh        = numeric scalar, T-value cut-off used for thresholding
% (default: 2)
% nP            = numeric scalar, number of permutations (default: 500).
% teststat      = string, test statistic to use, either 'nExtent' (number
% grid points above threshold forming cluster), 'maxT' (maximum T-value per
% cluster), 'sumT' (sum of all T-values above threshold forming cluster;
% default).
%
% OUTPUTS:
% corrp         = 4D matrix with corrected p-value for each data point.
% tg            = 4D matrix with t-value for each data point.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% Laurence Hunt, 2013.
% Johannes Algermissen, 2019.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/CueLockedAnalyses/TAfT/

%% Complete settings:

if nargin<2
  thresh = 2; % see Fieldtrip defaults
end

if nargin<3
  nP = 500; % number permutations
end

if nargin<4
    teststat = 'sumT'; % 'nExtent' or 'maxT' or 'sumT' 
end

%% Initialize and reshape input:

fprintf('Reshape input\n')

nS = size(c,1); %number of 'subjects' (recording sites)
nR = size(c,2); %number of regressors (e.g. channels -- anything spatially unclustered)
nF = size(c,3); %number of frequencies
nT = size(c,4); %number of timebins

cr = reshape(c,nS,nR*nF*nT); % concatenate regressors/frequencies/timebins into 2-D

%% 2) Build up a null distribution of the maximum cluster size for each permutation i and each regressor r:

fprintf('Build null distribution\n')
nulldist = nan(nP,nR); % initialize null distribution

for i = 1:nP % for each permutation

    %% 1) Design matrix:
    
    % a) Random with replacement:
    dm = ((rand(nS,1)>0.5)-0.5)*2; % draw from random uniform distribution between 0 and 1, threshold at 0.5 to 0 and 1, minus 0.5, times 2 gives -1 or +1
    % --> assigns subjects randomly to +1 or -1 (not necessarily 50:50 in each permutation)

    % b) Random without replacement:
%     dm = ones(nS,1); % assign all to 1
%     dm(1:ceil(nS/2)) = -1; % assign half to -1
%     dm = dm(randperm(nS)); % randomly permute order
    % --> ensures that exactly 50:50 is +1 and -1

    %% 2) T-test:
     
    [~,~,tg] = ols(cr,dm,eye(size(dm,2))); % compute t-value for each data point

    %% 3) Reshape:
    
    tg = reshape(tg,size(dm,2),nR,nF,nT); % reshape T-value output back into 4-D (subject, regressor, frequency, time)   
    %    figure;imagesc(squeeze(tg));colorbar % plot T-values

    for r = 1:size(tg,2) % for each regressor:
        
        % a) Threshold t-values
        tgr = squeeze(tg(1,r,:,:)); % extract frequency x time for relevant regressor
        tgr_thresh = abs(tgr)>thresh; % classify whether t-value above threshold or not        
        
        % b) Identify clusters of adjacent t-values above threshold:
        tgr_thresh = taft_prune(tgr_thresh); % delete single-voxel connections
        imlabel = bwlabel(tgr_thresh); % labels for connected components (i.e. clusters) above thres
        tmp = unique(imlabel); % remove doublings in cluster labels
        tmp(tmp==0) = []; % delete zeros
        
        % c) Loop over clusters and compute cluster statistic 
        nL = 0; % initialize test statistic vector
        
        if ~isempty(tmp) % if any connected clusters above threshold:
            for k = tmp' % For every cluster above threshold:
                
                % i) Extract relevant statistics:
                clusterT = abs(tgr(imlabel==k)); % extract all T values
                nExtent = length(clusterT(:)); % cluster size
                maxT = max(clusterT(:)); % maximum T
                sumT = sum(clusterT(:)); % sum of all Ts (cluster mass)
                
                % ii) Determine which one is used as test statistic:
                if strcmp(teststat,'nExtent')
                    nL(k) = nExtent; % cluster extent in # voxels
                elseif strcmp(teststat,'maxT')
                    nL(k) = maxT; % maximum T-value
                elseif strcmp(teststat,'sumT')
                    nL(k) = sumT; % cluster mass
                else
                    error('Invalid test statistic specification')
                end % end if
            end % end k
        end % end isempty
        
        % iii) Add to null distribution
        nulldist(i,r)=max(nL); % add maximal (!) cluster criterion of permutation
        
    end
end

%% 3) Run a one sample t-test on data and compare cluster criterion to permutation data:

% a) Compute t-values for actual data:

fprintf('Perform one-sample t-test, compare cluster criterion to null distribution\n')

dm = ones(nS,1); % design matrix: mean of all subjects (actual dm of interest)
[~,~,tg] = ols(cr,dm,eye(size(dm,2))); % compute t-value for each data point
tg = reshape(tg,size(dm,2),nR,nF,nT); % reshape back into 4-D

% b) Compare to null distribution:

corrp = ones(size(tg)); % initialize vector for corrected p-value for each data point

for r = 1:size(tg,2) % for each regressor

    tgr = squeeze(tg(1,r,:,:)); % extract frequency x time for relevant regressor
    fprintf('Regressor %d: Search space is %d voxels\n',r,size(tgr,1)*size(tgr,2));
    fprintf('Regressor %d: Maximum of null distribution: %f \n',r,max(nulldist(:,r)));
    fprintf('Regressor %d: %d voxels < -thres, %d voxels > thres\n',r,sum(tgr(:)< -1*thresh),sum(tgr(:)> 1*thresh));
    
    imlabel = bwlabel(abs(tgr)>thresh); % labels for connected components (i.e. clusters) above thres
    tmp = unique(imlabel); % remove doublings in cluster labels
    tmp(tmp==0) = []; % delete zeros
    
    if ~isempty(tmp) % if any clusters:
        
        fprintf('Regressor %d: Found %d clusters above threshold\n',r,length(tmp));
        cpimlabel = ones(size(imlabel)); % initialize 2-D matrix of corrected p-values per regressors

        for k = tmp' % For each cluster

            % i) Extract relevant statistics:
            nExtent = sum(sum(imlabel==k)); % cluster size
            clusterT = abs(tgr(imlabel==k)); % extract all T values of cluster
            maxT = max(clusterT(:)); % maximum T value
            sumT = sum(clusterT(:)); % sum of all T-values (cluster mass)
            
            % ii) Determine which one is used as test statistic: 
            if strcmp(teststat,'nExtent')
                cp = mean(nExtent<=nulldist(:,r)); % how often cluster size (or bigger size) in null distribution?
            elseif strcmp(teststat,'maxT')
                cp = mean(maxT<=nulldist(:,r)); % how often maximal T-value (or bigger T-value) in null distribution?
            elseif strcmp(teststat,'sumT')
                cp = mean(sumT<=nulldist(:,r)); % how often cluster mass (or bigger mass) in null distribution?
            else
                error('Invalid test statistic specification')
            end
            
            % iii) Output relevant statistics for that cluster
            fprintf('Cluster %d: %d voxels, max|T| = %.02f, sum|T| = %.02f, p = %.03f\n',k,nExtent,maxT,sumT,cp);
            cpimlabel(imlabel==k)=cp; % attach corrected p-value to all voxels in cluster
            
        end % end k
    end % end isempty

    % iv) Add to map of corrected p-values corrp:
    corrp(1,r,:,:) = cpimlabel; % add corrected p-value to each data point of regressor

end % end of for-loop for regressors.

end % end of function.