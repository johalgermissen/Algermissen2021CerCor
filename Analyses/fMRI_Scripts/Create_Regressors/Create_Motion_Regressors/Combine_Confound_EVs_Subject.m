function Combine_Confound_EVs_Subject(iGLM)

% Combine_Confound_EVs_Subject(iGLM).
% 
% Will take 
% - 6 realignment parameters
% - cerebrospinal fluid (CSF) signal
% - out-of-brain (OOB) signal
% - spikes in relative displacement
% and combine those into a single matrix.
% Will concatenate blocks to yield single matrix per subject.
%
% Mind setting root directory rootdir.
%
% INPUT:
% iGLM              = string, GLM identifier where to save matrix.
%
% OUTPUT:
% Will save matrix as tConfoundEVs.txt in timing_regressors folder of
% specific GLM of specific subject.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

if nargin < 1
    iGLM    = '1'; % new (subject-wise) GLM to write to
    fprintf('No GLM specified as target; use default %s \n',iGLM);
end

%% Set root directory:

fprintf('Set root directory\n')
rootdir = '/project/3017042.02';

%% Fixed settings:

% clear all; clc

nSub    = 36;
nBlock  = 6;
cutoff  = 2; % cutoff for relative displacement counted as spike

% Invalid subjects to skip:
invalidSubs = [15 25];
validSubs   = setdiff(1:nSub,invalidSubs);
fprintf('Skip subjects %s\n',strjoin(string(invalidSubs),', '));

%% Loop over subjects and blocks:

for iSub = validSubs % iSub = 1;

    fprintf('Starting with subject %d\n', iSub)
    
    % Initialize empty objects:
    allPars_block  = cell(nBlock,1);
    intercept     = cell(nBlock,1);
    relDisplacement = nan(0,0);

    for iBlock = 1:nBlock % iBlock = 1;

        % Set directory of block of subject:
        blockDir        = fullfile(rootdir,sprintf('Log/fMRI/sub-%03d/FEAT_Block%d.feat/AROMA',iSub,iBlock));        
        subjectDir      = fullfile(rootdir,sprintf('Log/fMRI/sub-%03d/GLM%s/timings_regressors',iSub,iGLM));
        if ~exist(subjectDir,'dir'); mkdir(subjectDir); end % make directory if necessary

        %% Load and concatenate realignment parameters, CSF, OOB:
        
        % 1) Load 6 (raw) realignment parameters:
        fprintf('Load realignment parameters\n');
        realignName         = fullfile(rootdir,'Log/fMRI/sub-%03d/FEAT_Block%d.feat/mc/prefiltered_func_data_mcf.par',iSub,iBlock);        
        realignPars         = textread(realignName); % read in file with relative displacement parameters
        realignPars(:,7)    = []; % last column always only zeros

        % 2) Load mean CSF signal:
        fprintf('Load CSF signal\n');
        CSF_Regressor       = load(fullfile(blockDir,'CSF_noise.txt')); % mean CSF signal

        % 3) Load mean OOB signal:
        fprintf('Load OOB signal\n');
        OOB_Regressor       = load(fullfile(blockDir,'OOB_noise.txt')); % mean OOB signal

        % 4) Load relative displacement parameters:
        fprintf('Load relative displacement parameters\n');
        realignName          = fullfile(rootdir,sprintf('Log/fMRI/sub-%03d/FEAT_Block%d.feat/mc/prefiltered_func_data_mcf_rel.rms',iSub,iBlock)); % motion parameters
        relMotion           = textread(realignName); % read in file with relative displacement parameters

        % Always states motion during upcoming volume, thus add 0 at the end of block, as block transitions not taken into account
        relDisplacement         = [relDisplacement; relMotion; 0]; % keep relative displacement from previous block

        % Combine into one matrix:
        fprintf('Combine vectors\n');
        allPars                 = [realignPars CSF_Regressor OOB_Regressor]; %

        % Demean all regressors:
        fprintf('Demean vectors\n');
        allPars_demeaned        = allPars - mean(allPars,1);
        % round(mean(Realign_Demean,1),10) % check: should be zero
        allPars_block{iBlock}   = allPars_demeaned; % save for this block

        % 5) Create intercept for this block (1 for every volume):
        fprintf('Generate intercept per block\n');
        intercept{iBlock}       = ones(size(allPars_block,1),1); % save number of volumes of this block

    end % end iBlock

    % ------------------------------------------------------------------- %
    % Concatenate blocks diagonally, first motion/CSF/OOB, than intercepts
    fprintf('Diagonally concatenate signals and intercepts\n');
    confoundEVs               = [blkdiag(allPars_block{:}) blkdiag(intercept{:})];
    confoundEVs(:,end)        = []; % delete intercept for last block, otherwise rank-deficient design matrix

    % ------------------------------------------------------------------- %
    %% B) Create spike matrix (if spikes present):

    nVolumes    = size(relDisplacement,1); % number of volumes
    nSpikes     = sum(relDisplacement > cutoff); % number of spikes (volumes above cutoff)
    fprintf('Found %d motion spikes\n', nSpikes);

    if nSpikes > 0
        spikeMatrix = zeros(nVolumes,nSpikes); % initialize matrix with separate column for each spike
        iSpike      = 1; % initialize counter for spikes

        for iVol = 1:nVolumes

            if relDisplacement(iVol) > cutoff
                spikeMatrix(iVol,iSpike) = 1;
                iSpike = iSpike+1; % increment so to use the next colum for the next spike
            end

        end
        % sum(spikeMatrix,1) % should be 1 per column

        confoundEVs = [confoundEVs spikeMatrix];

    else % otherwise no spike regressors

        confoundEVs = confoundEVs;

    end % end if nSpikes > 0
    
    %% Some diagnostic plots:
    
%     % Diagnose rank deficiency:
%     M = confoundEVs;
%     M = [ones(size(confoundEVs,1),1) confoundEVs];
%     size(confoundEVs,2)
%     rank(confoundEVs)

    % Cross-correlation matrix:
%    M = confoundEVs;
%    M = blkdiag(Realign_CSFOOB{:});
%     for m1 = 1:size(M,2); % First variable
%         for m2 = 1:size(M,2); % Second variable
%         Cor(m2,m1) = corr(M(:,m1),M(:,m2)); % Correlation
%         end
%     end
%     colormap('jet');
%     imagesc(Cor)
%     colorbar

%     % Eigenvalue decomposition:
%     covmat = cov(M);
%     [evecs,evals] = eig(covmat);
%     [evals,I] = sort(diag(evals),'descend');
% %     evecs = evecs(:,I);
%     pve = cumsum(evals)/sum(evals)*100;
%     figure; plot(1:length(evals),pve); 
%     axis([0 length(evals) 0 100])
%     xlabel('number of eigenvectors');
%     ylabel('percent variance explained');

    %% Save as text file:
%     round(mean(confoundEVs,1),10)'
    save(fullfile(subjectDir,'tConfoundEVs.txt'),'confoundEVs','-ascii')
    
end % end iSub

end % end of function.