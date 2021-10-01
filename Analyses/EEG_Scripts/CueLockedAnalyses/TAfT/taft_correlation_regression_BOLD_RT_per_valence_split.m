function taft_correlation_regression_BOLD_RT_per_valence_split(ROI2use,name,normRT,withoutInvalid)

% taft_correlation_regression_BOLD_RT_per_valence_split(ROI2use,name,normRT,withoutInvalid)
%
% For selected ROI, loop over subjects, initialize TAfT job, upsample
% volume-by-volume data, epoch to trial-by-trial data, z-standardize BOLD
% and RTs separately per valence condition, extract mean BOLD and RT per
% BOLD tertile, perform correlation/ regression, run one-sample t-test 
% across subjects, plot RT as function of BOLD:
% - separate plots per valence, points + SE for mean per bin, individual
% subject means as points
% - one plot, points + SE for mean per bin per valence, individual subject
% means as points
% Mind setting root directory.
% Mind to also adjust settings in taft_preprocess_initialize_job if
% necessary.
%
% INPUTs:
% ROI2use       = cell with 1 string, 1 ROI to upsample, epoch, and use for correlation.
% name          = string, name of ROI to use in plot header/ file name.
% nBin          = scalar integer, number of bins to use (default: 3)
% normRT        = Boolean, whether to normalize RTs in plots or not (default: false).
% withoutInvalid= Boolean, whether to exclude all subjects exclude in
% fMRI-informed EEG analyses (TAfT) or not (default: false).
%
% OUTPUTS:
% none, save to disk, provide results to console, plot.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Adapted from Tobias Hauser (https://github.com/tuhauser/TAfT).
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/CueLockedAnalyses/TAfT/

% BOLD split into tertiles on x-axis,
% mean RT per vmPFC BOLD tertile on y-axis,
% points and lines (quantile plot)

%% Fill in defaults:

if ~exist('ROI2use','var')
    ROI2use     = {'GLM1vmPFCValenceMan'}; name = 'vmPFC'; % 
%     ROI2use    = {'GLM1CingulateAnteriorValence'}; name = 'ACC'; % 
    fprintf('Use ROI %s, name %s\n',ROI2use{:}, name),
end

if ~exist('nBin','var')
    nBin        = 3; % number bins for quantile plot
    fprintf('Use %d bins\n',nBin),
end

if ~exist('normRT','var')
    normRT      = false;
    fprintf('Do not normalize RTs in plots\n'),
end

if ~exist('withoutInvalid','var')
    withoutInvalid = false;
    fprintf('Do not exclude subjects that were excluded in TAfT\n');
end

%% Fixed settings:

% Fixed settings:
nSub        = 36; % number subjects
nTrial      = 640; % number trials
nVal        = 2; % number valence conditions

%% Set directories:

dirs.root	    = taft_set_rootdir(); % /project/3017042.02
dirs.EEG 	    = fullfile(dirs.root,'Log','EEG','CueLockedResults');
dirs.TAfT           = fullfile(dirs.EEG,'TAfT_Betas');
dirs.RTcorrelation  = fullfile(dirs.TAfT,'TAfT_RT_correlation');
if ~exist(dirs.RTcorrelation,'dir'); mkdir(dirs.RTcorrelation); end

if normRT
    suffix = '_normRT';
else
    suffix = '';
end

%% Create fileName:

fileName  = sprintf(fullfile(dirs.RTcorrelation,sprintf('TAfT_TF_correlation_RT_%s_split_nBin%d%s.mat',...
    string(ROI2use),nBin,suffix)));
fprintf('File name to save/retrieve is \n%s\n',fileName);

%% Loop over subjects, extract fMRI and behavior

if ~exist(fileName,'file') % If file does not exist yet:

    fprintf('Unknown file %s, start creating\n',fileName);
    
    % Initialize objects:
    corMat.cor      = nan(nSub,1);
    corMat.N        = nan(nSub,1);
    corMat.b        = nan(nSub,4);
    corMat.BOLD     = nan(nSub,nBin,nVal);
    corMat.RT       = nan(nSub,nBin,nVal);
    
    %% Loop over subjects, extract fMRI and behavior
    
    for iSub = 1:nSub % iSub = 1;
        
        % 1) Initialize job:
        job             = taft_preprocess_initialize_job('TF',iSub,ROI2use);
               
        % 2) Select valid trials:
        job.goodTrlIdx  = 1:640; % select all trials
        
        % 3) Load, fit, combine, and reshape fMRI-beta data of different ROIs and blocks:
        fprintf('Subject %03d: Load fMRI data\n',iSub);
        X               = taft_preprocess_load_fMRI(job);

        % 2) Load behavior:
        behav           = taft_preprocess_load_behavior(job.subID);

        % --------------------------------------------------------------- %
        % 5) Retrieve indices of trial types:
        fprintf('Identify indices of trial types\n');
        
        % Identify trial indices:
        Goidx          = (~isnan(behav.RT) & behav.isgo==1); % valid RTs and performed action is Go
%         Goidx      = find(~isnan(Y) & behav.reqaction==1); % valid RTs and required action is Go
        G2Widx          = Goidx & (behav.valence==1); % Win trials
        G2Aidx          = Goidx & (behav.valence==2); % Avoid trials
        fprintf('Found %d trials with RTs, thereof %d Win trials, %d Avoid trials\n',sum(Goidx),sum(G2Widx),sum(G2Aidx));

        % --------------------------------------------------------------- %
        % 6) For analyses, normalize RTs and BOLD signal for Win and Avoid trials separately, put back together:
        fprintf('### For analyses: normalize RT for Win and Avoid trials separately, pu back together\n');

        % a) For RTs (Y):
        % Alternative A: Standardize separately per valence:
        Y               = behav.RT; % extract RTs
        Ysep            = nan(nTrial,1); % empty vector
        Ysep(G2Widx)    = z_standardize(Y(G2Widx)); % extract Win, normalize, put back
        Ysep(G2Aidx)    = z_standardize(Y(G2Aidx)); % extract Avoid, normalize, put back
        Ynorm           = Ysep(Goidx); % 	store Go trials only

        % Alternative B: Standardize overall:
%         Ynorm           = z_standardize(Y(Goidx)); % for RT: extract Go trials, normalize
        
        % b) BOLD Signal:
        % Alternative A: Standardize separately per valence:
        xSep            = nan(nTrial,1); % empty vector
        xSep(G2Widx)    = z_standardize(X(G2Widx)); % extract Win, normalize, put back
        xSep(G2Aidx)    = z_standardize(X(G2Aidx)); % extract Avoid, normalize, put back
        xNorm           = xSep(Goidx); % 	store Go trials only

        % Alternative B: Standardize overall:
%         Xnorm           = z_standardize(X(Goidx)); % for BOLD: extract Go trials, normalize
        
        if length(Ynorm) ~= length(xNorm); error('X and Y of different length!'); end
        
        % c) Compute correlation:
        fprintf('Compute correlation \n');
        corMat.cor(iSub)  = corr(xNorm,Ynorm); % correlation
        
        % d) Save cell size (Go trials):
        fprintf('Store number of Go trials\n');
        corMat.N(iSub)    = sum(Goidx); % number of Go trials (with RTs)
        
        % e) Perform linear regression with BOLD, valence, interaction:
        fprintf('Compute linear regression\n');
        % Valence regressor:
        val               = 2 - behav.valence(Goidx); % valence effect: Win is 1; Avoid is 0
        valNorm           = z_standardize(val); % standardize valence effect
        % Interaction:
        xVal              = xNorm .* valNorm; % interaction
        xValNorm          = z_standardize(xVal); % interaction valence x BOLD
        % Design matrix:
        DM                = [ones(length(xNorm),1) xNorm valNorm xVal]; % concatenate to design matrix
%         DM                = [ones(length(Xnorm),1) xNorm valNorm xValNorm]; % standardized interaction term
        % Linear regression:
        linReg            = DM\Ynorm; % pinv regression
        % Store regression coefficients per condition per subject:
        corMat.b(iSub,:)  = linReg; % regression coefficients
                
        % --------------------------------------------------------------- %
        % 7) For plots, normalize both RTs and BOLD first, 
        % extract mean BOLD and RT per BOLD quantile:        
        
        fprintf('### For plots: normalize both RT and BOLD first\n');
        
        % a) BOLD:
        xNorm           = nan(nTrial,1);    % initialize
        xNorm(Goidx)    = X(Goidx);         % extract BOLD on Go trials; rest stays NaN
        xNorm           = z_standardize(xNorm); % normalize
        XWin            = xNorm(G2Widx);    % extract Go2Win trials
        XAvoid          = xNorm(G2Aidx);    % extract Go2Avoid trials

        % b) RT:
        if normRT
            Ynorm           = z_standardize(Y); % normalize
        else    
            Ynorm           = Y; % for RT: keep raw scores
        end
        
        % c) Determine quantile cut-offs:
        
        fprintf('Determine quantile cut-offs \n');

        G2Wcutoff = nan(1,nBin+1); 
        G2Acutoff = nan(1,nBin+1);
        
        for iBin = 0:nBin
                G2Wcutoff(iBin+1) = quantile(XWin,iBin/nBin); % quantile of iBin / nBin
                G2Acutoff(iBin+1) = quantile(XAvoid,iBin/nBin); % quantile of iBin / nBin
        end
        
        if G2Wcutoff(1) ~= min(XWin); error('G2W min not proper'); end
        if G2Wcutoff(end) ~= max(XWin);  error('G2W max not proper'); end
        if G2Acutoff(1) ~= min(XAvoid); error('G2A min not proper'); end
        if G2Acutoff(end) ~= max(XAvoid); error('G2A max not proper'); end
        
        % d) Indices of trials within each quantile:
        fprintf('Determine indices within each quantile\n');

        G2WbinIdx = nan(length(xNorm),nBin); % initialize
        G2AbinIdx = nan(length(xNorm),nBin); % initialize

        for iBin = 1:nBin % 
            G2WbinIdx(:,iBin) = (G2Widx & xNorm > G2Wcutoff(iBin) & xNorm <= G2Wcutoff(iBin+1));
            G2AbinIdx(:,iBin) = (G2Aidx & xNorm > G2Acutoff(iBin) & xNorm <= G2Acutoff(iBin+1));
        end
                
        % Create decision variable:
        % Xsel = G2WbinIdx * [1:nBin]' + G2AbinIdx * [(nBin+1):(2*nBin)]';
        % tabulate(Xsel) % plot per bin

        % f) Compute mean within each bin:
        fprintf('Compute mean for RT and BOLD within each quantile\n');
        for iBin = 1:nBin % iBin = 1;
            corMat.BOLD(iSub,iBin,1)    = nanmean(xNorm(G2WbinIdx(:,iBin)==1));
            corMat.BOLD(iSub,iBin,2)    = nanmean(xNorm(G2AbinIdx(:,iBin)==1));
            corMat.RT(iSub,iBin,1)      = nanmean(Ynorm(G2WbinIdx(:,iBin)==1));
            corMat.RT(iSub,iBin,2)      = nanmean(Ynorm(G2AbinIdx(:,iBin)==1));
        end % end iBin
        
    end % end iSub
    
    % Save to disk:
    fprintf('Save to disk under %s\n', fileName);
    save(fileName,'corMat','-v7.3')
    fprintf('Finished saving :-)\n')
    
else % if exist: just load
    
    fprintf('Found existing file %s, load\n',fileName);
    load(fileName);
    fprintf('Finished loading :-)\n')

end

%% Select valid subjects:

% Select valid subjects:
% invalidSubs = []; % include all subjects
invalidSubs = [15 25]; % bad co-registrations
if withoutInvalid
    invalidSubs = [1 11 15 19 21 25]; suffix = 'withoutInvalid'; % invalid co-registrations and outliers
end
validSubs   = setdiff(1:nSub,invalidSubs);
nSubValid   = length(validSubs);

fprintf('Select valid subjects, exclude subjects %s\n',strjoin(string(invalidSubs),', '));

%% Correlation:

corVec      = corMat.cor(validSubs);
fprintf('Region: %s \n',string(ROI2use));

% Fisher-z transform:
fprintf('Correlation:\n');
fprintf('Min = %.02f, Max = %.02f, Mean = %.02f\n',min(corVec),max(corVec),mean(corVec)) % before Fisher-z transform
fprintf('Fisher z-transform\n');
corVec     = atanh(corVec);

% T-test:
fprintf('Perform t-test across subjects: \n');
[~,P,~,STATS] = ttest(corVec);
fprintf('t(%d) = %.02f, p = %.03f, d = %.02f\n',STATS.df,STATS.tstat,P/2,mean(corVec)/std(corVec)); % divide  by 2 because one-sided

%% Linear regression:

fprintf('Linear regression coefficients: intercept, BOLD main effect, valence main effect, interaction:\n');
for iCoef = 1:size(corMat.b,2)
    
    % Extract b-weights:
    b = corMat.b(validSubs,iCoef);
    fprintf('Regressor no. %d: Min = %.02f, Max = %.02f, Mean = %.02f\n',iCoef,min(b),max(b),mean(b)) % before Fisher-z transform
    b     = atanh(b); % Fisher-z transform
    fprintf('Fisher z-transform\n');

    % T-test:
    fprintf('Perform t-test across subjects: \n');
    [~,P,~,STATS] = ttest(b);
    fprintf('Regressor no. %d: t(%d) = %.02f, p = %.03f, d = %.02f\n',iCoef,STATS.df,STATS.tstat,P/2,mean(b)/std(b)); % divide by 2 because one-sided
    fprintf('--------------\n');
end

%% t-tests on quantiles against each other:

valNames = {'Win','Avoid'};

fprintf('Perform t-tests of conditions (valence, bins) against each other: \n');

for iVal = 1:2
    for iCol = 2:nBin
        fprintf('%s trials, compare bins %d and %d: ',valNames{iVal},iCol-1,iCol);
        difVec          = corMat.RT(validSubs,iCol,iVal)-corMat.RT(validSubs,iCol-1,iVal);
        [~,P,~,STATS]  = ttest(difVec);
        fprintf('t(%d) = %.02f, p = %.03f, d = %.02f\n',STATS.df,STATS.tstat,P/2,mean(difVec)/std(difVec)); % divide  by 2 because one-sided
    end
end

%% Plot two subplots, one per valence, point + errorbar per bin, individual subjecgt data points: 

% Fixed plot settings:
colMat      = [0 .6 .2; .8 0 0];
markerSize  = 8; lineWidthFig = 4; capSize = 8; fontSize = 36; pointSize = 180;
xLim        = [-1.5 1.5];

% Determine y-axis limits, ticks, and label:
if normRT
    yLim   = [-0.59 0.82]; yTick = -2:0.2:2; yLabel = 'RT (z)'; 
else
    yLim   = [0.52 1.01]; yTick = -2:0.1:2; yLabel = 'RT';
end

% Compute means and standard errors per condition (valence x bin):
meanBOLD    = squeeze(mean(corMat.BOLD(validSubs,:,:),1));
seBOLD      = squeeze(std(corMat.BOLD(validSubs,:,:),1)/sqrt(nSubValid));
meanRT      = squeeze(mean(corMat.RT(validSubs,:,:),1));
seRT        = squeeze(std(corMat.RT(validSubs,:,:),1)/sqrt(nSubValid));

% Start plot:
figure('Position',[100 100 1600 800]); hold on
for iVal = 1:2
subplot(1,2,iVal);

    fprintf('\nVal = %d\n',iVal);
    fprintf('Mean BOLD per bin: %s \n',strjoin(string(meanBOLD(:,iVal)),'; '));    
    fprintf('BOLD: min = %.3f; max = %.3f \n',min(corMat.BOLD(:)),max(corMat.BOLD(:)));    
    fprintf('Mean RT   per bin: %s \n',strjoin(string(meanRT(:,iVal)),'; '));    
    fprintf('RT  : min = %.3f; max = %.3f \n',min(corMat.RT(:)),max(corMat.RT(:)));    

    % Add mean per bin:   
    plot(meanBOLD(:,iVal),meanRT(:,iVal),'.','Color',colMat(iVal,:),'MarkerSize',markerSize,'linewidth',2); hold on

    % Vertical error bars (RTs):
    errorbar(meanBOLD(:,iVal),meanRT(:,iVal),seRT(:,iVal),'vertical','Color',colMat(iVal,:),'linestyle','none','linewidth',2,'Capsize',capSize); hold on % error bars

    % Horizontal error bars (BOLD): 
    errorbar(meanBOLD(:,iVal),meanRT(:,iVal),seBOLD(:,iVal),'horizontal','Color',colMat(iVal,:),'linestyle','none','linewidth',2,'Capsize',capSize); hold on % error bars

    % Add individual data points:
    for iBin = 1:nBin
        a = corMat.BOLD(validSubs,iBin,iVal)';
        b = corMat.RT(validSubs,iBin,iVal)';
        s = scatter(a,b,10,'k'); % no jitter
        set(s,'MarkerEdgeColor',colMat(iVal,:),'linewidth',2); % was 1
    end

    % Plot settings:
    set(gca,'xlim',xLim,'ylim',yLim,'xtick',-2:1:2,'ytick',yTick,...
        'fontsize',fontSize,'Linewidth',lineWidthFig) %,'ytick',-5:0.1:4)
    xlabel(sprintf('%s BOLD (z)',name),'fontweight','normal','fontsize',fontSize);
    ylabel(yLabel,'fontweight','normal','fontsize',fontSize); 

end % end iVal

% Save plot:

saveas(gcf,sprintf('/project/3017042.02/Log/CueLockedPaperPlots/RT_%s_split_subplot_nBin%d%s.png',string(ROI2use),nBin,suffix));
% pause(1)
% close gcf
fprintf('Finished plot :-)\n');

%% Plot both valences in one plot, point + errorbar per bin per valence, individual subject data points: 

% Fixed plot settings:
colMat      = [0 .6 .2; .8 0 0];
markerSize  = 16; lineWidthFig = 4; lineWidthLines = 3; capSize = 12; fontSize = 44; pointSize = 250;
fprintf('BOLD: min = %.03f; max = %.03f \n',min(corMat.BOLD(:)),max(corMat.BOLD(:)));    
fprintf('RT: min = %.03f; max = %.03f \n',min(corMat.RT(:)),max(corMat.RT(:)));    
xLim        = [-1.45 1.3]; % Minimal and maximal BOLD

% Determine y-axis limits, ticks, and label:
if normRT % if standardized
    yLim   = [-0.59 0.82]; yTick = -2:0.2:2; yLabel = 'RT (z)'; 
else % minimal and maximal RT for any bin if non-standardized
    yLim   = [0.515 1.015]; yTick = -2:0.1:2; yLabel = 'RT';
end

% Prepare data: average over subjects, keep nBin and nVal dimensions
meanBOLD    = squeeze(mean(corMat.BOLD(validSubs,:,:),1));
seBOLD      = squeeze(std(corMat.BOLD(validSubs,:,:),1)/sqrt(nSubValid));
meanRT      = squeeze(mean(corMat.RT(validSubs,:,:),1));
seRT        = squeeze(std(corMat.RT(validSubs,:,:),1)/sqrt(nSubValid));
rng(19913010) % keep jitter constant

% Start figure:
figure('Position',[100 100 1200 800]); hold on

% Line connections:
for iBin = 1:nBin 
    plot(meanBOLD(iBin,:),meanRT(iBin,:),'k.','linewidth',lineWidthLines); hold on
end

for iVal = 1:2
    
    % Condition means as points:
    plot(meanBOLD(:,iVal),meanRT(:,iVal),'-','Color',colMat(iVal,:),'MarkerSize',markerSize,'linewidth',lineWidthLines); hold on
    
    % Vertical error bar (RT):
    errorbar(meanBOLD(:,iVal),meanRT(:,iVal),seRT(:,iVal),'vertical','Color',[0 0 0],'linestyle','none','linewidth',lineWidthLines,'Capsize',capSize); hold on % error bars
    
    % Horizontal error bar (BOLD):
    errorbar(meanBOLD(:,iVal),meanRT(:,iVal),seBOLD(:,iVal),'horizontal','Color',[0 0 0],'linestyle','none','linewidth',lineWidthLines,'Capsize',capSize); hold on % error bars
    
    % Individual data points:
    for iBin = 1:nBin
        
        % Jitter:
        a = corMat.BOLD(validSubs,iBin,iVal)'; % BOLD values (x-axis)
        b = corMat.RT(validSubs,iBin,iVal)'; % RT values (y-axis)
        s = scatter(a,b,[],'k'); % no jitter
        set(s,'MarkerEdgeColor',colMat(iVal,:),'linewidth',3); % was 1
        s.MarkerEdgeAlpha = .5; % somewhat transparent
        
    end
end

% Plot settings:
set(gca,'xlim',xLim,'ylim',yLim,'xtick',-2:1:2,'ytick',yTick,...
    'fontsize',fontSize,'Linewidth',lineWidthFig) %,'ytick',-5:0.1:4)
xlabel(sprintf('%s BOLD (z)',name),'fontweight','normal','fontsize',fontSize);
ylabel(yLabel,'fontweight','normal','fontsize',32); 

% Save:
if withoutInvalid
    saveas(gcf,fullfile(dirs.root,sprintf('Log/CueLockedPaperPlots/RT_%s_split_oneplotpoints_nBin%d_withoutInvalid.png',string(ROI2use),nBin)));
else
    saveas(gcf,fullfile(dirs.root,sprintf('Log/CueLockedPaperPlots/RT_%s_split_oneplotpoints_nBin%d%s.png',string(ROI2use),nBin,suffix)));
end
% pause(1)
% close gcf

fprintf('Finished plot :-)\n');

end % end of function.
