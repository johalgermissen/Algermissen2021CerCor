% Figure02.m

% Plots for Figure 2 in manuscript.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/FiguresCueLocked

% Set root directory:
rootdir = figures_set_rootdir(); % '/project/3017042.02';

%% Figure 2B-D. SLICE DISPLAY MAPS:

addpath(fullfile(rootdir,'Analyses/fMRI_Scripts/Display_Results/'));

dirs.root           = rootdir;

dirs.save           = fullfile(dirs.root,'Log/CueLockedPaperPlots');
job.lockSettings    = 'stimlocked';
job.type            = 'standard';
job.zLim            = [0 3]; % 
job.GLMID           = '1';
% job.GLMID           = '1without6';

% Sagittal:
job.iCope = 1; job.iView = 'sagittal'; job.iSlice = 2; job.cLim = 100; % 
EEGfMRIPav_sliceDisplay_cope(dirs,job)

job.iCope = 2; job.iView = 'sagittal'; job.iSlice = 2; job.cLim = 100; % 
EEGfMRIPav_sliceDisplay_cope(dirs,job)

job.iCope = 3; job.iView = 'sagittal'; job.iSlice = 2; job.cLim = 30; % 
EEGfMRIPav_sliceDisplay_cope(dirs,job)

% Coronal:
job.postFix = '_striatum';
job.iCope = 1; job.iView = 'coronal'; job.iSlice = -10; job.cLim = 100; % 
EEGfMRIPav_sliceDisplay_cope(dirs,job)
job = rmfield(job,'postFix');

job.postFix = '_striatum';
job.iCope = 1; job.iView = 'coronal'; job.iSlice = 4; job.cLim = 100; % 
EEGfMRIPav_sliceDisplay_cope(dirs,job)
job = rmfield(job,'postFix');

job.postFix = '_ACC_JBL';
job.iCope = 3; job.iView = 'coronal'; job.iSlice = 4; job.cLim = 30; % 
EEGfMRIPav_sliceDisplay_cope(dirs,job)
job = rmfield(job,'postFix');

%% Figure 2A, E-G. fMRI PARAMETER ESTIMATES:

dirs.root   = rootdir; 

% GLM selection:
GLMID       = '1';
dirs.data   = fullfile(dirs.root,sprintf('Log/fMRI/fMRI_ROIs/GLM%s_ROIData',GLMID));
points      = true;

% Subject selection:
% 15 and 25 already excluded, so now 1 11 18 (actually 19), 20 (actually 21), 24 (actually 26)
nSub            = 34;
invalidSubs     = []; % subject already excluded in input file
% invalidSubs     = [1 11 19 21]; % outliers in TAfT and fMRI (15 and 25 already excluded)
fprintf('Exclude subjects %s\n',num2str(invalidSubs));
validSubs       = setdiff(1:nSub,invalidSubs);
nSubValid       = length(validSubs);

% Create hypothesis data:
data        = [0.08;-0.04;0.04;-0.08];
dlmwrite(fullfile(dirs.data,'Hypothesis.txt'),data);

% Load data:
ROInames    = {'GLM1LeftPutamenValence','GLM1MedialCaudateValence','GLM1JBLIncongruency'};
ROItitles   = {'Left Putamen (Valence)','Medial Caudate (Valence)','Pre-SMA (Incongruency)'};

for iROI = 1:length(ROInames)
    ROIname     = ROInames{iROI};
    ROItitle    = ROItitles{iROI};
    data        = load(fullfile(dirs.data,sprintf('%s.txt',ROIname)));
    
    % Compute standard errors:
    nCond       = size(data,1);
    if size(data,2)==1; validSubs = 1; end % for Hypothesis field
    subMean     = nanmean(data(:,validSubs),1); % average over conditions --> one value per subject
    grandMean   = nanmean(subMean); % average over subjects --> grand mean (1 value)
    
    % Correct SE by subtracting overall mean per subject and adding back grand average; correct for # conditions
    condMean    = nan(nCond,1);
    condSE      = nan(nCond,1);
    for iCond = 1:nCond
        condMean(iCond) = nanmean(data(iCond,validSubs));
        condSE(iCond)   = nCond/(nCond-1)*nanstd(data(iCond,validSubs)-subMean+repmat(grandMean,1,nSubValid))./sqrt(nSubValid);
    end
    
    % General Settings:
    xLoc = [0.5 1.5 3.0 4.0];
    ColMat = [0 .6 .2; .8 0 0;0 .6 .2; .8 0 0];
    lineWidth = 4; capSize = 12; fontSize = 36;
    fprintf('ROI %s: min = %.03f, max = %.03f\n',ROIname,min(data(:)),max(data(:)));
    if strcmp(ROIname,'Hypothesis')
        yTick = [];
    else
        yTick = -1:.05:1;
    end
    if strcmp(ROIname,'GLM1LeftPutamenValence') || strcmp(ROIname,'GLM1MedialCaudateValence')
        yLim = [-0.32 0.32];
        yTick = -1:.2:1;
    elseif strcmp(ROIname,'GLM1JBLIncongruency')
        yLim = [-0.2 0.75];
        yTick = -1.2:.4:1.2;
    else
        yLim = [-0.1 0.1];
    end
    if points && ~isempty(invalidSubs)
        yTick = [-1:.2:1];
        yLim = [round(min(data(:)),2)-.01 round(max(data(:)),2)+.01];

    end
    % Figure:
    figure('Position',[100 100 1000 800]); hold on
    % Bar plots:
    for iCond = 1:nCond
        p{iCond} = bar(xLoc(iCond),condMean(iCond),.75,'FaceColor',ColMat(iCond,:)); % bar plot
        if points
            s = scatter(repmat(xLoc(iCond),1,nSubValid),data(iCond,validSubs),[],'k', 'jitter','on', 'jitterAmount',0.15); hold on % was 0.05
            set(s,'MarkerEdgeColor',[0.4 0.4 0.4],'linewidth',3); % was 1 
        end
    end
    % Error bars:
    if ~strcmp(ROIname,'Hypothesis')
        for iCond = 1:nCond
            errorbar(xLoc(iCond),condMean(iCond),condSE(iCond),'k','linestyle','none','linewidth',lineWidth,'Capsize',capSize); % error bars
        end
    end
    % Settings:
    set(gca,'FontName','Arial','FontSize',fontSize,'Linewidth',lineWidth);
    set(gca,'xlim',[0 max(xLoc)+0.6],'xtick',[1 3.5],'xticklabel',{'Go','NoGo'})
    xlabel('Performed action','FontSize',fontSize,'FontName','Arial','Color',[0 0 0]);
    set(gca,'ylim',yLim)
    set(gca,'ytick',yTick)
    ylabel('Param. estimates (a.u.)','FontSize',fontSize,'FontName','Arial','Color',[0 0 0]);
    title(sprintf('%s',ROItitle),'FontSize',fontSize,'FontName','Arial','Color',[0 0 0]);
    legend([p{1} p{2}],{'Win','Avoid'},'fontsize',fontSize); legend boxoff;
    % Save:
    outputfile = sprintf('/project/3017042.02/Log/CueLockedPaperPlots/Matlab_Barplot_%s',ROIname);
    if ~isempty(invalidSubs)
        outputfile  = sprintf('%s_withoutInvalid',outputfile);
    end
    if points
        outputfile  = sprintf('%s_points',outputfile);
    end
    saveas(gcf,sprintf('%s.png',outputfile));
    close gcf
end

%% Figure 2H. fMRI-RT correlation plots

addpath(fullfile(rootdir,'Analyses/EEG_Scripts/CueLockedAnalyses/TAfT'));

taft_correlation_regression_BOLD_RT_per_valence_split({'GLM1vmPFCValenceMan'},'vmPFC')
taft_correlation_regression_BOLD_RT_per_valence_split({'GLM1CingulateAnteriorValence'},'ACC')
taft_correlation_regression_BOLD_RT_per_valence_split({'GLM1vmPFCValenceMan'},'vmPFC',true)
taft_correlation_regression_BOLD_RT_per_valence_split({'GLM1CingulateAnteriorValence'},'ACC',true)

% END.
