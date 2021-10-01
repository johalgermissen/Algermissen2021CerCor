% Figure04.m

% Plots for Figure 4 in manuscript.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/FiguresCueLocked

% Set root directory:
rootdir = figures_set_rootdir(); % '/project/3017042.02';
addpath(fullfile(rootdir,'Analyses/EEG_Scripts/CueLockedAnalyses/TAfT'));

%% Load data:

EEGdomain   = 'TF';

% ROIs2use    = {}; % empty
ROIs2use    = {'GLM1StriatumAction','GLM1CingulateAnteriorAction','GLM1LeftMotorHand','GLM1RightMotorHand','GLM1vmPFCValenceMan'}; % final paper
% With ACC Valence contrast (instead of Action contrast) ROI:
% ROIs2use    = {'GLM1StriatumAction','GLM1CingulateAnteriorValence','GLM1LeftMotorHand','GLM1RightMotorHand','GLM1vmPFCValenceMan'};
% With striatum split into regions coding positive/negative valence effect:
% ROIs2use    = {'GLM1LeftPutamenValence','GLM1MedialCaudateValence','GLM1CingulateAnteriorAction','GLM1LeftMotorHand','GLM1RightMotorHand','GLM1vmPFCValenceMan'};
% With relative displacement (realignment parameters summary metric based on Fellner et al., 2016):
% ROIs2use    = {'GLM1StriatumAction','GLM1CingulateAnteriorAction','GLM1LeftMotorHand','GLM1RightMotorHand','GLM1vmPFCValenceMan','mp_Fellner'};

% 2) Behavioral regressors:
% behav2use   = {}; % none
behav2use   = {'isgo','iswin','isW2G'}; % main effects and interaction --> final paper

% 3) Perform only on selected trials:
selTrials   = 'all';

[job, dirs, betas] = taft_postprocess_load_job(EEGdomain,ROIs2use,behav2use,selTrials);

%% Figure 4 TAfT: TF PLOT (T-VALUES WITH SIGNIFICANT PATCHES HIGHLIGHTED)

job.invalidSubs = [1 11 15 19 21 25]; % invalid co-registrations and spike outliers (> 5)
selChans = {'FCz','Cz'};

% Select ROI:
for iROI = 1:length(job.regNames)
    
    % Create tg and corrp files:
    rng(20190822) % set random number generator for constant p-values
    [sortBetas,~]   = taft_postprocess_TF_selectData(job,betas,iROI); % select data, align channels across subjects
    sortBetas       = taft_postprocess_FisherZ(sortBetas); % Perform Fisher z-transform:
    [tg,corrp]      = taft_postprocess_TF_TFplot(job,dirs,sortBetas,iROI,selChans,2,1000,2,false); % job,dirs,sortBetas,iROI,selChans,thresh,nP,zlim,isNormalize,isSave
    close gcf

    % Settings:
    zlim    = 3;
    pCrit   = 0.06;
    figure('Position',[0 100 1000 800]); hold on
    timeIdx = size(tg,4);
    contourf(sortBetas{1}.time(1:timeIdx),sortBetas{1}.freq,real(squeeze(tg)),50,'linestyle','none'); % 1.3 sec at index 53
    if strcmp(job.lock,'stim')
        set(gca,'xlim',[0 1.3],'ylim',[1 15],'clim',[-1*zlim  1*zlim],...
            'xtick',[0 0.5 0.815 1 1.3 2],'xtickLabel',{'Cue','0.5','AvgRT','1','+','2'},... % yscale log
            'ytick',[2 4 8 12 15],'yscale','lin',...
            'fontsize',32,'Linewidth',3) % -.25 1.3
        plot([0 0],get(gca,'ylim'),':k','LineWidth',3);
        plot([0.815 0.815],get(gca,'ylim'),':k','LineWidth',3);
    else
        set(gca,'xlim',[-1 0.5],'ylim',[1 14.8],'fontsize',32,'clim',[-1*zlim  1*zlim]);
        plot([0 0],get(gca,'ylim'),':k');
    end
    xlabel('Time (in s)','FontSize',32,'FontName','Arial','fontweight','bold');
    ylabel('Frequency (in Hz)','FontSize',32,'FontName','Arial','fontweight','bold');
    plot([0 0],get(gca,'ylim'),':k','LineWidth',3);
    plot([0.815 0.815],get(gca,'ylim'),':k','LineWidth',3);
    colorbar('Ticks',[(-1*zlim):(zlim/2):zlim],'Fontsize',32)
    if sum(corrp(:) < pCrit) > 0 % if anything significant
        % see https://nl.mathworks.com/matlabcentral/answers/250279-how-to-make-contour-plots-transparent-in-matlab-r2015a#answer_211204
       pause(1)
       isSig = double(corrp < pCrit); % map significant clusters
       hold on % on top of old plot
       [~, hContour]  = contourf(sortBetas{1}.time(1:timeIdx),sortBetas{1}.freq,squeeze(isSig),1);
       hContour.LineWidth = 5;
       drawnow;  % this is important, to ensure that FacePrims is ready in the next line!
       hFills = hContour.FacePrims;  % array of TriangleStrip objects
       [hFills.ColorType] = deal('truecoloralpha');  % default = 'truecolor'
       for idx = 1:numel(hFills)
           hFills(idx).ColorData(4) = 1;   % default=255
       end
       hold off
    end
    % Save as:
    fig = gcf;
    saveas(fig,sprintf('/project/3017042.02/Log/CueLockedPaperPlots/TAfT_TFplot_%s_%s_%s_zlim%.1f.png',job.regNames{iROI},job.lock,strjoin(selChans,''),zlim)); % automatic ROI name
%     saveas(fig,sprintf('/project/3017042.02/Log/CueLockedPaperPlots/TAfT_TFplot_%s_withmpFellner_short_%s_%s_zlim%.1f.png',job.regNames{iROI},job.lock,strjoin(selChans,''),zlim)); % automatic ROI name
%     saveas(fig,sprintf('/project/3017042.02/Log/CueLockedPaperPlots/TAfT_TFplot_%s_withmpFellner_long_%s_%s_zlim%d.png',job.regNames{iROI},job.lock,strjoin(selChans,''),zlim)); % manual name
%     saveas(fig,sprintf('/project/3017042.02/Log/CueLockedPaperPlots/TAfT_TFplot_%s_withmpFellner_long_outliers_%s_%s_zlim%d.png',job.regNames{iROI},job.lock,strjoin(selChans,''),zlim)); % manual name
    pause(1)
    close gcf
end

%% Compute sigLine:

job.invalidSubs = [1 11 15 19 21 25]; % invalid co-registrations and spike outliers (> 5)
timeStop    = 53;
sigLine = zeros(size(betas{1}.powspctrm,1),timeStop);

for iROI = [1 5] % iROI = 5;
    selChans        = {'FCz','Cz'};
    rng(20190822) % set random number generator for constant p-values
    [sortBetas,~]   = taft_postprocess_TF_selectData(job,betas,iROI); % select data, align channels across subjects
    sortBetas       = taft_postprocess_FisherZ(sortBetas); % Perform Fisher z-transform:
    [tg,corrp]      = taft_postprocess_TF_TFplot(job,dirs,sortBetas,iROI,selChans,2,1000,2,false); % iROI, selChans, thresh,nP
    close all
    
    % Settings:
    pCrit       = 0.05;

    % a) Sum of columns, normalize later:
    sigLine(iROI,:) = mean(squeeze(tg .* double(corrp<pCrit)),1); % signed power above threshold
    
    % Normalize by peak: 
    sigLine(iROI,:) = sigLine(iROI,:) / max(abs(sigLine(iROI,:))) *100;
    
end
fprintf('Finished :-)\n');

%% Load raw theta signal:

job_tmp     = job; % store job settings from TafT, get back later

% only use 28 subjects:
addpath(fullfile(rootdir,'Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Helpers'));
addpath(fullfile(rootdir,'Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel'));

dirs        = set_dirs(rootdir);
par         = set_par();

job = []; % initialize empty job

job.dirs = dirs; % add directories

% Data settings:
job.nSub                = 36; % necessary for validSubs before loading data
job.sub2exclude         = [11 12 15 23 25 26 30]; % exclusion in TAfT-outcome-locked

job.lockSettings        = 'stimlocked'; % stimlocked resplocked

job.baselineSettings    = 'trend'; % 'all' or 'condition' or 'trial' or 'ft_trial' or 'trend'    
job.baselineTimings     = [-0.250 -0.050];
job.baselineTimings     = 0;

job.responseSettings    = 'Go'; % 'Go' or 'Hand'  
job.actionSettings      = 'exeAct'; % 'reqAct' or 'exeAct'
job.accSettings         = 'correct'; % 'bothAcc' or 'correct' or 'incorrect' or 'noAcc' or 'reqAct'
% for bothAcc, always use reqAct
job.TFtype              = 'hanning'; % hanning morlet
job.nFreq               = 15; % 15
% Channel settings:
job.chanArea            = 'midfrontal'; % 'midfrontal' or 'frontal' or 'central' or 'leftmotor' or 'rightmotor'
% Band settings:
job.band                = 'theta'; % 'theta' or 'alpha' or 'beta' or 'broad' 
% Contrast of interest:
job.contrastType        = 'Go'; % 'Congruency' or 'Go' or 'GoLeftRight' or 'GoValence' or 'GoAccuracy' or 'Accuracy' or 'CongruentAccuracy' or 'IncongruentAccuracy'
% ERP-corrected:
job.stimERPcor = false;
job.respERPcor = false;

% Load and prepare data:
job         = TF_update_job(job); % initialize job
[job, data] = TF_load_data(job); % load data
[job, data] = TF_prepare_generic_data(job,data); % prepare generic objects
job         = TF_update_job(job); % update job
[job, data] = TF_prepare_contrast_data(job,data); % prepare objects specific to contrast

% ---------------------------------------------------- %
% Relevant conditions:
iCond   = [1 2 3 4];

% job.invalidSubs initialized above:
job.invalidSubs = [1 11 15 19 21 25]; % invalid co-registrations and spike outliers (> 5)
job.validSubs = setdiff([1:job.nSub],job.invalidSubs);

% Baseline:
iBaseline = find(data.mu.time == 0);  % iBaseline = 41;
baseline = min(squeeze(nanmean(data.SubCondTime(job.validSubs,iCond,iBaseline))));

% Extract relevant conditions:
time        = data.mu.time;
thetaPower  = squeeze(nanmean(data.SubCondTime(job.validSubs,iCond,:),1))-baseline;
job         = job_tmp; % get TAfT job settings back

%% Plot sigLine:

lineWidth   = 8;
fontSize    = 44;
legendNames = {'vmPFC','Striatum'}; condVec = [5 1];
colMat      = [0.75 0.75 0; 1 0 0; .8 0 0; .75 .56 0; 0 0 1]; % 5 colors
xWidth      = 1800; % 1200 1800
includeEEG  = true;

close all
figure('Position',[0 0 xWidth 800]); hold on % same as TFplot for standard EEG
iCond = 0;

% Plot TAfT lines:
p = {};

% yyaxis left
for iROI = condVec % iROI = 5;
    iCond = iCond + 1;
    p{iCond} = plot(sortBetas{1}.time(1:timeStop),sigLine(iROI,:),'color',colMat(iROI,:),'linewidth',lineWidth);
end

% Plot raw theta lines:
if includeEEG
    % yyaxis right
    maxTheta = max(thetaPower(:));
    job.colMat      = [0 0.6 0.2;0.8 0 0;0 0.6 0.2;0.8 0 0];
    job.lineStyle   = {'-','-','--','--'};
    for iCond = 1:size(thetaPower,1)
        plot(time,thetaPower(iCond,:)/maxTheta*100,'-','linewidth',lineWidth/4*2,'Color',job.colMat(iCond,:),'Linestyle',job.lineStyle{iCond}); % Go2Win: green, solid
    end
end

% Add vertical line at zero:
plot(sortBetas{1}.time,zeros(length(sortBetas{1}.time)),'k-','linewidth',lineWidth);

% Add horizontal lines:
plot([0 0],get(gca,'ylim'),':k','LineWidth',3);
plot([0.81 0.81],get(gca,'ylim'),':k','LineWidth',3);

% Settings:
set(gca,'xlim',[0 1.3],...
    'xtick',[0 0.5 0.81 1 1.3],'xtickLabel',{'Cue','0.5','AvgRT','1','+'},... % 'ytick',[2 4 8 16 32],... 
    'fontsize',fontSize,'Linewidth',lineWidth) % -.25 1.3
xlabel('Time (s)','FontSize',fontSize,'FontName','Arial','fontweight','normal');
ylabel('% of maximum','FontSize',fontSize,'FontName','Arial','fontweight','normal'); 

% Save:
if includeEEG
    saveas(gcf,sprintf('/project/3017042.02/Log/CueLockedPaperPlots/MeanTlines_%s_theta_xWidth%d.png',strjoin(legendNames,'_'),xWidth)); % automatic ROI name
else
    saveas(gcf,sprintf('/project/3017042.02/Log/CueLockedPaperPlots/MeanTlines_%s_thresh_t_xWidth%d.png',strjoin(legendNames,'_'),xWidth)); % automatic ROI name
end
pause(1)
close gcf

%% Figure 4. TAfT: TOPOPLOT TF DATA 

% Select ROI, timing, frequency:
% iROI = 1; iTime = [0.6 1.1]; iFreq = [1 10]; % Striatum
% iROI = 3; iTime = [0.5 1.3]; iFreq = [10 15]; % LeftMotor
% iROI = 4; iTime = [0.6 1.3]; iFreq = [10 15]; % RightMotor
% iROI = 5; iTime = [0.1 0.5]; iFreq = [1 7]; % vmPFC

zlim = 3;

[sortBetas,Tvalues] = taft_postprocess_TF_selectData(job,betas,iROI);

figure('units','normalized','outerposition',[0 0 1 1]); hold on % fullscreen--no risk of implicitly resizing any axes
cfg = []; cfg.ylim = iFreq; cfg.zlim = [-1*zlim 1*zlim]; cfg.marker='on'; cfg.style='straight';
cfg.layout = 'easycapM11.mat'; cfg.comment = 'no'; cfg.xlim = [iTime(1) iTime(end)];
cfg.colorbar = 'yes'; % want 'no', i.e. do it yourself % --> add externally
ft_topoplotTFR(cfg,Tvalues);
title(sprintf('%s, %.2f to %.2f sec', job.ROIs(iROI).ROIname, iTime(1), iTime(end)),'fontsize',32)

% Save:
saveas(gcf,sprintf('/project/3017042.02/Log/CueLockedPaperPlots/TAfT_Topoplot_%s_%s_%.1f-%.1fsec_%d-%dHz_zlim%.1f.png',job.ROIs(iROI).ROIname,job.lock,iTime(1),iTime(end),iFreq(1),iFreq(end),zlim)); % automatic ROI name
close gcf

% END
