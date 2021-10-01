% FigureS15.m

% Plots for Figure S15 in manuscript.
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

% For motor cortex plots:
ROIs2use    = {'GLM1StriatumAction','GLM1CingulateAnteriorAction','GLM1LeftMotorHand','GLM1RightMotorHand','GLM1vmPFCValenceMan'}; % final paper

% For results with left putamen and medial caudate instead of striatum:
% ROIs2use    = {'GLM1LeftPutamenValence','GLM1MedialCaudateValence','GLM1CingulateAnteriorAction','GLM1LeftMotorHand','GLM1RightMotorHand','GLM1vmPFCValenceMan'};

% 2) Behavioral regressors:
behav2use   = {'isgo','iswin','isW2G'}; % main effects and interaction --> final paper

% 3) Perform only on selected trials:
selTrials   = 'all';

[job, dirs, betas] = taft_postprocess_load_job(EEGdomain,ROIs2use,behav2use,selTrials);

%% Figure S15 TAfT: TF PLOT (T-VALUES WITH SIGNIFICANT PATCHES HIGHLIGHTED)

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

%% Figure S15. TAfT: TOPOPLOT TF DATA 

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
