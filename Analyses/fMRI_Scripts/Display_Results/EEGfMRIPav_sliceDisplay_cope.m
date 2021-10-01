function EEGfMRIPav_sliceDisplay_cope(dirs,job)

% EEGfMRIPav_sliceDisplay_cope(dirs,job)
%
% Plots slice of given contrast for given GLM using slice display toolbox.
%
% Before using: mind to set directories/paths in slice_display_set_dirs to respective toolboxes needed!
%
% INPUTS:
% dirs              = structure with various directories (none required,
% defaults set automatically)
% job               = structure with various settings,
% .lockSettings     = string, either 'stimlocked' or 'outcomelocked' (this version only supports 'stimlocked')
% .type             = string, either 'standard' or 'conjunction' (this version only supports 'standard')
% .GLMID            = string, name of GLM
% .sign             = string, either 'pos' or 'neg'
% .iCope            = scalar integer, number of cope to plot
% .iView            = string, either 'sagittal' or 'coronal' or 'axial'
% .iSlice           = scalar integer, value of slice to display
% .cLim             = vector of two floats, c-value range (color intensity)
% .zLim             = vector of two floats, z-value range (opacity)
% .isSave           = Boolean, save plot or not
%
% OUTPUTS:
% generate plot, save if specified.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

%% 1) Set directory and job:

dirs    = slice_display_set_dirs(dirs);
job     = slice_display_set_job(job);

% Load labels (based on GLMID):
tmp = load(fullfile(dirs.label,sprintf('GLM%slabels.mat',job.GLMID))); % plot and axis labels
GLMlabels = tmp.GLMlabels;

% Load color maps:
fprintf('Load slice display color maps\n');
load(fullfile(dirs.sliceDisplay,'colormaps.mat')); % Get custom colormaps from slice display

%% 2) Determine data directory:

% Directory with cope/zstat/thresh files:
dirs.data = fullfile(dirs.root,sprintf('Log/fMRI/GLM%s_FEAT_Combined_all',job.GLMID));

% if job.levels == 2
% 
% elseif job.levels == 3
%     dirs.data  = fullfile(dirs.root,sprintf('Log/fMRI/GLM%s_FEAT_Combined',job.GLMID));
% end

%% 3) Create plot:

fprintf('Create plot for GLM %s, contrast %d, %s view, slice %d \n',job.GLMID, job.iCope, job.iView, job.iSlice);
job.sign = GLMlabels{job.iCope,4}; % retrieve sign of contrast: pos or neg
            
% Step 1: Initialize empty layers and settings variables:
% ------------------------------------------------------------------------------
layers                              = sd_config_layers('init',{'truecolor','dual','contour'});
settings                            = sd_config_settings('init');

% Step 2: Define layers
% ------------------------------------------------------------------------------

% Layer 1: Anatomical map:
% layers(1).color.file                = fullfile(spm('Dir'),'canonical','single_subj_T1.nii');
layers(1).color.file                = fullfile(dirs.sliceDisplay,'MNI152_T1_2mm_brain.nii');
layers(1).color.map                 = gray(256);

% Layer 2: Dual-coded layer (contrast estimates color-coded; z-statistics opacity-coded):
layers(2).color.file                = fullfile(dirs.data,sprintf('cope1_cope%d_%s.nii',job.iCope,job.sign)); % beta-map input: cope1.nii
layers(2).color.map                 = CyBuGyRdYl; % one of the custom color-maps in slice_display
layers(2).color.label               = GLMlabels{job.iCope,3}; % '\beta_{left} - \beta_{right} (a.u.)'; % title for legend
layers(2).color.range               = [-1*job.cLim job.cLim]; % y-axis range (opacity) of legend (z/t-values)

% Layer 2 settings:
layers(2).opacity.file              = fullfile(dirs.data,sprintf('zstat1_cope%d_%s.nii',job.iCope,job.sign)); % t-map/ z-map input: zstat1.nii
layers(2).opacity.label             = '| z |'; % y-axis label (opacity) of legend
layers(2).opacity.range             = job.zLim; % y-axis range (opacity) of legend (z/t-values)

% Layer 3: Contour of significantly activated voxels (thresholded map: thresh_zstat1.nii.gz)
if ~isfield(job,'postFix')
    job.postFix = '';
else
    fprintf('Add postfix %s to thresholded map\n',job.postFix); 
end
layers(3).color.file                = fullfile(dirs.data,sprintf('thresh_zstat1_cope%d_all%s.nii',job.iCope,job.postFix));
layers(3).color.map                 = [0 0 0]; % gray is [0.3 0.3 0.3]; # black is [0 0 0];
layers(3).color.line_width          = 2; % default 3
layers(3).color.line_style          = '-';

% Image settings:
settings.fig_specs.n.slice_column   = 1; % # columms for display
settings.fig_specs.title            = GLMlabels{job.iCope,2}; % 'left - right'; % caption
settings.slice.orientation          = char(job.iView); % axial sagittal coronal
settings.slice.disp_slices          = job.iSlice; % single slice

% Display the layers:
[settings,p] = sd_display(layers,settings);
frame_h = get(handle(gcf),'JavaFrame');
% set(frame_h,'Maximized',1);
% set(frame_h,[500 800 10 600],1);
            
% Save figure:
fileName = sprintf('SD_GLM%s_Cope%02d_c%d_z%d-%d_%s_slice%d',job.GLMID,job.iCope,...
    layers(2).color.range(2),layers(2).opacity.range(1),layers(2).opacity.range(2),...
    settings.slice.orientation,int8(settings.slice.disp_slices(1)));
saveas(gcf,fullfile(dirs.save,[fileName '.png']));
saveas(gcf,fullfile(dirs.save,[fileName '.eps']),'epsc');
pause(1)
close gcf

% END
