function plot_PEs_ROI_cue_run()

% plot_PEs_ROI_cue_run()
%
% Calls plot_PEs_ROI_cue() to create bar plots of selected ROIs.
%
% INPUTS:
% none
%
% OUTPUTS:
% none, calls plot_PEs_ROI_cue() which creates and saves plots under 
% 'Log/fMRI/fMRI_ROIs/GLM%s_ROIPlots'
%
% Mind to change rootdir to your own directory structure
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

%% Add path to directory:

fprintf('Add path to Plot_PEs\n');
rootdir = '/project/3017042.02'; 
addpath(fullfile(rootdir,'/Analyses/fMRI_Scripts/Run_ROI/Plot_PEs'));

%% Create bar plots with scatter points for selected ROIs:

plot_PEs_ROI_cue('GLM1LeftPutamenValence','Left Putamen (Valence)',true,[-0.32 0.32]);
plot_PEs_ROI_cue('GLM1MedialCaudateValence','Medial Caudate (Valence)',true,[-0.32 0.32]);
plot_PEs_ROI_cue('GLM1JBLIncongruency','Pre-SMA (Incongruency)',true,[-0.2 0.75]);

end % end of function.
