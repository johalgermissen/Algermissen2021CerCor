function dirs = slice_display_set_dirs(dirs)

% dirs = slice_display_set_dirs(dirs)
% 
% Initialize directories, add paths to necessary toolboxes.
% All diretories have defaults (none necessary).
% 
% INPUT:
% 
% dirs              = structure with options (optional settings)
%
% OPTIONAL INPUTS:
% .root             = string, root directory of project
% .save             = string, where to save plots
% .panel            = string, path of Panel toolbox version
% .sliceDisplay     = string, path to slice display toolbox version
% .spm              = string, path of SPM version
% 
% OUTPUTS:
% dirs              = same object with settings filled in.
%
% Mind to set dirs.root, dirs.panel, dirs.sliceDisplay, and dirs.spm to your personal settings!
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

%% Root and save directories:

% -------------------------------------------------- %
% Set root directory:
if ~isfield(dirs,'root')
    fprintf('No root directory specified, set default\n');
    dirs.root = '/project/3017042.02'; % root directory--needs to be adapted to users' folder structure
end

% -------------------------------------------------- %
% Set directory where labels are stored:

if ~isfield(dirs,'labels')
    fprintf('No label directory specified, set to default\n');
    dirs.label = fullfile(dirs.root,'Analyses/fMRI_Scripts/Display_Results/Labels/');
end

% -------------------------------------------------- %
% Set directory where to save plots:

if ~isfield(dirs,'save')
    fprintf('No save directory specified, set to default\n');
    dirs.save = fullfile(dirs.root,'Log/CueLockedPaperPlots');
    if ~exist(dirs.save,'dir'); mkdir(dirs.save); end
end

%% 2) Add necessary toolboxes:

% -------------------------------------------------- %
% Add your own panel version (https://www.mathworks.com/matlabcentral/fileexchange/20003-panel):
if ~isfield(dirs,'panel')
    dirs.panel = '/home/action/johalg/matlabScripts/panel-2.12';
end

fprintf('Add paths to panel toolbox\n');
addpath(dirs.panel); 

% -------------------------------------------------- %
% Add your own slice display directory:
if ~isfield(dirs,'sliceDisplay')
    fprintf('No path to slice display specified, set default\n');
    dirs.sliceDisplay = '/home/action/johalg/matlabScripts/slice_display-master'; 
end

fprintf('Add paths to slice display\n');
addpath(dirs.sliceDisplay);

% -------------------------------------------------- %
% Add your own version of SPM:
if ~isfield(dirs,'spm')
    dirs.spm = '/home/common/matlab/spm8';
end

fprintf('Add paths to SPM\n');
addpath(dirs.spm)

end % END of function.
