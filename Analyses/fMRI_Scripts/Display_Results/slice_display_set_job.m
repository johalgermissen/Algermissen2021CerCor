function job = slice_display_set_job(job)

% job = slice_display_set_job(job)
% 
% Update job settings for slice display.
% 
% INPUT:
% 
% job               = structure with options (necessary and optional settings)
%
% NECESSARY INPUTS:
% .lockSettings     = string, either 'stimlocked' or 'outcomelocked'
% .type             = string, either 'standard' or 'conjunction' (this version only supports 'standard')
% .GLMID            = string, name of GLM
% .iCope            = scalar integer, number of cope under type 'standard'
% .sign             = string, either 'pos' or 'neg'
% .iView            = string, either 'sagittal' or 'coronal' or 'axial'
% .iSlice           = scalar integer, value of slice to display
%
% OPTIONAL INPUTS:
% .levels           = scalar integer, number of levels in GLM (depends on
% lockSettings)
% .cLim             = vector of two floats, c-value range (color intensity)
% .zLim             = vector of two floats, z-value range (opacity)
% .isSave           = Boolean, save plot or not
% 
% OUTPUTS:
% job               = same object with settings filled in.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

%% 1) Necessary settings:

if ~isfield(job,'lockSettings')
    error('No lockSettings specified')
end

if ~isfield(job,'type')
    error('No type specified')
end

if ~isfield(job,'GLMID')
    error('No GLM specified')
end

if ~isfield(job,'iCope')
    error('No iCope specified')
end

if ~isfield(job,'sign')
    error('No sign specified')
end

if ~isfield(job,'iView')
    error('No iView specified')
end

if ~isfield(job,'iSlice')
    error('No iSlice specified')
end

%% 2) Optional settings: contrast

% ------------------------------------------ %
% Number of levels:
if ~isfield(job,'levels')
    job.levels     = 2; % 
end

%% Optional settings: c and z range

% ------------------------------------------ %
% c range:
if ~isfield(job,'cLim')
    job.zLim    = [-100 100]; % 
    fprintf('No zLim specified, set default to %d-%d\n',job.zLim(1),job.zLim(2))
end

% ------------------------------------------ %
% z range:
if ~isfield(job,'zLim')
    job.zLim    = [0 3]; % 
    fprintf('No zLim specified, set default to %0.1f-%0.1f\n',job.zLim(1),job.zLim(2))
end

%% Save or not:

if ~isfield(job,'isSave')
    job.isSave = true;
    fprintf('No isSave specified, set default to %s\n',job.isSave)  
end

end % END of function.
