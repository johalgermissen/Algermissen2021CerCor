function output = regress_baseline(input,x,startIdx,stopIdx,corMode)

% output = regress_baseline(input,x,startIdx,stopIdx,corMode)
%
% takes EEG data, extracts average baseline in specified
% time window, computes regression over trials, corrects for trial-by-trial
% baseline predicted by that regression.
%
% INPUTS:
% input         = matrix of n-dimensional format: 
%               First dimension is trials (do regression over)
%               Last dimension is time (extract baseline period)
%               Dimensions in between don't matter, will be concatenated
% x             = single regressor in regression -- indices of trials
% provided (default: 1 until number of trials provided)
% Note that some trials might be rejected, so x might be discontinuous
% startIdx      = scalar integer, time index (in last dimension) from when 
% onwards to do baseline extraction
% stopIdx       = scalar integer, time index (in last dimension) until when 
% to do baseline extraction
% corMode       = string, type of baseline correction, either 'subtraction'
% or 'division'
%
% OUTPUTS:
% output        = input matrix corrected for baseline
%
% Recommendations:
% Do 'subtraction' type baseline correction for time-domain data
% Do 'division' type baseline correction for frequency-domain data, perform
% decibel-conversion (10*log10(output)) afterwards
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

%% Fill in input settings:

if nargin < 2
    x = 1:size(input,1);
    fprintf('No trial indices for regressor provided---assume 1 until number of input trials\n')
end
if nargin < 3
    startIdx = 1;
    fprintf('No start index for baseline extraction provided---assume %d\n',startIdx)
end
if nargin < 4
    stopIdx = 1;
    fprintf('No stop index for baseline extraction provided---assume %d\n',stopIdx)
end
if nargin < 5
    corMode = 'subtraction';
    fprintf('No mode of baseline correction provided---perform %s\n',corMode)
end

%% Check if input is cell:

if iscell(input)
    back2cell = true;
    fprintf('Input is of type cell, reshape into matrix\n')
    nDim  = length(size(input{1}))+1; % Determine which would be new extra dimension
    input = permute(cat(nDim,input{:}),[nDim 1:(nDim-1)]); % concatenate into final dimension, bring final dimension to start
else
    back2cell = false;
end

%% Assess matrix dimensions:

dimVec  = size(input); % number of dimensions

nTrial  = dimVec(1); % number trials
nTime   = dimVec(end); % number time bins
nRest   = prod(dimVec)/nTrial/nTime; % number other bins

% Check x: must be column vector:
if size(x,1) == 1
    x = x'; % transpose so now column vector
end    

if length(x) ~= nTrial
    error('Input vector x of different length than number of trials (first dimension input data)');
end

%% Reshape into 3D matrix: trial x rest x time

input_3D = reshape(input,nTrial,nRest,nTime);

%% Perform regression over trials:

fprintf('Perform regression-based baseline-correction over %d trials for %d separate bins\n',nTrial,nRest);

% Initialize empty template for baseline:
baseline = nan(nTrial,nRest,nTime);

% Loop through nRest, extract baseline in specified window, fit regression over trials:

for iRest = 1:nRest % iRest = 1;

    % Extract data in specified baseline window, average within time window:
    y = nanmean(input_3D(:,iRest,startIdx:stopIdx),3);  % average in baselineTimings
    
    % Fit regression over trials:
    p = polyfit(x,y,1); % fit baseline as function of trial number
    
    % Predict baseline: 
    f = polyval(p,x); % predict values from regression for x
    
    % Fill into baseline template:
    baseline(:,iRest,:) = repmat(f,1,1,nTime); % save for all trials
    
end % end of iRest

%% Subtract fitted baseline from actual data:

if strcmp(corMode,'subtraction')
    output_3D = input_3D - baseline; % element-wise subtraction
elseif strcmp(corMode,'division')
    output_3D = input_3D ./ baseline; % element-wise division
else
    error('Unknown correction mode');
end

% Reshape back into original format:
output = reshape(output_3D,size(input));

%% Shape back to cell if input was cell:

if back2cell
    fprintf('Input was of type cell, reshape back into cell\n')
    tmp = output; % reassign
    output = arrayfun(@(iTrial) squeeze(tmp(iTrial,:,:)), 1:size(tmp,1),'UniformOutput', 0);
end

fprintf('Finished baseline correction\n');

end % end of function.