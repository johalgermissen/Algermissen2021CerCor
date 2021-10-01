function output = withinSE(input)

% output = withinSE(input)
%
% Compute SEs following Cousineau (2005) and Morey (2008) over subjects
% separately for conditions and any other dimensions; correct for number of
% conditions at the end.
% 
% INPUT: 
% input 	= matrix of n (>= 2) dimensions;
%   first dimension assumed to be 'subject';
%   second dimension assumed to be 'condition'.
% For all other dimensions, separate SEs will be given
%
% OUTPUT:
% output    = matrix of n-1 dimensions with corrected standard errors, 
%   first dimension is 'condition';
%   all other dimensions maintained.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

%% Retrieve dimensions:

dim         = size(input);
if length(dim) < 2
    error('Input must be at least two dimensions (subjects x conditions)')
end

nSub        = dim(1);
nCond       = dim(2);
nRest       = prod(dim)/nSub/nCond;

%% Reshape into nSub x nCond x nRest:

input       = reshape(input,nSub,nCond,nRest);

%% Compute mean per subject (average over conditions, keep rest):

subMean     = nan(size(input));
for iSub = 1:nSub
    subMean(iSub,:,:) = nanmean(input(iSub,:,:),2);
end

%% Compute grand mean per sample (average over conditions, keep rest):

grandMean   = nanmean(subMean,1);

%% Substitute subject mean (per subject over conditions) by grand mean (over subjects and conditions):

output      = input - subMean + repmat(grandMean,nSub,1,1);

%% Compute standard-deviation (over subjects):

output      = nanstd(output,1); % standard deviation across first dimension (subjects)

%% Divide by sqrt(N) to get SE:

output      = output/sqrt(nSub);

%% Correct for number of conditions (see Morey, 2008):

output      = nCond/(nCond-1)*output;

end % end of function.