function Y = z_standardize(X)

% Y = z_standardize(X)
% 
% z-standardize a numeric vector X that might presence of NaNs.
% Based on https://nl.mathworks.com/matlabcentral/answers/105736-zscore-of-an-array-with-nan-s
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/CueLockedAnalyses/TAfT/

if any(isnan(X(:))) % if any NaNs
    
    xmu=nanmean(X); % extract mean without NAs
    xsigma=nanstd(X); % extra std without NAs
    Y=(X-repmat(xmu,length(X),1))./repmat(xsigma,length(X),1); % repeat into vector, object-wise division
else
    
    [Y,xmu,xsigma]=zscore(X); % otherwise just use standard z-score function

end

end % end of function.
