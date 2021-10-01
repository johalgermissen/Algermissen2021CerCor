function evaluate_stat(stat,pCrit)

% [p] = evaluate_stat(stat,pCrit)
%
% Locate clusters above threshold for given p-value threshold.
% 
% INPUTS:
% stat          = Fieldtrip object, output of ft_timelockstatistics or
% ft_freqstatistics.
% pCrit         = numeric scalar, critical p-value theshold for clusters
% (default: 0.05).
%
% OUTPUTS:
% None, print to console.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel/

if nargin < 2
    pCrit = 0.05;
end

%% Overview all p-values of any cluster:

if ~isempty(stat.posclusters); fprintf('P-values of positive clusters:\n'); disp([stat.posclusters.prob]); end 
if ~isempty(stat.negclusters); fprintf('P-values of negative clusters:\n'); disp([stat.negclusters.prob]); end 

%% Complete overview:

if isfield(stat,'freq') % if TF input
    % Positive clusters:
    if ~isempty(stat.posclusters)
        for iCluster = 1:length(stat.posclusters) % loop through clusters
            if stat.posclusters(1,iCluster).prob < pCrit % determine if cluster below critical threshold
                if size(stat.posclusterslabelmat,1) > 1 && size(stat.posclusterslabelmat,2) > 1 % if 3D: channels x frequencies x time
                    sigChan = stat.label(squeeze(any(any(stat.posclusterslabelmat==iCluster,2),3))); % average over time and frequencies
                    sigFreq = stat.freq(squeeze(any(any(stat.posclusterslabelmat==iCluster,1),3))); % average over time and channels
                    sigTime = stat.time(squeeze(any(any(stat.posclusterslabelmat==iCluster,1),2))); % average over frequencies and channels
                    fprintf('Positive cluster %d: p = %.03f, %.03f - %.03f seconds, %d - %d Hz, channels %s\n',...
                        iCluster,stat.posclusters(1,iCluster).prob,sigTime(1),sigTime(end),sigFreq(1),sigFreq(end),strjoin(string(sort(sigChan)),'/'));
                elseif size(stat.posclusterslabelmat,1) > 1 % if 2D: channels x time
                    sigChan = stat.label(any(squeeze(stat.posclusterslabelmat==iCluster),2)); % average over time
                    sigTime = stat.time(any(squeeze(stat.posclusterslabelmat==iCluster),1)); % average over channels
                    fprintf('Positive cluster %d: p = %.03f, %.03f - %.03f seconds, channels %s\n',...
                        iCluster,stat.posclusters(1,iCluster).prob,sigTime(1),sigTime(end),strjoin(string(sort(sigChan)),'/'));
                elseif size(stat.posclusterslabelmat,2) > 1 % if 2D: frequencies x time
                    sigFreq = stat.freq(any(squeeze(stat.posclusterslabelmat==iCluster),2)); % average over time
                    sigTime = stat.time(any(squeeze(stat.posclusterslabelmat==iCluster),1)); % average over frequencies
                    fprintf('Positive cluster %d: p = %.03f, %.03f - %.03f seconds, %d - %d Hz\n',...
                        iCluster,stat.posclusters(1,iCluster).prob,sigTime(1),sigTime(end),sigFreq(1),sigFreq(end));
                else % if only time
                    sigTime = stat.time(squeeze(any(stat.posclusterslabelmat==iCluster,1))); % average over frequencies
                    fprintf('Positive cluster %d: p = %.03f, %.03f - %.03f seconds\n',...
                        iCluster,stat.posclusters(1,iCluster).prob,sigTime(1),sigTime(end));
                end
            end
        end
    end
    % ------------------------------------------------------------------- %
    % Negative clusters:
    if ~isempty(stat.negclusters)
        for iCluster = 1:length(stat.negclusters) % loop through clusters
            if stat.negclusters(1,iCluster).prob < pCrit % determine if cluster below critical threshold
                if size(stat.negclusterslabelmat,1) > 1 && size(stat.negclusterslabelmat,2) > 1 % if 3D: channels x frequencies x time
                    sigChan = stat.label(squeeze(any(any(stat.negclusterslabelmat==iCluster,2),3))); % average over time and frequencies
                    sigFreq = stat.freq(squeeze(any(any(stat.negclusterslabelmat==iCluster,1),3))); % average over time and channels
                    sigTime = stat.time(squeeze(any(any(stat.negclusterslabelmat==iCluster,1),2))); % average over frequencies and channels
                    fprintf('Negative cluster %d: p = %.03f, %.03f - %.03f seconds, %d - %d Hz, channels %s\n',...
                        iCluster,stat.negclusters(1,iCluster).prob,sigTime(1),sigTime(end),sigFreq(1),sigFreq(end),strjoin(string(sort(sigChan)),'/'));
                elseif size(stat.negclusterslabelmat,1) > 1 % if 2D: channels x time
                    sigChan = stat.label(any(squeeze(stat.negclusterslabelmat==iCluster),2)); % average over time
                    sigTime = stat.time(any(squeeze(stat.negclusterslabelmat==iCluster),1)); % average over channels
                    fprintf('Negative cluster %d: p = %.03f, %.03f - %.03f seconds, channels %s\n',...
                        iCluster,stat.negclusters(1,iCluster).prob,sigTime(1),sigTime(end),strjoin(string(sort(sigChan)),'/'));
                elseif size(stat.negclusterslabelmat,2) > 1 % if 2D: frequencies x time
                    sigFreq = stat.freq(any(squeeze(stat.negclusterslabelmat==iCluster),2)); % average over time
                    sigTime = stat.time(any(squeeze(stat.negclusterslabelmat==iCluster),1)); % average over frequencies
                    fprintf('Negative cluster %d: p = %.03f, %.03f - %.03f seconds, %d - %d Hz\n',...
                        iCluster,stat.negclusters(1,iCluster).prob,sigTime(1),sigTime(end),sigFreq(1),sigFreq(end));
                else % if only time
                    sigTime = stat.time(squeeze(any(stat.negclusterslabelmat==iCluster,1))); % average over frequencies
                    fprintf('Negative cluster %d: p = %.03f, %.03f - %.03f seconds\n',...
                        iCluster,stat.negclusters(1,iCluster).prob,sigTime(1),sigTime(end));
                end
            end
        end
    end
% ----------------------------------------------------------------------- %
% Time domain:
else % if time input
    % Positive clusters:
    if ~isempty(stat.posclusters) % if any positive clusters available
        for iCluster = 1:length(stat.posclusters) % loop through clusters
            if stat.posclusters(1,iCluster).prob < pCrit % determine if cluster below critical threshold
                if size(stat.posclusterslabelmat,1) > 1 % if multiple channels
                    sigTime = stat.time(any(stat.posclusterslabelmat==iCluster));
                    sigChan = stat.label(any(stat.posclusterslabelmat==iCluster,2));
                    fprintf('Positive cluster %d: %.03f - %.03f seconds, channels %s, p = %.04f\n',...
                        iCluster,sigTime(1),sigTime(end),strjoin(string(sort(sigChan)),'/'),stat.posclusters(1,iCluster).prob);
                else % if averaged over channels
                    sigTime = stat.time(stat.posclusterslabelmat==iCluster); 
                    fprintf('Positive cluster %d: %.03f - %.03f seconds, p = %.04f\n',...
                        iCluster,sigTime(1),sigTime(end),stat.posclusters(1,iCluster).prob);
                end
            end
        end
    end
    % ------------------------------------------------------------------- %
    % Negative clusters:
    if ~isempty(stat.negclusters) % if any negative clusters available
        for iCluster = 1:length(stat.negclusters) % loop through clusters
            if stat.negclusters(1,iCluster).prob < pCrit % determine if cluster below critical threshold
                if size(stat.negclusterslabelmat,1) > 1 % if multiple channels
                    sigTime = stat.time(any(stat.negclusterslabelmat==iCluster)); 
                    sigChan = stat.label(any(stat.negclusterslabelmat==iCluster,2));
                    fprintf('Negative cluster %d: %.03f - %.03f seconds, channels %s, p = %.04f\n',...
                        iCluster,sigTime(1),sigTime(end),strjoin(string(sort(sigChan)),'/'),stat.posclusters(1,iCluster).prob);
                else % if averaged over channels
                    sigTime = stat.time(stat.negclusterslabelmat==iCluster); 
                    fprintf('Negative cluster %d: %.03f - %.03f seconds, p = %.04f\n',...
                        iCluster,sigTime(1),sigTime(end),stat.negclusters(1,iCluster).prob);
                end
            end
        end
    end
end

end % end of function.