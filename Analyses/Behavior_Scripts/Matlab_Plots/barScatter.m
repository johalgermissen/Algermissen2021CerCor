function [hb,he]=barScatter(data, errors, bw_xtick, bScatter, bSED, bColor, xpos)
% [hb,he]=barScatter(data, errors, bw_xtick, bScatter)
%
% barScatter takes either an m-by-n matrix of data to be plotted 
% errorbars can be custom specified, or are computed as the SEM across n
% bScatter (true/false) indexes whether scatter points for each subject should be included
% bSEM (true/false) indexes whether error bars should be included
% bColor is a m-by-3 matrix with the colors for the bars in the rows
% xpos is a vector of length m with the positions of the bars on the x-axis
% barScatter calls
% the MATLAB bar function and plots
% - m bars
% - a scatterplot of the individual datapoints, if bScatter = true
% - errorbars of length error (vector of m long), or, if empty, using SD
% - xtick labels on the x axis. 
% NB note that the grouping occurs by row, so each row is a separate group
%
% See the MATLAB functions bar and errorbar for more information
%
% [hb,he]=barScatter(...) will give handle HB to the bars and HE to the
% errorbars.
% 
% ------------------------------------------------------------------------
% Written by Hanneke den Ouden 2015 <h.denouden@gmail.com>
% Donders Center for Cognitive Neuroimaging
% Donders Center for Brain, Cognition and Behavior
% Radboud University Nijmegen

% Adapted by Johannes Algermissen to include color coding option and x-axis
% position for bars and thicker dots
% 2018 <johannes.algermissen@gmail.com>
% Donders Center for Cognition
% Donders Center for Brain, Cognition and Behavior
% Radboud University Nijmegen
% ------------------------------------------------------------------------

% Get function arguments:
% Each function argument requires that all previous arguments have been specified
% All remaining arguments are here set to defaults:
if nargin < 1
    error('Must have at least the first argument:  barweb(data, width, groupnames, bw_title, bw_xlabel, bw_ylabel, bw_colormap, gridstatus, bw_legend)');
end
if nargin <3
    bw_xtick = [];end
if nargin < 4 % include scatters
    bScatter = true;end
if nargin < 5 % include error bars
    bSED = true;end
if nargin < 6 % if no color specified: use Matlab default color (dark blue; find out with get(gca,'ColorOrder'))
    bColor = vec2mat(repelem([0 0.4470 0.7410],size(data,1)),3);end
if nargin < 7 % if no x-axis position specified: just numerate from 1 to number of groups
    xpos = 1:size(data,1);end

barvalues = squeeze(nanmean(data,2)); % mean per group
nbars = size(data,1); % determine number of groups

% Plot bars
 for t = 1:nbars
    hb=bar(xpos(t),barvalues(t,:),'FaceColor',bColor(t,:));
    hold on;
end

% Calculate errors
set(gca, 'fontsize', 10, 'fontweight', 'bold');
he=[];
if ~exist('errors','var')|| isempty(errors) || nargin==1;
    for t = 1:nbars % for each group
        tmp = data(t,~isnan(data(t,:))); % extract rows without missings
        errors(t) = std(tmp); % set standard error
        if bSED
         errors(t) = errors(t)./sqrt(length(tmp)); % compute standard error of the mean
        end
    end
elseif length(barvalues) ~= length(errors)
    error('barvalues and errors vectors must be of same dimension');
end

% Make scatter dots for each subject (column of input data)
if bScatter
    for t = 1:nbars
        tmp = data(t,~isnan(data(t,:))); % drop rows with missing data
        s = scatter(xpos(t)*ones(1,length(tmp)),tmp,[],'k', 'jitter','on', 'jitterAmount',0.15); % was 0.05
        set(s,'MarkerEdgeColor',[0.4 0.4 0.4],'linewidth',3); % was 1 
    end
end
th=errorbar(xpos,barvalues, errors, 'r', 'linestyle', 'none','linewidth',2.5);
set(th,'Linewidth',2.5,'color',[0 0 0]); % error bars in black
he=[he;th];
removeErrorBarEnds(th)
xlim([0.5 nbars+0.5]);

if~isempty(bw_xtick)
    set(gca,'xtick',1:nbars,'xticklabel',bw_xtick)
end

hold off;

return