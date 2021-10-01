function EEGfMRIPav_0_markInterpolation(rootdir)

% Set directories based on root directory, add paths for Fieldtrip
% 
% INPUT:
% rootdir           =  root directory relative to which paths are initialized
% (default '/project/3017042.02')
% 
% OUTPUT:
% saves file 'channels2repair' to disk in dirs.data.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Preprocessing/

%% Setup. 

if ~exist('rootdir','var')
    rootdir = '/project/3017042.02';
    fprintf('rootdir unspecified, assume %s\n',rootdir)
end

fprintf('Add helper files to path\n');
addpath(fullfile(rootdir,'Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Helpers'));

% Set directories and parameters:
dirs    = set_dirs(rootdir);

%% Mark the channels for interpolations.
% Channels of interest: Fz, FCz, Cz

channel2repair              = {
     1 [3, 6, 17, 22, 34, 35, 39, 42]; % F3, C4, Fz, FC2, F2, C1, AF3, FC4 % --> maybe exclude entire subject due to motion?
     2 [6, 10, 21, 26, 61]; % C4, 02, FC1, FC6, FT9
     3 [40]; % AF4
     4 [2, 3, 4, 33, 41, 42, 61]; % FP2, F3, F4, F1, FC3, FC4, FT9
     5 [NaN]; % none
     6 [35, 61]; % C1, FT9
     7 [64]; % CPz
     8 [3, 6, 12, 21, 34, 53, 61]; % F3, C4, F8, FC1, F2, AF7, FT9
     9 [23, 27, 33, 34, 47]; % CP1, CP5, F1, F2, F5
    10 [NaN]; % none
    11 [41]; % FC3
    12 [18, 23, 24, 28, 36, 49, 64]; % Cz, CP1, CP2, CP6, C2, C5, CPz % --> maybe exclude entire subject due to motion???
    13 [10, 15, 23, 24, 28, 36, 58]; % 02, P7, CP1, CP2, CP6, C2, TP8
    14 [4, 21, 35, 61]; % F4, FC1, C1, FT9 
    15 [8]; % P4 
    16 [40]; % AF4
    17 [8, 10, 19, 20, 61]; % P4, 02, Pz, Oz, FT9 % --> MR cleaning didn't always work, consider excluding?? 
    18 [10, 14, 18, 26, 30, 50, 58, 64]; % O2, T8, Cz, FC6, TP10, C6, TP8, CPz
    19 [5, 6, 21, 40, 42, 43, 50, 61]; % C3, C4, FC1, AF4, FC4, CP3, C6, FT9 % --> MR cleaning didn't always, work consider excluding
    20 [4, 22, 31, 34, 37, 42]; % F4, FC2, POz, F2, P1, FC4
    21 [3, 4, 17, 19, 22, 26, 28, 36, 39, 41, 42, 54, 61]; % F3, F4, Fz, Pz, FC2, FC6, CP6, C2, AF3, FC3, FC4, AF8, FT9
    22 [19, 24, 38]; % Pz, CP2, P2
    23 [17, 19]; % Fz, Pz % --> check FCz when recovered, might be bad as well
    24 [4, 18, 22, 39, 61]; % F4, Cz, FC2, AF3, FT9
    25 [4, 36, 40, 41, 42]; % F4, C2, AF4, FC3, FC4 % --> also check FCz if recovered
    26 [29, 30, 36, 37, 64]; % TP9, TP10, C2, P1, CPz
    27 [22, 39, 41, 42]; % FC2, AF3, FC3, FC4 % --> MR cleaning didn work on all trials, consider excluding
    28 [21, 36, 40, 61, 64]; % FC1, C2, AF4, FT9, CPz
    29 [4, 17, 21, 22, 23, 24, 26]; % F4, Fz, FC1, FC2, CP1, CP2, FC6
    30 [5, 9, 18, 20, 21, 35, 43, 49, 49, 54, 61]; % C3, O1, Cz, Oz, FC1, C1, CP3, C5, AF8, FT9
    31 [4, 14, 18, 22, 24, 30, 31, 57, 64]; % F4, T8, Cz, FC2, CP2, TP10, POz, TP7, CPz
    32 [27, 28, 61, 64]; % CP5, CP6, FT9, CPz
    33 [22, 26, 33, 42, 61]; % FC2, FC6, F1, FC4, FT9 % --> MR cleaning didn't work for all trials, consider excluding
    34 [21, 35, 36]; % FC1, C1, C2
    35 [12, 21, 22, 23, 25, 29, 30]; % F8, FC1, FC2, CP1, FC5, TP9, TP10
    36 [23, 58]; % CP1, TP8
    };

%% Save channels2repair:

fprintf('Save channels to exclude\n');

save(fullfile(dirs.data,'channel2repair.mat'),'channel2repair');

end % end of function.