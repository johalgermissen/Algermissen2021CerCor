function [par] = set_par()

% Set parameters for EEG pre-processing and analyses.
%
% INPUT:
% None
% 
% OUTPUT:
% par      - structure with settings for pre-processing
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Should work in Matlab 2018b.

%% Set task and session name, session number, trigger number, epoching times, condition labels and titles:

% ASCII response codes: 65 - 72, 97 - 104
% https://intranet.donders.ru.nl/index.php?id=lab-response-fiberoptic&no_cache=1&no_cache=1&sword_list%5B%5D=fiber

%% Set global parameter settings:

par.nSub    = 36;
par.nBlock  = 6;

% Sampling rate in Hz:
fprintf('Set sampling rate\n');
par.sRate                   = 1000; % CTF sampling rate

%% Channel information:

par.chan.nChan              = 67;
par.chan.recordChan         = 1:66;
par.chan.nRecordChan        = length(par.chan.recordChan); % number all recorded channels
par.chan.EEG                = [1:31 33:64 67]; % 67 to-be recovered by re-referencing later
par.chan.nEEGchan           = length(par.chan.EEG); % number EEG channels
par.chan.REF                = {'Ref'}; %  FCz is reference, called 'Ref' on Easycap sheet, but not reported in data --> create after re-referencing
par.chan.ECG                = 32;
par.chan.HR                 = 65;
par.chan.RESP               = 66;
par.chan.layoutfile         = 'easycapM11.mat'; % 2D layout
par.chan.positionfile       = 'easycap-M1.txt'; % 3D coordinates

% % Check layout file:
% cfg                         = [];
% cfg.layout                  = par.chan.layoutfile;
% par.layout                  = ft_prepare_layout(cfg);
% figure; ft_plot_lay(par.layout); title(cfg.layout);

% 2D image:
% the ref position is ignored and should be added; otherwise electrode positions will be inaccurate
% channelpositions            = ft_read_sens(par.chan.positionfile);
% figure;  ft_plot_sens(channelpositions,'label','label'); title('channel positions');

%% Epoching:

fprintf('Trial length for epoching\n');

for iCue = 1:8
    par.epoch.eventCode{iCue}     = sprintf('S%d',110+iCue); % s in ft
end
par.epoch.epochtime               = [-1.75 2.8]; % 1.5s before baseline & after fb offset.
par.epoch.time4TR                 = [-0.25 1.3]; % critial time window for artifact rejection.   

% Trial time course:
%   Cue + response window: 1300 ms
%   Fixation cross: 1400 - 1600 ms
%   Feedback:               750 ms
%   ITI:            1250 - 2000 ms
%   --> Total trial duration: 4700 - 5600 ms

% -------------------------------------------------- %
% Re-referencing:
par.epoch.implicitRef           = 'FCz';

% -------------------------------------------------- %
% Filtering:

% Band pass filter settings:
fprintf('Set filter settings\n');
% High pass filter in Hz: filter < 1 Hz
par.epoch.hpf                     = 0.5;
% Low pass filter in Hz: filter > 35 Hz
par.epoch.lpf                     = 15;

% -------------------------------------------------- %
% Linear baseline correction before ICA: 

fprintf('Set time window for baseline correction\n');
par.epoch.base4correction         = [-.2 0];

%% Time frequency decomposition:

fprintf('Set TF decomposition parameters\n');

% General settings:
par.TF.TFtype           = 'hanning'; % morlet, hanning

% --------------------------------------- % 
% Frequencies:
par.TF.nFreqs           = 15;
par.TF.minFreq          = 1; % lowest frequency
par.TF.stepFreq         = 1; % distance between frequencies
par.TF.maxFreq          = par.TF.minFreq + par.TF.nFreqs * par.TF.stepFreq; % highest frequency

% Frequencies to decompose:
par.TF.freq4tf           = par.TF.minFreq:par.TF.stepFreq:par.TF.maxFreq; % linearly spaced
% par.TF.freq4tf         = logspace(log10(par.TF.minFreq),log10(par.TF.maxFreq),par.TF.nFreq); % logarithmically spaced

% --------------------------------------- %
% b) Time points to decompose on:

% Set adaptively based on lockSettings in actual TF decomposition file
% if strcmp(lockSettings,'stimlocked')
%     parTF.toi4tf          = -1:.025:2; % from 1 sec before till 2 sec after cue presentation; windows of 0.025 sec
% elseif strcmp(lockSettings,'resplocked')
%     parTF.toi4tf          = -1.5:.025:1;
% else
%     error('Invalid lock settings')
% end

% --------------------------------------- %
% c) Baseline:
par.TF.baselinetime      = [-0.25 -.05]; % Suggestion Jenn
% par.TF.baselinetime    = [-1 -.50]; % Alternative: further away from zero
% par.TF.baselinetime    = 0; % instantaneously at zero

% ---------------------------------------
% d) Morlet wavelet specific settings:
par.TF.width4tf          = 4; % use 4 cycles

% ---------------------------------------
% e) Hanning taper specific settings:
% par.TF.ftimwin           = 1/par.TF.freqRes; % compute length of time window based on frequency resolution
par.TF.ftimwin           = 0.4; % or set manually
par.TF.pad               = 8; % pad up to 4 sec. to get well-behaved frequency bins
% parTf.pad               = ceil(max(cellfun(@numel,data.time)/data.fsample)); % determine adaptively

%% Add default color map:
% needs directory as input

% Add folder with color maps:

set(0, 'DefaultFigureColormap', jet(64)); par.cmap = 'jet';  
% set(0, 'DefaultFigureColormap', parula(64)); par.cmap = 'parula';  

% Add custom color maps:
% addpath(fullfile('/project/3017042.02/colorMaps'));
% set(0, 'DefaultFigureColormap', redblue(64)); par.cmap = 'redblue'; 
% set(0, 'DefaultFigureColormap', turbo); par.cmap = 'turbo'; 

fprintf('Add color map %s \n',par.cmap)
% END