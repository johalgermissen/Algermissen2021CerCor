function taft_runonce_relignment_mpFellner()

% taft_runonce_relignment_mpFellner()
% 
% Loop over subjects and blocks, load realignment parameters from FSL
% computed during motion correction (prefiltered_func_data_mcf.par),
% standardize, compute difference to previous trial, square, 
% compute summary across all six parameters (mean relative displacement),
% save under subject-/ block-specific 'AROMA' directory in fMRI directory.
% Mind setting the root directory in taft_set_rootdir().
% 
% Metric suggested by:
% Fellner, M. C., Volberg, G., Mullinger, K. J., Goldhacker, M., Wimber, M.,
% Greenlee, M. W., & Hanslmayr, S. (2016). Spurious correlations in 
% simultaneous EEG-fMRI driven by in-scanner movement. 
% Neuroimage, 133, 354-366.
%
% EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
% J. Algermissen, 2018-2021.
% Adapted from Tobias Hauser (https://github.com/tuhauser/TAfT).
% Should work in Matlab 2018b.

% we are here:
% cd /project/3017042.02/Analyses/EEG_Scripts/CueLockedAnalyses/TAfT/

%% Set directories:

dirs.root   = taft_set_rootdir(); % /project/3017042.02
dirs.log    = fullfile(dirs.root,'Log');
dirs.behav  = fullfile(dirs.log,'Behavior/Data_beh_mat');    
dirs.fMRI   = fullfile(dirs.log,'fMRI');    

%% Settings:

nSub        = 36;
nBlock      = 6;

%% Loop through subjects and blocks:

for iSub = 1:nSub % iSub = 1;

    for iBlock = 1:nBlock % iBlock = 1;
    
        fprintf('Subject %03d Block %d: Load data\n',iSub,iBlock);
        
        % 1) Load realignment parameters:
        motionname          = fullfile(dirs.fMRI,sprintf('sub-%03d/FEAT_Block%d.feat/mc/prefiltered_func_data_mcf.par',iSub,iBlock));        
        realignPars         = textread(motionname); % read in file with relative displacement parameters
        realignPars(:,7)    = []; % delete last column (always only zeros)

        % 2) Initialize output object:
        mp                  = nan(size(realignPars));
        
        % 3) Standardize realignment parameters:
        mean_mp         = mean(realignPars,1); % mean of entire time series
        std_mp          = std(realignPars,1); % sd of entire time series

            
        % 4) Loop over realignment parameters, loop over trials,
        % standardize, derivative, squared:
        for iParam = 1:size(realignPars,2)
            for iTrial = 1:(size(realignPars,1)-1)
                
                % Subtract mean of entire time series, divide by std of
                % entire time series.
                % Next trial minus this trial, squared
                mp(iTrial,iParam) = ((realignPars(iTrial+1,iParam)-mean_mp)/std_mp - (realignPars(iTrial,iParam)-mean_mp)/std_mp).^2;

            end % end iTrial.            
        end % end iParam.
        
        % 5) Set last trial to zero (no further volume to determine motion):
        mp(end,:)           = mp(end-1,:); % fill in second-to-last trial

        % 6) Sum across all 6 parameters:
        sum_mp              = sum(mp,2); % mean over parameters
        
        % 7) Save:
        fprintf('Subject %03d Block %d: Save data\n',iSub,iBlock);
        dirs.target = fullfile(dirs.fMRI,sprintf('sub-%03d/FEAT_Block%d.feat/AROMA',iSub,iBlock));        
        save(fullfile(dirs.target,'mp_Fellner.txt'),'sum_mp','-ascii')

    end % end iBlock
end % end iSub

end % end of function.
