# README General

## Content ##
This collection contains the following files:
- **Analyses**: Scripts for all behavioral, EEG, and fMRI analyses required to reproduce the reported results
- **Log**: behavioral, EEG, and fMRI data required to reproduce the reported results

Under *rootdir/Analyses/FiguresCueLocked*, we provide special scripts to reproduce the respective **results/ figures** in the paper.

**Code** for the entire paper will be **maintained** under https://github.com/johalgermissen/Algermissen2021CerCor, with a permanent copy of the code at the time of publication under https://github.com/denoudenlab/Algermissen2021CerCor.

**Raw data** for this project are available under https://doi.org/10.34973/pezs-pw62. In line with requirements of the Ethics Committee and the Radboud University security officer, potentially identifying data (such as imaging data) can only be shared to identifiable researchers. Hence, researchers requesting access to the data have to register and accept a data user agreement; access will then automatically be granted via a click-through procedure (without involvement of authors or data stewards).

**Pre-processed data and results** for this project are available under https://doi.org/10.34973/2t72-bj41.

For more details, see the modality-specific READMEs.

## Root directory ##
Note that for most files, the **root directory** of this project needs to be adjusted. The folders in *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses* have one central file where the root directory can be set centrally (*XX_set_rootdir.m*). In all other files, check whether the root directory must be adjusted at the beginning.

## Steps to do to reproduce the respective results and figures of the manuscript and the Supplementary material
- Figure 1: Behavior: 
    - Mixed-effect models:
        - Run *rootdir/Analyses/Behavior_Scripts/EEGfMRIPav_extract_rawdata.m* to extract and save behavioral data from raw Matlab/ PsychToolbox files. 
        - Run *rootdir/Analyses/Behavior_Scripts/EEGfMRIPav_1_Read_Initial.R* to read data into R, format, save again.
        - Run *rootdir/Analyses/Behavior_Scripts/EEGfMRIPav_3_MixedModels.R* for fitting linear mixed-effects models.
    - Plots: 
	- Run *rootdir/Analyses/Behavior_Scripts/Matlab_Plots/EEGfMRIPav_plot_behavior.m*.
    	- See also *rootdir/Analyses/FiguresCueLocked/Figure01.m*.
- Figure 2. fMRI: 
    - Figure 2A:
        - Run *rootdir/Analyses/fMRI_Scripts/Run_ROI/Plot_PEs/plot_PEs_ROI_cue_run.m* to create model plot.
    - Figure 2B-D:
        - Preprocessing:
            - Run *rootdir/Analyses/fMRI_Scripts/Create_FEAT/feat_create_preAroma.sh* to create pre-processing files.
            - Run *rootdir/Analyses/fMRI_Scripts/Run_FEAT/feat_run_preAroma.sh* to execute pre-processing files.
            - Run *rootdir/Analyses/fMRI_Scripts/Run_ICA_AROMA/ICA_Aroma_run.sh* to execute ICA-AROMA.
            - Run *rootdir/Analyses/fMRI_Scripts/Create_FEAT/concat_flirt_blocks.sh* to concatenate blocks of each subject.
        - Regressors:
            - Run *rootdir/Analyses/fMRI_Scripts/Create_Regressors/Create_Task_Regressors/EEGfMRIPav_GLM1.m* to create task regressors.
            - Run *rootdir/Analyses/fMRI_Scripts/Create_Regressors/Create_Motion_Regressors/create_regressor_CSF.sh* and *rootdir/Analyses/fMRI_Scripts/Create_Regressors/Create_Motion_Regressors/create_regressor_OOB.sh* to create nuisance regressors.
            - Run *rootdir/Analyses/fMRI_Scripts/Create_Regressors/Create_Motion_Regressors/Combine_Confound_EVs_Subject.m* to create single matrix of nuisance regressors on subject level.
        - GLM:
            - Run *rootdir/Analyses/fMRI_Scripts/Create_GLM/feat_create_GLM_Subject_1stlevel.sh* to create first-level GLM files.
            - Run *rootdir/Analyses/fMRI_Scripts/Create_GLM/feat_set_empty_regressors_subject.sh* to set empty regressors for each subject.
            - Run *rootdir/Analyses/fMRI_Scripts/Run_GLM/feat_run_GLM_Subject_1stlevel.sh* to fit first-level GLM 1.
            - Run *rootdir/Analyses/fMRI_Scripts/Run_GLM/feat_run_GLM_Subject_2ndlevel.sh* to fit second-level GLM 1.
        - Plots:
            - Run *rootdir/Analyses/fMRI_Scripts/Display_Results/gather_unzip_contrasts_2level.sh* to collect all relevant files in one single directory.
            - Run *rootdir/Analyses/fMRI_Scripts/Display_Results/gather_add_SVC_2level.sh.sh* to add clusters obtained in GLM with small-volume correction to other clusters obtained with whole-brain GLM.
            - Run *rootdir/Analyses/fMRI_Scripts/Display_Results/EEGfMRIPav_sliceDisplay_run.m* to create plots.
    - Figure 2E-G: 
        - Pre-process fMRI and fit GLM 1 as specified under Figure 2B-D. 
        - Run *rootdir/Analyses/fMRI_Scripts/Run_ROI/Extract_ROIs/extract_ROI_HOA_GLM1_SVC.sh* to extract selected ROI masks.
        - Run *rootdir/Analyses/fMRI_Scripts/Run_ROI/invert_highres2standard_warp_subject.sh* to bring mask to subject-level.
        - Run *rootdir/Analyses/fMRI_Scripts/Run_ROI/run_ROI_PEs_2ndlevel.sh* to extract mean PEs in selected ROIs.
        - Run *rootdir/Analyses/fMRI_Scripts/Run_ROI/Plot_PEs/preprocess_PEs_cueLocked.R* to extract summary statistics and save in matrix.
        - Run *rootdir/Analyses/fMRI_Scripts/Run_ROI/Plot_PEs/plot_PEs_ROI_cue.m* to create plots.
    - Figure 2H: 
        - See steps for Figure 4.
        - Run *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/TAfT/taft_correlation_regression_BOLD_RT_per_valence_split.m* for t-values/ p-values and plots.
    - See also *rootdir/Analyses/FiguresCueLocked/Figure02.m*.
- Figure 3. EEG: 
    - Run all files in *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Preprocessing/* in numerical order.
    - Run *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel/EEGfMRIPav_CueLocked_8_TF_grouplevel_create.m* based on *respSettings = 'Go'*,*actionSettings = 'exeAct'*, *accSettings = 'correct'*, 
    - Obtain permutation test results using *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel/EEGfMRIPav_CueLocked_8_TF_grouplevel_test.m*
    - Create plots using *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel/EEGfMRIPav_CueLocked_8_TF_grouplevel_plot.m*
    - See also *rootdir/Analyses/FiguresCueLocked/Figure03.m*
- Figure 4. TAfT: 
    - After EEG and fMRI preprocessing and running GLM1, run *rootdir/Analyses/fMRI_Scripts/Run_ROI/Extract_ROIs/extract_ROI_HOA_GLM1_TAfT.sh* to extract selected ROI masks.
    - Run *rootdir/Analyses/fMRI_Scripts/Run_ROI/invert_highres2standard_warp_block.sh* to bring masks to block-level.
    - Run *rootdir/Analyses/fMRI_Scripts/Run_ROI/run_ROI_rawdata_TAfT.sh* to extract volume-by-volume BOLD (first eigenvariate) in selected ROIs.
    - Run *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/TAfT/taft_runonce_create_onsets.m*, *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/TAfT/taft_runonce_relignment_mpFellner.m*, and *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/TAfT/taft_runonce_select_trials.m* to create necessary input files for trial onset, realignment parameters, subsets of trials.
    - Run *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/TAfT/taft_preprocess_main.m* to create TF-domain object with betas from multiple linear regression.
    - Run *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/TAfT/taft_postprocess_TF_main.m* to load TF-domain object with betas from multiple linear regression, create plots, perform cluster-based permutation tests. 
    - See also *rootdir/Analyses/FiguresCueLocked/Figure04.m*.
- Figure S01. Results without 6 subjects exclude from TAfT: 
    - See steps under Figures 1, 2, and 3, but exclude subjects specified as under "witoutInvalid" or "without6".
- Figure S02. fMRI Masks: 
    - Created with FSLeyes based on masks specified under *rootdir/Log/fMRI/fMRI_Masks/masksHarvardOxford/* and *rootdir/Log/fMRI/fMRI_Masks/masksTAfTCueLocked/*.
- Figure S03. Regressors and contrast for GLM1: 
    - No plots, only tables, see files for creating regressors for GLMs under *rootdir/Analyses/fMRI_Scripts/Create_Regressors/Create_Task_Regressors/EEGfMRIPav_GLM1.m*.
- Figure S04. List of significant clusters: 
    - Check also under *rootdir/Log/fMRI/GLM1_FEAT_Combined_pos.gfeat/* and *rootdir/Log/fMRI/GLM1_FEAT_Combined_neg.gfeat/*.
- Figure S05. fMRI over time: 
    - See steps for Figure 4.
    - Run *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/TAfT/taft_save_BOLD_HRF_per_trial.m* to fit and save trial-by-trial HRF amplitude for selected ROIs.
    - Run *rootdir/Analyses/fMRI_Scripts/Run_ROI/analyze_BOLD_effect_over_time.R* to obtain results and plots.
    - See also *rootdir/Analyses/FiguresCueLocked/FigureS05.sh*.
- Figure S06. fMRI GLM correct trials only: 
    - See instructions under Figure 2, but perform first-level and second-level GLM based on scripts for GLM2.
   - See also *rootdir/Analyses/FiguresCueLocked/FigureS06.m*.
- Figure S07. EEG ERPs: 
    - For pre-processing, see instructions Figure 3.
    - Run *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel/EEGfMRIPav_CueLocked_8_time_grouplevel_create.m* based on *respSettings = 'Go'*,*actionSettings = 'exeAct'*, *accSettings = 'correct'*. 
    - Obtain permutation test results using *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel/EEGfMRIPav_CueLocked_8_time_grouplevel_test.m*.
    - Create plots using *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel/EEGfMRIPav_CueLocked_8_time_grouplevel_plot.m*.
    - See also *rootdir/Analyses/FiguresCueLocked/FigureS07.m*.
- Figure S08. Alpha power wit ERPs subtracted: 
    - For pre-processing, see instructions Figure 3.
    - Run *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel/EEGfMRIPav_CueLocked_8_TF_grouplevel_create.m* based on *respSettings = 'Go'*,*actionSettings = 'exeAct'*, *accSettings = 'correct'*, *stimERPcor = 'true'*.
    - Obtain permutation test results using *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel/EEGfMRIPav_CueLocked_8_time_grouplevel_test.m*.
    - Create plots using *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel/EEGfMRIPav_CueLocked_8_time_grouplevel_plot.m*.
    - See also *rootdir/Analyses/FiguresCueLocked/FigureS08.m*.
- Figure S09. Alpha power: Correct vs. incorrect
    - For pre-processing, see instructions Figure 3.
    - Run *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel/EEGfMRIPav_CueLocked_8_TF_grouplevel_create.m* based on *respSettings = 'Go'*,*actionSettings = 'exeAct'*, *accSettings = 'bothAcc'*.
    - Export per subject per condition data using *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel/EEGfMRIPav_CueLocked_8_time_grouplevel_test.m* with contrastType = 'Congruency'.
    - Obtain F-test and t-test results using *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel/EEGfMRIPav_CueLocked_9_RM-ANOVA_t-test.R*
    - Plot
    - Check how RM-ANOVA
- Figure S10. Theta power: Correct vs. incorrect: 
    - For pre-processing, see instructions Figure 3.
    - Run *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel/EEGfMRIPav_CueLocked_8_TF_grouplevel_create.m* based on *respSettings = 'Go'*,*actionSettings = 'exeAct'*, *accSettings = 'bothAcc'*.
    - Obtain permutation test results using *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel/EEGfMRIPav_CueLocked_8_time_grouplevel_test.m*.
    - Create plots using *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel/EEGfMRIPav_CueLocked_8_time_grouplevel_plot.m*.
   - See also *rootdir/Analyses/FiguresCueLocked/FigureS10.m*.
- Figure S11. Evidence accumulation: 
   - For pre-processing, see instructions Figure 3.
   - Run *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel/EEGfMRIPav_CueLocked_8_TF_grouplevel_create.m* based on *job.bin.Type = 'RT', job.bin.Num = 3*.
   - Obtain t-test results using *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel/EEGfMRIPav_CueLocked_8_TF_grouplevel_test_peak.m*.
   - Create plots using *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel/EEGfMRIPav_CueLocked_8_TF_grouplevel_plot.m*.
   - See also *rootdir/Analyses/FiguresCueLocked/FigureS11.m*.
- Figure S12. Head motion: 
   - See steps for Figure 4. 
   - Run *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/TAfT/taft_preprocess_main.m* with *mp_Fellner* as additional regressor; also include behavioral regressors *isgo*, *iswin*, *isG2W* by default.
    - Run *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/TAfT/taft_postprocess_TF_main.m* to load TF-domain object with betas from multiple linear regression, create plots, perform cluster-based permutation tests.
    - Run *rootdir/Analyses/EEG_Scripts/TAfT/taft_realignment_per_condition.m* to obtain t-values and p-values comparing realignment parameters between conditions.
    - See also *rootdir/Analyses/FiguresCueLocked/FigureS12*.
- Figure S13. Increase relative to baseline: 
    - For pre-processing, see instructions Figure 3.
    - Run rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel/EEGfMRIPav_CueLocked_8_TF_grouplevel_create.m* based on *respSettings = 'Go'*,*actionSettings = 'exeAct'*, *accSettings = 'correct'*. 
    - Create plots using *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel/EEGfMRIPav_CueLocked_8_TF_grouplevel_plot.m*.
    - - See also *rootdir/Analyses/FiguresCueLocked/FigureS13.m*.
- Figure S14. Handedness: 
   - For pre-processing, see instructions Figure 3.
   - Run *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel/EEGfMRIPav_CueLocked_8_TF_grouplevel_create.m* based on *respSettings = 'Hand'*,*actionSettings = 'exeAct'*, *accSettings = 'correct'*.
    - Obtain permutation test results using *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel/EEGfMRIPav_CueLocked_8_TF_grouplevel_test.m*.
    - Create plots using *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel/EEGfMRIPav_CueLocked_8_TF_grouplevel_plot.m*.
    - See also *rootdir/Analyses/FiguresCueLocked/Figure14.m*.
- Figure S15. TF-domain TAfT other ROIs: 
   - See steps for Figure 4.
   - Save plots for 'GLM1LeftMotorHand' and 'GLM1RightMotorHand'.
   - Also run with 'GLM1LeftPutamenValence' and 'GLM1MedialCaudateValence' instead of 'GLM1StriatumAction'.
   - See also *rootdir/Analyses/FiguresCueLocked/FigureS15.m*.
- Figure S16. Time-domain TAfT: 
    - See steps for Figure 4.
    - Run *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/TAfT/taft_preprocess_main.m* to create time-domain object with betas from multiple linear regression.
    - Run *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/TAfT/taft_postprocess_time_main.m* to load time-domain object with betas from multiple linear regression, create plots, perform cluster-based permutation tests. 
    - See also *rootdir/Analyses/FiguresCueLocked/FigureS16.m*.
- Figure S17. EEG-informed fMRI analyses. 
   - For EEG pre-processing, see instructions Figure 3.
   - Run *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel/EEGfMRIPav_CueLocked_8_TF_grouplevel_create.m* based on *respSettings = 'Go'*,*actionSettings = 'exeAct'*, *accSettings = 'correct'*.
    - Obtain clusters using *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel/EEGfMRIPav_CueLocked_8_TF_grouplevel_test.m* for contrastType = *'Congruency'* and contrastType = 'Go', save clusters as masks.
    - Run *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel/EEGfMRIPav_CueLocked_8_extract_singleTrial.m* on respective masks to create trial-by-trial summary metric as BOLD regressor.
   - See fMRI instructions under Figure 2, but perform first-level and second-level GLM based on scripts for GLM3.
    - See also *rootdir/Analyses/FiguresCueLocked/FigureS17.m*.

END of file.
