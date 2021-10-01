# Readme fMRI

Analysis scripts for analyzing fMRI data and some relevant intermediate files (concatenated pre-processed blocks, regressors per GLM, masks used for ROI analyses) and results (maps from 3 GLMs). 

**Code** for the entire paper will be **maintained** under https://github.com/johalgermissen/Algermissen2021CerCor, with a permanent copy of the code at the time of publication under https://github.com/denoudenlab/Algermissen2021CerCor.

For fMRI **raw data** itself, see the separate collection under https://doi.org/10.34973/pezs-pw62.
**Pre-processed data and results** for this project are available under https://doi.org/10.34973/2t72-bj41.

Mind that in most files, you have to adjust the **root directory** to your own folder structure before running them.

## Scripts ##

We provide analyses scripts under *rootdir/Analyses/fMRI_Scripts*. 

- Folder with **template** .fsf files are:
	- *rootdir/Analyses/fMRI_Scripts/FEAT_Templates* --- contains .fsf template files which allow to create individual .fsf files for each subject.	
 	- *rootdir/Analyses/fMRI_Scripts/FEAT_GLM_Sample_Scripts* --- contains .fsf files for 2nd-level GLMs.
- For **submitting .fsf jobs** on a high performancecomputing cluster, use Run_GLM/fsl_sub_DCCN.sh (see in respective .sh files).
- For **pre-processing**, run scripts in the following order:
	- Run *rootdir/Analyses/fMRI_Scripts/Create_FEAT/feat_create_preAroma.sh* to create preprocessing files for each block for each subject
	- Run *rootdir/Analyses/fMRI_Scripts/Run_Feat/feat_run_preAroma.sh* to run pre-processing files for each block for each subject
	- Run *rootdir/Analyses/fMRI_Scripts/Run_ICA_AROMA/ICA_Aroma_run.sh* to run ICA-AROMA for each block for each subject
	- Run *rootdir/Analyses/fMRI_Scripts/Run_Feat/concat_flirt_blocks.sh* to concatenate blocks for each subject.

- For **first-level (subject level GLMs)**, run scripts in the following order:
	- **Nuisance regressors**: Run scripts in *rootdir/Analyses/fMRI_Scripts/Create_Regressors/Create_Motion_Regressors* to:
		- *create_regressor_CSF.sh*: Create mean CSF signal per subject per block.
		- *create_regressor_OOB.sh*: Create mean OOB signal per subject per block.
		- *Combine_Confound_EVs_Subject.m*: Combine CSF and OOB per subject across blocks, add realignment parameters, add spike regressor when relative displacement above cutoff, add intercept per block.
	- **Task regressors**: use *rootdir/Analyses/fMRI_Scripts/Create_Regressors/Create_Motion_Regressors/EEGfMRIPav_GLMX.m* to create tasks regressor per GLM X for each subject:
		- *GLM1*: Simple GLM with 4 regressors for contrasts valence, action, conflict, reported in main text.
		- *GLM2*: Similar to GLM1, but regressors split per accuracy (correct/incorrect), reported in Supplementary Material S06.
		- *GLM3*: Similar to GLM2, but trial-by-trial EEG regressors (midfrontal alpha power from conflict contrast, midfrontal theta power from action contrast) added, reported in Supplementary Material S17.
			- Requires you to first run *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel/EEGfMRIPav_CueLocked_8_TF_grouplevel_test.m* with contrastType = *'Congruency'* and contrastType = 'Go' to save  clusters as masks, and then *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel/EEGfMRIPav_CueLocked_8_extract_singleTrial.m* on respective masks to create trial-by-trial summary metric as BOLD regressor.
	- **Fitting first-level GLM**: 
	    - Run *rootdir/Analyses/fMRI_Scripts/Create_GLM/feat_create_GLM_Subject_1stlevel.sh* to create .fsf file for 1st-level GLM for each subject.
	    - Run *rootdir/Analyses/fMRI_Scripts/Create_GLM/feat_set_empty_regressors_subject.sh* to insert empty regressor settings into GLM .fsf files.
	    - Run *rootdir/Analyses/fMRI_Scripts/Run_GLM/feat_run_GLM_Subject_1stlevel.sh* to run 1st-level GLM for each subject.
- For **second-level (sample-level) GLMs**, run scripts in the following order:
	- Run *rootdir/Analyses/fMRI_Scripts/Run_GLM/Run_GLM/feat_run_GLM_Sample_2ndlevel.sh* to run 2nd-level GLM for each subject.
	- Run GLM1: anatomical masks for GLMs with small volume correction (SVC) are provided under *rootdir/Log/fMRI/fMRI_Masks/masksHarvardOxford* or can be alternatively created via executing *rootdir/Analyses/fMRI_Scripts/Run_ROI/Extract_ROIs/extract_ROI_HOA_GLM1_SVC.sh*.	
- For creating **dual coded fMRI images** with the slice display toolbox, run scripts in the following order:
	- Run *rootdir/Analyses/fMRI_Scripts/Display_Results/gather_unzip_contrasts_2level.sh* to assemble necessary files from all contrasts from a given GLM.
	    - For GLM1: Run *rootdir/Analyses/fMRI_Scripts/Display_Results/run gather_add_SVC_2ndlevel.sh* to add clusters from GLMs with SVC to relevant thresh_zstat1 file (creates extra new thresh_zstat1 file).
	- Make plots with Slice display toolbox:
		- Adjust directories to your own directory structure/ add paths to relevant toolboxes in *rootdir/Analyses/fMRI_Scripts/Display_Results/slice_display_set_dirs.m*.
		- Run *rootdir/Analyses/fMRI_Scripts/Display_Results/EEGfMRIPav_sliceDisplay_run.m* to create plots for each GLM for relevant slices.
- For creating **bar plots** of mean parameter estimates per regressor in selected ROI, 
	- Masks are provided under *rootdir/Log/fMRI/fMRI_Masks/masksTAfTCueLocked* or can be alternatively created via executing *rootdir/Analyses/fMRI_Scripts/Run_ROI/Extract_ROIs/extract_ROI_HOA_GLM1_TAfT.sh*.
	- Run *rootdir/Analyses/fMRI_Scripts/Run_ROI/invert_highres2standard_warp_subject.sh* to invert the native-space-to-MNI registration per subject.
	- Run *rootdir/Analyses/fMRI_Scripts/Run_ROI/run_ROI_PEs_2ndlevel.sh* to extract summary statistics on PE within mask.
	- Run *rootdir/Analyses/fMRI_Scripts/Run_ROI/Plot_PEs/preprocess_PEs_cueLocked.R* to extract, rearrange, and save mean PE per regressor.
	- Run *rootdir/Analyses/fMRI_Scripts/Run_ROI/Plot_PEs/plot_PEs_ROI_cue_run.m* (which calls *plot_PEs_ROI_cue.m)* to create bar plots of parameter estimates per regressor.
		- Mind to adjust root directory in both *rootdir/Analyses/fMRI_Scripts/plot_PEs_ROI_cue.m* and *rootdir/Analyses/fMRI_Scripts/plot_PEs_ROI_cue_run.m*.
- For creating files to be used as regressors in **fMRI-informed EEG analyses (TAfT)**, run scripts in the following order: 
	- Masks are provided under *rootdir/Log/fMRI/fMRI_Masks/masksTAfTCueLocked* or can be alternatively created via executing *rootdir/Analyses/fMRI_Scripts/Run_ROI/Extract_ROIs/extract_ROI_HOA_GLM1_TAfT.sh*.
	- Run *rootdir/Analyses/fMRI_Scripts/Run_ROI/invert_highres2standard_warp_block.sh* to invert the native-space-to-MNI registration per block.
	- Run *rootdir/Analyses/fMRI_Scripts/Run_ROI/run_ROI_rawdata_TAfT.sh* to extract mean time series from selected ROIs for fMRI-informed EEG analyses (will get stored in AROMA folder of each block of each subject).

## Special notes on fMRI data
- For subject 004, no fieldmap was collected due to time constraints.
- For subjects 015 and 025, (structural-to-standard space) registration failed (likely due to a problem with the T1), so they are excluded from all analyses involving fMRI.
- For subject 025, no DTI data was collected due to time constraints.
- For subject 032, DTI and T1 scan were acquired in a separate session.

End of file.
