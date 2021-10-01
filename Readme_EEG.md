# Readme EEG

Analysis scripts for analyzing EEG data and some relevant data files (rejected channels, rejected trials, independent components) during pre-processing.

**Code** for the entire paper will be **maintained** under https://github.com/johalgermissen/Algermissen2021CerCor, with a permanent copy of the code at the time of publication under https://github.com/denoudenlab/Algermissen2021CerCor.

For EEG **raw data** itself, see the separate collection under https://doi.org/10.34973/pezs-pw62.
**Pre-processed data and results** for this project are available under https://doi.org/10.34973/2t72-bj41.

All folders containing analysis scripts have a *XXX_set_rootdir.m* script in which you can centrally set the **root directory** (*rootdir* or *dirs.root*) to your own folder structure before running them. Note that *.sh* and *.R* scripts cannot call these files, so you have to set the root directory manually in these scripts.

## Scripts ##

We provide analyses scripts under *Analyses/EEG_Scripts/CueLockedAnalyses/*. Run scripts in the following order:
- **Preprocessing** under *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Preprocessing/*:
    - *EEGfMRIPav_0_markInterpolation.m*: Interactive script to visually inspect the data and mark channels for interpolation.
    - *EEGfMRIPav_0_Polhemus_Template.m*: Interactive script to create Polhemus template for subject for which no (valid) Polhemus data is available.
    - *EEGfMRIPav_1_preICA.m*: Script for initial preprocessing up to trial rejection (e.g. epoching, re-referencing, filtering, baseline correction)..
    - *EEGfMRIPav_2_ICA.m: Script to perform ICA.
    - *EEGfMRIPav_3_visualICs.m*: Script for visualization of ICs.
    - *EEGfMRIPav_4_rejICs.m*: Script to for removal of marked ICs.
    - *EEGfMRIPav_5_TR.m*: Script to mark trials for rejection.
    - *EEGfMRIPav_6_finalPP.m*: Script to reject marked trials, interpolate, and apply spatial filter (Laplacian/ current source density).
    - *EEGfMRIPav_CueLocked_7_TF.m*: Script for TF decomposition (both stimulus- and response-locked).
    - *preprocessing_set_rootdir.m*: Set root directory.
- **Group-level time-domain and TF-domain analyses** under *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Grouplevel*:
    - *EEGfMRIPav_CueLocked_8_TF_grouplevel_create.m*: Script for creating object combining relevant TF-domain data across subjects for defined conditions.
    - *EEGfMRIPav_CueLocked_8_TF_grouplevel_create_RT.m*: Script for creating object combining relevant TF-domain data across subjects split by cue valence and by RT bin.
    - *EEGfMRIPav_CueLocked_9_TF_grouplevel_plot.m*: Interactive script to create relevant plots based on TF-domain data
    - *EEGfMRIPav_CueLocked_9_TF_grouplevel_test.m*: Interactive script to perform permutation and other statistical tests on TF-domain data; also for saving cluster above threshold and for exporting data per subject per condition in defined time/frequency/channel ROI.
    - *EEGfMRIPav_CueLocked_9_TF_grouplevel_test_peak.m*: Interactive script to do t-tests on peak height and latency of TF-domain data per subject per condition in defined time/frequency/channel ROI.
    - *EEGfMRIPav_CueLocked_8_time_grouplevel_create.m*: Script for creating object combining relevant time-domain data across subjects for defined conditions.
    - *EEGfMRIPav_CueLocked_8_time_grouplevel_create_BOLD.m*: Script for creating object combining relevant time-domain data across subjects split by BOLD tertile (retain highest/ lowest tertile).
    - *EEGfMRIPav_CueLocked_9_time_grouplevel_plot.m*: Interactive script to create relevant plots based on time-domain data.
    - *EEGfMRIPav_CueLocked_9_time_grouplevel_test.m*: Interactive script to perform permutation and other statistical tests on time-domain data.
    - *EEGfMRIPav_CueLocked_9_RM-ANOVA_t-test.R*: Perform RM-ANOVAs and t-tests on exported data.
    - *EEGfMRIPav_CueLocked_8_extract_singleTrial.m*: Extract trial-by-trial mean in given time/frequency/channel mask exported above, save as regressor for given fMRI GLM.
    - *TF_update_job.m*: Initialize *job* object to load TF-domain data, aggregate for certain contrast in certain time/frequency/channel range.
    - *TF_load_data.m*: Load TF-domain data previously aggregated per subject/ condition.
    - *TF_prepare_contrast_data.m*: Aggregate TF-domain data for certain contrast in certain time/frequency/channel range.
    - *TF_prepare_generic_data.m*: Aggregate TF-domain data per subject across conditions and vice versa irrespective of contrast.
    - *time_update_job.m*: Initialize *job* object to load time-domain data, aggregate for certain contrast in certain time/frequency/channel range.
    - *time_load_data.m*: Load time-domain data previously aggregated per subject/ condition.
    - *time_prepare_contrast_data.m*: Aggregate time-domain data for certain contrast in certain time/frequency/channel range.
    - *time_prepare_generic_data.m*: Aggregate time-domain data per subject across conditions and vice versa irrespective of contrast.
    - *twoLinePlot.m*: Make line plot with separate line (mean plus shade for SEM) for each side of chosen contrast for mean of given frequency band and given set of channels.
    - *multiTopoPlot.m*: Make a series of topoplots for given set of time ranges contrasting two sides of chosen contrast for mean of given time range and frequency range.
    - *multiLinePlot.m*: Make line plot with separate line (mean plus shade for SEM) for each condition for mean of given frequency band and given set of channels.
    - *singleTopoPlot.m*: Make a single topoplot contrasting two sides of chosen contrast for mean of given time range and frequency range.
    - *TFplot.m*: Make time-frequency plot (heat map) contrasting two sides of chosen contrast for mean of given set of channels. 
    - *evaluate_stat.m*: Automatically return time/ frequency/  channel location of clusters above threshold.
    - *permutation_test.m*: take two vectors X and Y and permute their elements N number of times; return number of samples where X $>$ Y as p-value.
    - *grouplevel_set_rootdir.m*: Set root directory.
- **Helper files** under rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/CueLocked_Helpers/: 
    - *barScatter.m*: Create scatter plot.
    - *boundedline.m*: Create line plot with shades for error bars.
    - *helpers_set_rootdir.m*: Set root directory.
    - *load_recode_behav.m*: Load and recode behavior for given subject.
    - *mycsvwrite.m*: Export matrix as .csv file.
    - *redefine_resplock.m*: Re-lock stimulus-locked time-domain data to reaction times (RTs) (for NoGo responses, take mean RT per subject per condition).
    - *regress_baseline.m*: Perform baseline correction of time-domain or TF-domain data by performing linear regression per time/frequency/channel bin over trials and correct for value predicted by linear regression.
    - *set_dirs.m*: Set EEG directories.
    - *set_pars.m*: Set parameters for EEG pre-processing and analyses.
    - *subtract_ERP.m*: Compute mean ERP per condition over trial and subtract it from trial-by-trial time-domain data of respective condition.
    - *withinSE.m*: Compute within-subject standard error per condition using Cousineau-Morey method of correcting for between-subjects differences.
- **fMRI-informed EEG analyses** (TAfT) under *rootdir/Analyses/EEG_Scripts/CueLockedAnalyses/TAfT*: 
    - See specialized *README_TAfT.md* file in that folder.
    - *taft_realignment_per_condition.m*: Compute t-tests to compare mean realignment parameters between certain task conditions.
    - *taft_correlation_regression_BOLD_RT_per_valence_split.m*: Compute linear regression and correlation between trial-by-trial BOLD and RTs separate for Win and Avoid cues for each subject, perform t-tests across subjects, create various plots.

## Special notes on EEG data

- For subject 008, blocks 1 and 2 (trials 1-220) were accidentally collected together without interruption; thus this subject only has 5 different data files.
- For subject 004, Polhemus data seems to be invalid because 71 (instead of 70) channel positions were recorded; not clear how to resolve single channel positions being very off; consider using the mean of all other valid Polhemus files as a template.
- For subject 005, Polhemus data seems to be invalid because 69 (instead of 70) channel positions were recorded; the RESP channel has probably been forgotten, so channel FCz is at position 66 instead of 67.
- For subjects 032, 033, and 034, no Polhemus data were recorded due to unavailability of the recording system after a software update; consider using the mean of all other valid Polhemus files as a template.

End of file.
