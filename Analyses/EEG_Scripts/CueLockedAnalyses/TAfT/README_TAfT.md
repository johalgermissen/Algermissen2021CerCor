README_TAfT.md

# TAfT - Temporal Analysis of fMRI data Toolbox

The Temporal Analysis of fMRI data (TAfT) Toolbox is a toolbox first developed by Tobias U. Hauser (https://github.com/tuhauser/TAfT) as an SPM toolbox and used in [Hauser, Hunt et al., 2015, Journal of Neuroscience](https://www.jneurosci.org/content/35/32/11209). This toolbox has been ported to Fieldtrip (still using SPM for the HRF template and pseudo-inversion of the design matrix) by Johannes Algermissen.

It has been used in [Algermissen et al., 2021, Cerebral Cortex](https://www.biorxiv.org/content/10.1101/2020.09.11.292870v3.abstract) and Algermissen et al., 2021, bioRxiv. See the respective implementation under https://doi.org/10.34973/2t72-bj41.

Code will be **maintained** under https://github.com/johalgermissen/TAfT.
Code for the entire paper will be maintained under https://github.com/johalgermissen/Algermissen2021CerCor, with a permanent copy of the code at the time of publication under https://github.com/denoudenlab/Algermissen2021CerCor.

This toolbox performs the following steps: 
- Use fMRI time series data (volume-by-volume mean of a selected ROI) as inputs;
- Use behavioral regressors (trial-by-trial) as inputs;
- Correct fMRI time series data for nuisance regressors, apply high-pass filtering;
- Upsample fMRI time series data to new TR, epoch into trials;
- Fit HRF to fMRI data of each trial (separately per trial or per block), obtain HRF amplitude (b-weight) per trial;
- Run separate linear regression for each channel/ (frequency)/ time-bin of EEG data across trials using fMRI and behavioral trial-by-trial regressors; separately for each subject;
- Run one-sample t-test across subjects, perform cluster-based sign-flipping permutation test across selected time/(frequency) range for mean signal of selected channels.

The toolbox is structured as follows:

## Adjust ## 

Scripts to be adjusted at the beginning are:
* *taft_set_rootdir.m*: set root directory of where files (code and data) are located.
 
## Runonce ## 

Scripts to be only run once to create files needed for every subsequent creation of a TafT object:
- *taft_runonce_create_onsets.m*: Create vectors containing the onset of each trial (cue-, response-, or outcome-locked), necessary for epoching trials later; for each block of each subject.
- *taft_runonce_relignment_mpFellner.m*: Create vector of volume-by-volume relative displacement (summary metric of realignment parameters as used by [Fellner et al., 2016, NeuroImage](https://www.sciencedirect.com/science/article/pii/S1053811916002366)) to be used as additional regressor; for each block of each subject.
- *taft_runonce_select_trials.m*: Create vector of indices of trials falling into certain behavioral category (Go/NoGo/Win/Avoid trials), to be used to run regression selectively on those trials; for each subject.

Mind setting your own root directory and the paths to SPM and Fieldtrip! 

## Preprocess ##

Create a TAfT object (3D or 4D map of beta weights in Fieldtrip object) for each subject:
- *taft_preprocess_main.m*: Interactive script to specify fMRI and behavioral regressors and selected trials and run wrappers around the other functions.
- *taft_preprocess_initialize_job.m*: Specify all other settings for creating relevant TAfT object; mind setting the root directory and paths to SPM and Fieldtrip! 
- *taft_preprocess_load_EEG.m*: Load EEG data of one subject; select subset of data as necessary:
- *taft_preprocess_load_behavior.m*: Load behavioral data of one subject; recode accuracy, RTs etc.
- *taft_preprocess_load_fMRI.m*: Load fMRI data, call other fMRI functions, concatenate blocks, combine ROIs in one single design matrix.
- *taft_preprocess_wrapper_upsample_fit.m*: For each ROI for each block, call functions to upsample and epoch and fit HRF amplitude.
- *taft_preprocess_filter_upsample_epoch.m*: Regress out nuisance parameters, upsample, epoch.
- *taft_preprocess_fit_HRF_trial.m*: Fit HRF to each trial separately.
- *taft_preprocess_fit_HRF_block.m*: Fit HRF to each trial in one single GLM for entire block.
- *taft_preprocess_combine_EEG_fMRI.m*: Add behavioral regressors to design matrix, perform multiple linear regression across trials for each channel/ (frequency)/ time bin separately.

Mind setting your own root directory and the paths to SPM and Fieldtrip! 

## Postprocess ##  

Load a previously created TAfT object, select and align data, perform cluster-based sign-flipping permutation test across subjects: 
- *taft_postprocess_TF_selectData.m*: Select time/frequency data for selected ROI, align channels across subjects, select valid subjects, perform one-sample t-test across subjects.
- *taft_postprocess_TF_TFplot.m*: For selected ROI for given channels/ time range, perform cluster-based permutation test for mean signal across selected channels, output p-value, create TF-plot of T-values with with "significant clusters" highlighted.
- *taft_postprocess_TF_topoplot.m*: Plot T-values from t-test across subjects as topoplot for given frequency and time range.
- *taft_postprocess_time_selectData.m*: Select time data for selected ROI, align channels across subjects, select valid subjects, perform one-sample t-test across subjects.
- *taft_postprocess_time_ERplot.m*: For selected ROI for given channels/ time range, perform cluster-based permutation test for mean signal across selected channels, output p-value, create line-plot of T-values with with "significant clusters" highlighted.
- *taft_postprocess_time_topoplot.m*: Plot T-values from t-test across subjects as topoplot for given time range.
- *taft_postprocess_subtract.m*: Subtract two beta maps from each other for each subject.
- *taft_postprocess_TF_allSubjects.m*: Create TF plot for beta-map of each subject in subplot grid; handy for detecting outliers.

Mind setting your own root directory and the paths to SPM and Fieldtrip! 

## Other files:
- *taft_set_rootdir.m*: Set root directory for project.
- *taft_save_BOLD_HRF_per_trial.m*: Export trial-by-trial HRF amplitude for given ROI for given subject.
- *clustertf.m*: Perform cluster-based permutation test for 4D (subject x regressors x frequency x time) data. Original code written by Laurence Hunt.
- *ols.m*: Perform OLS regression with t-tests and F-tests on parameters. Original code written by Tim Behrens and Laurence Hunt.
- *taft_findc.m*: Return index value in a matrix closest to a given value. Useful for epoching trials. Original code written by Laurence Hunt.
- *taft_pinv.m*: Compute the pseudo-inverse (pinv) of a design matrix using SPM's pinv function (which checks for rank deficiency).
- *taft_prune.m*: Can be used to prune bridges (voxels only connected by 2 other neighbors to a certain cluster) within a given cluster for greater biological plausibility.
- *z_standardize.m*: Z-standardize vectors in presence of NaN values.

## Maintenance

This toolbox will be maintained under https://github.com/johalgermissen/TAfT.
A permanent fork, containing the files used for [Algermissen et al., 2021, Cerebral Cortex](https://www.biorxiv.org/content/10.1101/2020.09.11.292870v3.abstract), will be under https://github.com/denoudenlab/TAfT.

## Acknowledgments
Thanks to Tobias U. Hauser and Laurence T. Hunt for sharing code!

End of file.
