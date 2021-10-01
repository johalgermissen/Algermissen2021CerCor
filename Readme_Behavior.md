# README behavior

Analysis scripts for analyzing behavioral data and some model outputs.

**Code** for the entire paper will be **maintained** under https://github.com/johalgermissen/Algermissen2021CerCor, with a permanent copy of the code at the time of publication under https://github.com/denoudenlab/Algermissen2021CerCor.

For behavioral **raw data** itself, see the separate collection under https://doi.org/10.34973/pezs-pw62.
**Pre-processed data and results** for this project are available under https://doi.org/10.34973/2t72-bj41.

Note to set the **root directory** in respective script to your own folder structure before running them.

## Scripts
Analyses scripts are provided under *rootdir/Analyses/Behavior_Scripts*.
- Steps to recreate behavioral results:
	- *EEGfMRIPav_extract_rawdata.m*: Run to convert .mat files under *rootdir/Log/Behavior/Data_beh_mat* into .csv files under rootdir/Log/Behavior/Data_beh_csv
	- *EEGfMRIPav_1_Read_Initial.R*: Read all csv files, recode variables, concatenate subjects and store them as EEGfMRIPav_all.csv under Log/Behavior/Behavior_concatenated/EEGfMRIPav_all.csv
	- *EEGfMRIPav_2_Preprocess_Automated.R*: When opening R, refresh variable coding (e.g. ordered factors) and create more convenience variables.
	- *EEGfMRIPav_3_LMEMs.R*: Fit mixed-effects models.

End of file.
