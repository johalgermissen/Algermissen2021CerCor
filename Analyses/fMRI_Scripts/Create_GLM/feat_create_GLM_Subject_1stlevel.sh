#!/bin/bash

# Automatically creates .fsf files for running 1st-level subject-wise GLM based on GLM-specific template.
# 
# Make executable:
# chmod a+x feat_create_GLM_Subject_1stlevel.sh # make executable
# 
# Execute:
# ./feat_create_GLM_Subject_1stlevel.sh
#
# EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
# J. Algermissen, 2018-2021.

rootdir=/project/3017042.02 # root directory--needs to be adapted to users' folder structure

GLMID=1 # specify GLM number

# Create directory if it doesn't exist yet 
mkdir -p ${rootdir}/Analyses/fMRI_Scripts/FEAT_GLM${GLMID}_Subject_Scripts

# Loop over subjects
for subject in {1..36}; do

	subjectID=`zeropad $subject 3` # subject ID with 3 digits
	echo "Start subject ${subjectID}"

	# Get number of volumes for this file: use fslinfo, first five rows, last row of these, print number of volumes to variable npts
	npts=`fslinfo ${rootdir}/Log/fMRI/sub-${subjectID}/postAROMA.nii.gz | head -n 5 | tail -n 1 | awk '{print $2}'`
	
	# Specify template, number of time points, number of subject, output variable
	cat ${rootdir}/Analyses/fMRI_Scripts/FEAT_Templates/feat_template_GLM${GLMID}_subject.fsf | sed s/TIMEPOINTS/$npts/g | sed s/SUBJECT/$subjectID/g  > ${rootdir}/Analyses/fMRI_Scripts/FEAT_GLM${GLMID}_Subject_Scripts/feat_GLM${GLMID}_Subject_Sub${subjectID}.fsf

done # end subject loop
