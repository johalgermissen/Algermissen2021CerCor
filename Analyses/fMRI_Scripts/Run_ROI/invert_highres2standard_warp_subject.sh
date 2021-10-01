#!/bin/bash

# It can help to centrally invert these matrices for all subjects/ blocks (takes 2 min per matrix), which saves time because they are often used to bring (different) ROIs to native space
# Input: highres2standard_warp in GLM run
# Output: standard2highres_warp in same GLM
# 
# Make executable:
# chmod a+x invert_highres2standard_warp_subject.sh
#
# Execute:
# ./invert_highres2standard_warp_subject.sh # run
# qsub -N "invert_highres2standard_warp_subject" -l walltime=2:00:00,mem=10gb /project/3017042.02/Analyses/fMRI_Scripts/Run_ROI/invert_highres2standard_warp_subject.sh
#
# EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
# J. Algermissen, 2018-2021.

rootdir=/project/3017042.02 # root directory--needs to be adapted to users' folder structure

fsldir=/opt/fsl/6.0.0/bin # FSL's directory--needs to be adapted to users' own directory

GLMID=8H # don't use registration from pre-processing, but new registration computed at 1st-level GLM

# set subject ID, loop
for (( subject=1; subject<=36 ; subject++ )); do

	subjectID=`zeropad $subject 3` # subject ID with 3 digits
	echo "Started subject ${subjectID}"

	# Invert transform from structural to standard space (takes ~2 min.)
	$fsldir/invwarp -w ${rootdir}/Log/fMRI/sub-${subjectID}/GLM${GLMID}/GLM${GLMID}_sub${subjectID}.feat/reg/highres2standard_warp -o ${rootdir}/Log/fMRI/sub-${subjectID}/GLM${GLMID}/GLM${GLMID}_sub${subjectID}.feat/reg/standard2highres_warp -r ${rootdir}/Log/fMRI/sub-${subjectID}/GLM${GLMID}/GLM${GLMID}_sub${subjectID}.feat/reg/highres

	# standard2example_func.mat wrong because only uses FLIRT (not FNIRT)
	# Invert from functional to structural space (highres2example_func.mat) already computed automatically 

	echo "Completed subject ${subjectID}"

done # done with subject
