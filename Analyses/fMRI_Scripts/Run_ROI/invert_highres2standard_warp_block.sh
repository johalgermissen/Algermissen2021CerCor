#!/bin/bash

# It can help to centrally invert these matrices for all subjects/ blocks (takes 2 min per matrix), which saves time because they are often used to bring (different) ROIs to native space
# 
# Make executable:
# chmod a+x invert_highres2standard_warp_block.sh # make executable
#
# Execute:
# ./invert_highres2standard_warp_block.sh # run
# 
# Submit as job to cluster:
# qsub -N "invert_highres2standard_warp_block" -l walltime=10:00:00,mem=10gb ${rootdir}/Analyses/fMRI_Scripts/Run_ROI/invert_highres2standard_warp_block.sh
#
# EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
# J. Algermissen, 2018-2021.

rootdir=/project/3017042.02 # root directory--needs to be adapted to users' folder structure

fsldir=/opt/fsl/6.0.0/bin # FSL's directory--needs to be adapted to users' own directory

# set subject ID, loop
for (( subject=1 ; subject<=36 ; subject++ )); do

	subjectID=`zeropad $subject 3` # subject ID with 3 digits
	echo "Started subject ${subjectID}"

	# Set block ID, loop
	for (( block=1 ; block<=6 ; block++ )); do

		blockID=`zeropad $block 1` # block ID with 1 digit
		echo "Started subject ${subjectID} block ${blockID}"

		# Invert transform from structural to standard space (takes ~2 min.)
		echo "Invert warp"
		$fsldir/invwarp -w ${rootdir}/Log/fMRI/sub-${subjectID}/FEAT_Block${blockID}.feat/reg/highres2standard_warp -o ${rootdir}/Log/fMRI/sub-${subjectID}/FEAT_Block${blockID}.feat/reg/standard2highres_warp -r ${rootdir}/Log/fMRI/sub-${subjectID}/FEAT_Block${blockID}.feat/reg/highres

		# standard2example_func.mat wrong because only uses FLIRT (not FNIRT)
		# Invert from functional to structural space (highres2example_func.mat) already computed automatically 

	echo "Completed subject ${subjectID} block ${blockID}"

	done # end block loop

echo "completed subject ${subjectID}"

done # end subject loop
