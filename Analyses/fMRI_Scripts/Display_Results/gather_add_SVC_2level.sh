#!/bin/bash

# Copy significant cluster from thresh_zstat1_cope.nii.gz file of GLM with SVC (postfix) to thresh_zstat1_cope.nii.gz file of GLM without SVC (no postfix)
# 
# Make executable:
# chmod a+x gather_add_SVC_2level.sh
# 
# Execute:
# ./gather_add_SVC_2level.sh
#
# EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
# J. Algermissen, 2018-2021.

GLMID=5F # GLM 
copeID=1 # number of cope for which to combine thresholded zmaps
signs=pos # sign of contrast, either pos or neg
postfix=striatum # postfix used to identify GLM feat folder with SVC, e.g. striatum or ACC_JBL (without intial underscore)

# Define directory where to copy niftis too:
rootdir=/project/3017042.02 # root directory--needs to be adapted to users' own directory structure
contrast_dir=${rootdir}/Log/fMRI/GLM${GLMID}_FEAT_Combined_all

fsldir=/opt/fsl/6.0.0/bin # FSL's directory--needs to be adapted to users' own directory structure

# Make directory if it doesn't exist yet:
mkdir -p $contrast_dir
echo "Target directory is $contrast_dir"

# ------------------------------------------------------------------- #
# 1) Copy and rename thresholded map:
for signID in ${signs}; do # for both positive and negative (if available):

	echo "copy cope $copeID $signID"
	cope_dir=${rootdir}/Log/fMRI/GLM${GLMID}_FEAT_Combined_${signID}_${postfix}.gfeat/cope${copeID}.feat

	# Rename thresholded map after contrast:
	cp $cope_dir/thresh_zstat1.nii.gz $contrast_dir/thresh_zstat1_cope${copeID}_${signID}_${postfix}.nii.gz

done

# ------------------------------------------------------------------- #
# 2) Unzip all nii.gz files:

echo "Unzip all"

for file in $contrast_dir/*.nii.gz; do

	gunzip $file

done

# ------------------------------------------------------------------- #
# 3) Combine thresholded maps:

echo "Combine thresholded maps"

$fsldir/fslmaths thresh_zstat1_cope${copeID}_${signID}.nii -add $contrast_dir/thresh_zstat1_cope${copeID}_${signID}_${postfix}.nii $contrast_dir/thresh_zstat1_cope${copeID}_all_${postfix}.nii
gunzip $contrast_dir/thresh_zstat1_cope${copeID}_all_${postfix}.nii

echo "Finished :-)"

# END
