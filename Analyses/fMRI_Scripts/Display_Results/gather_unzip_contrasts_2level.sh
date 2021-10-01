#!/bin/bash

# Gather thresh_zstat1, cope1, tstat1, and zstat1 files of positive and negative group-level contrast of given GLM,
# unzip all, combine and unzip thresholded masks
# 
# Make executable:
# chmod a+x gather_unzip_contrasts_2level.sh
# 
# Execute:
# ./gather_unzip_contrasts_2level.sh # run
#
# EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
# J. Algermissen, 2018-2021.

GLMID=1  # GLM name plus suffices, like "1_without6"
maxCope=6 # number of copes in this GLM

# Adjust number of copes in for-loop, both for retrieval and for thresholded-map combination!
# if additional postfixes: adjust pos and neg and start of signID for-loop below

# ------------------------------------------------------------------- #
# Define directory where to copy niftis too:
rootdir=/project/3017042.02 # root directory--needs to be adapted to users' own directory structure
targetdir=$rootdir/Log/fMRI/GLM${GLMID}_FEAT_Combined_all

fsldir=/opt/fsl/6.0.0/bin # FSL's directory--needs to be adapted to users' own directory structure

# ------------------------------------------------------------------- #
# Make directory if it doesn't exist yet:
mkdir -p $targetdir
echo "Target directory is $targetdir"

# 1) Copy all relevant files into one common directory:

for ((copeID=1; copeID<=$maxCope; copeID++)); do

	for signID in pos neg; do

		echo "copy cope $copeID $signID"
	
		fmridir=${rootdir}/Log/fMRI/GLM${GLMID}_FEAT_Combined_${signID}.gfeat/cope${copeID}.feat

		# ------------------------------------------------------------------- #
		# Copy & rename thresholded map after contrast:
		cp $fmridir/thresh_zstat1.nii.gz $targetdir/thresh_zstat1_cope${copeID}_${signID}.nii.gz

		# ------------------------------------------------------------------- #
		# Copy & rename cope, tstat, zstat:
		cp $fmridir/stats/cope1.nii.gz $targetdir/cope1_cope${copeID}_${signID}.nii.gz
		cp $fmridir/stats/tstat1.nii.gz $targetdir/tstat1_cope${copeID}_${signID}.nii.gz
		cp $fmridir/stats/zstat1.nii.gz $targetdir/zstat1_cope${copeID}_${signID}.nii.gz

	done
done

# ------------------------------------------------------------------- #
# 2) Unzip all nii.gz files:

echo "Unzip all"

for file in $targetdir/*.nii.gz; do

	gunzip $file

done

# ------------------------------------------------------------------- #
# 3) Combine thresholded maps (and unzip again):

echo "Combine thresholded maps"
cd $targetdir # cd to output directory

for ((copeID=1; copeID<=$maxCope; copeID++)); do

	$fsldir/fslmaths $targetdir/thresh_zstat1_cope${copeID}_pos.nii -add $targetdir/thresh_zstat1_cope${copeID}_neg.nii $targetdir/thresh_zstat1_cope${copeID}_all.nii
	gunzip $targetdir/thresh_zstat1_cope${copeID}_all.nii

done 

echo "Finished :-)"
# END
