#!/bin/bash

# Extract volume-by-volume cerebrospinal fluid (CSF) signal as suggested by from Maarten Mennes.
# Saved as CSF_noise.txt in AROMA folder of specific block of specific subject.
# 
# Make executable:
# chmod a+x create_regressor_CSF.sh
# 
# Submit as job to the cluster:
# qsub -N "create_regressor_CSF_01_36" -l walltime=03:00:00,mem=10gb create_regressor_CSF.sh
# takes ~ 4 min. for segmentation, then 20 seconds per block, so 6 min in total; 1.1 GB --> around 3.5h for all subjects --# Tip: do in runs of 10 subjects each to speed up
#
# EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
# J. Algermissen, 2018-2021.

rootdir=/project/3017042.02 # root directory--needs to be adapted to users' own directory structure
fsldir=/opt/fsl/6.0.0/bin # FSL's directory--needs to be adapted to users' own directory structure

#GLMID=1
for (( subject=1; subject<=36; subject++)); do 
	subjectID=`zeropad $subject 3`;

	# 1) Create temporary segmentation folder to-be-deleted later:
	mkdir -p ${rootdir}/Log/fMRI/sub-${subjectID}/NIFTI/Segmentation 

	# 2) Copy T1 into segmentation folder:
	cp -a ${rootdir}/Log/fMRI/sub-${subjectID}/NIFTI/sub-${subjectID}_T1_co_brain.nii.gz ${rootdir}/Log/fMRI/sub-${subjectID}/NIFTI/Segmentation/sub-${subjectID}_T1_co_brain.nii.gz

	# 3) Do probabilistic (!) Segmentation on T1: 
	# creates as output sub-001_T1_co_brain_prob_0.nii.gz sub-001_T1_co_brain_prob_1.nii.gz sub-001_T1_co_brain_prob_2.nii.gz

	echo "Start segmentation Subject ${subjectID}"
	$fsldir/fast --channels=1 --type=1 --class=3 -p ${rootdir}/Log/fMRI/sub-${subjectID}/NIFTI/Segmentation/sub-${subjectID}_T1_co_brain.nii.gz

	echo "Finish segmentation"

	# Loop over blocks, create CSF mask in native space, extract mean signal:
	for (( block=1 ; block<=6 ; block++ )); do 
		blockID=`zeropad $block 1`; 

		echo "Start Subject ${subjectID} Block ${blockID}"

		# Create temporary folder Nuisance_CSF to-be-deleted at the end:
		mkdir -p ${rootdir}/Log/fMRI/sub-${subjectID}/FEAT_Block${blockID}.feat/Nuisance_CSF 

		# 1) Threshold CSF segmentation (binarize into 1 for CSF/ 0 for non-CSF)
		# infile: sub-001_T1_co_brain_prob_0.nii.gz (prob_0 is CSF)
		# outfile: th_csf_mask.nii.gz
		echo "Start Subject ${subjectID} Block ${blockID}: Start segmentation of T1"
		$fsldir/fslmaths ${rootdir}/Log/fMRI/sub-${subjectID}/NIFTI/Segmentation/sub-${subjectID}_T1_co_brain_prob_0.nii.gz -thr 0.95 ${rootdir}/Log/fMRI/sub-${subjectID}/FEAT_Block${blockID}.feat/Nuisance_CSF/th_csf_mask.nii.gz

		# 2) Bring thresholded CSF segmentation from structural to functional space using highres2example_func
		# infile: th_csf_mask.nii.gz
		# outfile: th_csf_mask_func.nii.gz
		echo "Start Subject ${subjectID} Block ${blockID}: Bring thresholded CSF segmentation from structural to functional space"
		$fsldir/flirt -in ${rootdir}/Log/fMRI/sub-${subjectID}/FEAT_Block${blockID}.feat/Nuisance_CSF/th_csf_mask.nii.gz -ref ${rootdir}/Log/fMRI/sub-${subjectID}/FEAT_Block${blockID}.feat/example_func -applyxfm -init ${rootdir}/Log/fMRI/sub-${subjectID}/FEAT_Block${blockID}.feat/reg/highres2example_func.mat -out ${rootdir}/Log/fMRI/sub-${subjectID}/FEAT_Block${blockID}.feat/Nuisance_CSF/th_csf_mask_func.nii.gz

		# 3) Bring csfprior from standard to functional space, using standard2highres_warp (created myself) and highres2example_func.mat (automatically created):
		# infile: avg152T1_csf (from MMs personal folder)
		# outfile: native_csf_prior.nii.gz (into Nuisance_CSF folder)
		# uses mean_func and highres2example_func.mat from pre-processing output
		echo "Start Subject ${subjectID} Block ${blockID}: Bring csfprior to functional space"
		$fsldir/applywarp --in=${fsldir}/data/standard/tissuepriors/avg152T1_csf --out=${rootdir}/Log/fMRI/sub-${subjectID}/FEAT_Block${blockID}.feat/Nuisance_CSF/native_csf_prior.nii.gz --ref=${rootdir}/Log/fMRI/sub-${subjectID}/FEAT_Block${blockID}.feat/example_func --warp=${rootdir}/Log/fMRI/sub-${subjectID}/FEAT_Block${blockID}.feat/reg/standard2highres_warp --postmat=${rootdir}/Log/fMRI/sub-${subjectID}/FEAT_Block${blockID}.feat/reg/highres2example_func.mat

		# 4) Multiply functional space prior with functional space segmentation, threshold and binarize to create final mask
		# infile: th_csf_mask.nii.gz
		# outfile: post_th_csf_mask_func.nii.gz
		echo "Start Subject ${subjectID} Block ${blockID}: Combine prior and mask from segmentation"
		$fsldir/fslmaths ${rootdir}/Log/fMRI/sub-${subjectID}/FEAT_Block${blockID}.feat/Nuisance_CSF/th_csf_mask_func.nii.gz -mul ${rootdir}/Log/fMRI/sub-${subjectID}/FEAT_Block${blockID}.feat/Nuisance_CSF/native_csf_prior.nii.gz -thrP 95 -bin ${rootdir}/Log/fMRI/sub-${subjectID}/FEAT_Block${blockID}.feat/Nuisance_CSF/post_th_csf_mask_func.nii.gz

		# 5) Extract time series of first eigenvariate values (mean of 0):
		# Input filtered functional data and mask, get CSF_noise.txt as output: one column with mean CSF activity per volume
		echo "Start Subject ${subjectID} Block ${blockID}: Extract eigenvariate"
		$fsldir/fslmeants -i ${rootdir}/Log/fMRI/sub-${subjectID}/FEAT_Block${blockID}.feat/AROMA/denoised_func_data_nonaggr -o ${rootdir}/Log/fMRI/sub-${subjectID}/FEAT_Block${blockID}.feat/Nuisance_CSF/CSF_noise.txt -m ${rootdir}/Log/fMRI/sub-${subjectID}/FEAT_Block${blockID}.feat/Nuisance_CSF/post_th_csf_mask_func.nii.gz

		# 5) Copy output into AROMA folder (usable across GLMs)
		echo "Start Subject ${subjectID} Block ${blockID}: Copy to AROMA folder"
		cp  ${rootdir}/Log/fMRI/sub-${subjectID}/FEAT_Block${blockID}.feat/Nuisance_CSF/CSF_noise.txt ${rootdir}/Log/fMRI/sub-${subjectID}/FEAT_Block${blockID}.feat/AROMA/CSF_noise.txt

		# 6) Optional: delete Nuisance_CSF folder again to save space
		echo "Start Subject ${subjectID} Block ${blockID}: Delete nuisance CSF folder"
		rm -rf ${rootdir}/Log/fMRI/sub-${subjectID}/FEAT_Block${blockID}.feat/Nuisance_CSF
		echo "Finished Subject ${subjectID} Block ${blockID}"

	done; # end of block loop

	# 7) Optional: delete segmentation folder again to save space
	echo "Start Subject ${subjectID} Block ${blockID}: Delete segmentation folder"
	rm -rf ${rootdir}/Log/fMRI/sub-${subjectID}/NIFTI/Segmentation

done; # end of subject loop
echo "Finished :-)"

