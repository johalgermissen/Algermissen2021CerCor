#!/bin/bash

# Extract mean volume-by-volume signal (time series) within mask for selected ROIs to be used as regressors in fMRI-informed EEG analyses (TAfT).
# Steps involved:
# 1) Run invert_highres2standard_warp_block.sh first to invert FNIRT registration matrix.
# 2) Brings given mask from standard to functional space for each subject for each block.
# 3) Uses fslmeants to get mean volume-by-volume signal within mask (first eigenvariate).
#
# Make executable:
# chmod a+x run_ROI_rawdata_TAfT.sh
# 
# Execute:
# ./run_ROI_rawdata_TAfT.sh # run
#
# Submit as job to the cluster like this:
# qsub -N "JBLIncongruency" -l walltime=4:00:00,mem=5gb /project/3017042.02/Analyses/fMRI_Scripts/Run_ROI/run_ROI_rawdata_TAfT.sh
# Takes roughly 1 min. (5 min. with warp) per subject; 40 (200) minutes in total; use 5 GB RAM

# Directories:
rootdir=/project/3017042.02 # root directory--needs to be adapted to users' own directory structure
fsldir=/opt/fsl/6.0.0/bin # FSL's directory--needs to be adapted to users' own directory structure

# Select mask for which to extract time series: 
region=GLM1vmPFCValenceMan # name for output file
#region=GLM1CingulateAnteriorValence # name for output file
#region=GLM1MedialCaudateValence # name for output file
#region=LeftPutamenValence # name for output file
#region=GLM1CingulateAnteriorAction # name for output file
#region=GLM1StriatumAction # name for output file
#region=GLM1LeftMotorHand # name for output file
#region=GLM1RightMotorHand # name for output file

maskName=$region # name for input file
lockName=CueLocked

# Set subject ID, loop
for (( subject=1 ; subject<=36; subject++ )); do

	subjectID=`zeropad $subject 3` # subject ID with 3 digits
	echo "Started subject ${subjectID}"

	# Set block ID, loop
	for (( block=1 ; block<=6 ; block++ )); do

		blockID=`zeropad $block 1` # block ID with 1 digit
		echo "started subject ${subjectID} block ${blockID}"

		## a) Bring mask from standard to functional space 
		# standard to structural, i.e. previously created invert_highres2standard_warp_block.sh
		# structural to functional, which is given highres2example_func.mat
		# (takes 30 sec.)

		FILE=${rootdir}/Log/fMRI/sub-${subjectID}/FEAT_Block${blockID}.feat/AROMA/${maskName}.txt
		if [ -f "$FILE" ]; then
			echo "File $FILE already exists"
		else
			echo "Apply warp"
			$fsldir/applywarp -i ${rootdir}/Log/fMRI/fMRI_Masks/masksTAfT${lockName}/${maskName} -r ${rootdir}/Log/fMRI/sub-${subjectID}/FEAT_Block${blockID}.feat/reg/example_func -o ${rootdir}/Log/fMRI/sub-${subjectID}/FEAT_Block${blockID}.feat/reg/${region}_func -w ${rootdir}/Log/fMRI/sub-${subjectID}/FEAT_Block${blockID}.feat/reg/standard2highres_warp --postmat=${rootdir}/Log/fMRI/sub-${subjectID}/FEAT_Block${blockID}.feat/reg/highres2example_func.mat
		fi

		## b) Extract mean signal per region (or better first eigenvariate) 
		# (takes 10 sec.)		
		$fsldir/fslmeants -i ${rootdir}/Log/fMRI/sub-${subjectID}/FEAT_Block${blockID}.feat/AROMA/denoised_func_data_nonaggr -o ${rootdir}/Log/fMRI/sub-${subjectID}/FEAT_Block${blockID}.feat/AROMA/${region}.txt -m ${rootdir}/Log/fMRI/sub-${subjectID}/FEAT_Block${blockID}.feat/reg/${region}_func --eig

		echo "Completed subject ${subjectID} block ${blockID}"

	done # done with block

	echo "Completed subject ${subjectID}"

done # end subject loop
