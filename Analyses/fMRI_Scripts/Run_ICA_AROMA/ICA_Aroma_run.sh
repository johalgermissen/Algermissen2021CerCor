#!/bin/bash

# Performs ICA AROMA on pre-processed data looping over blocks of subjects (submits single jobs to the computing cluster).
#
# Make executable:
# chmod a+x ICA_Aroma_run.sh
# 
# Submit job like this:
# qsub -N "ICA_Aroma_run" -l walltime=0:00:10,mem=1gb ICA_Aroma_run.sh # run
#
# EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
# J. Algermissen, 2018-2021.

rootdir=/project/3017042.02 # needs to be adapted to users' own directory structure
fsldir=/opt/fsl/6.0.0/bin # FSL's directory--needs to be adapted to users' own directory structure
submitdir=${rootdir}/Analyses/fMRI_Scripts/Run_GLM # where fsl_sub_DCCN.sh is
AROMAdir=/home/action/johalg  # needs to be adapted to users' own directory structure (download ICA-AROMA from https://github.com/maartenmennes/ICA-AROMA)
pythondir=/opt/python/2.7.8/bin/python # needs to be adapted to users' own directory structure

for (( subject=1 ; subject<=36 ; subject++ )); do # loop throught these subjects

	subjectID=`zeropad $subject 3`; # 3 digit subject ID

	for (( block=1 ; block<=6 ; block++ )); do # loop throught these blocks

		blockID=`zeropad $block 1`; # 1 digit block ID

		# 1) Make new target directory:
		mkdir -p $rootdir/Log/fMRI/sub-${subjectID}/FEAT_Block${blockID}.feat/AROMA

		# 2) Create extra mask (as recommended in ICA-AROMA Manual):
		$fsldir/bet $rootdir/Log/fMRI/sub-${subjectID}/FEAT_Block${blockID}.feat/example_func $rootdir/Log/fMRI/sub-${subjectID}/FEAT_Block${blockID}.feat/mask_aroma.nii.gz -f 0.3 -n -m -R
	
		# 3) Run ICA-AROMA: just specify FEAT directory as input, directory created above as output, extra mask, nonaggressive denoising 
		$submitdir/fsl_sub_DCCN.sh -N "Aroma_Sub${subjectID}_Bl${blockID}" -T 300 -u 15gb $pythondir/python2.7 $AROMAdir/ICA-AROMA-master/ICA_AROMA.py -feat $rootdir/Log/fMRI/sub-${subjectID}/FEAT_Block${blockID}.feat -out $rootdir/Log/fMRI/sub-${subjectID}/FEAT_Block${blockID}.feat/AROMA -m $rootdir/Log/fMRI/sub-${subjectID}/FEAT_Block${blockID}.feat/mask_aroma.nii.gz -den nonaggr;

		# Alternatively submit as job to the cluster:
		# ~/fsl_sub_DCCN.sh -N "Aroma_Sub${subjectID}_Bl${blockID}" -T 200 -u 15gb $pythondir/python2.7 $AROMAdir/ICA-AROMA-master/ICA_AROMA.py -feat $rootdir/Log/fMRI/sub-009/FEAT_Block3.feat -out $rootdir/Log/fMRI/sub-009/FEAT_Block3.feat/AROMA -m $rootdir/Log/fMRI/sub-009/FEAT_Block3.feat/mask_aroma.nii.gz -den nonaggr;

	done # end block loop

done # end subject loop

echo "Finished :-)"
# END
