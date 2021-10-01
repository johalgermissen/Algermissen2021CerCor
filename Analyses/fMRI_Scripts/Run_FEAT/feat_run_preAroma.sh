#!/bin/bash

# Automatically run .fsf files for preprocessing based on .fsf created per block per subject (see folder "Create_FEAT")
#
# Make executable:
# chmod a+x feat_run_preAroma.sh
# 
# Execute:
# ./feat_run_preAroma.sh
#
# EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
# J. Algermissen, 2018-2021.

rootdir=/project/3017042.02 # root directory--needs to be adapted to users' folder structure
submitdir=${rootdir}/Analyses/fMRI_Scripts/Run_GLM # where fsl_sub_DCCN.sh is

for (( subject=1 ; subject<=36 ; subject++ )); do # loop throught these subjects

	subjectID=`zeropad $subject 3`; # 3 digit subject ID

	for (( block=1 ; block<=6 ; block++ )); do # loop throught these blocks

		blockID=`zeropad $block 1`; # 1 digit block ID

		# use fsl_sub_DCCN.sh (in my home directory), run on -T, minutes, storage space, feat command, file to run (with complete path)
		$submitdir/fsl_sub_DCCN.sh -N "FEAT_Sub${subjectID}_Bl${blockID}" -T 90 -u 7gb feat /project/3017042.02/Analyses/fMRI_Scripts/FEAT_preAroma_Scripts/feat_Sub${subjectID}_Block${blockID}_preAroma.fsf;

	done # end block loop

done # end subject loop

# END
