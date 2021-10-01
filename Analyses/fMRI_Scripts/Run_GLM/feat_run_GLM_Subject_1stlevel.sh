#!/bin/bash

# Runs first-level GLM per subject based on .fsf file (see folder "Create_GLM")
#
# Make executable:
# chmod a+x feat_run_GLM_Subject_1stlevel.sh
#
# Loop over subjects, submit job per subject:
# ./feat_run_GLM_Subject_1stlevel.sh
#
# EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
# J. Algermissen, 2018-2021.

rootdir=/project/3017042.02 # root directory--needs to be adapted to users' folder structure
submitdir=${rootdir}/Analyses/fMRI_Scripts/Run_GLM # where fsl_sub_DCCN.sh is

GLMID=1
# no 15, no 25
#for variants of for-loops in bash, see: https://linuxize.com/post/bash-for-loop/
# For all subjects:
for subject in {1..36};
do # loop throught these subjects
	subjectID=`zeropad $subject 3`; # 3 digit subject ID
#	subjectID=002
	# Copy registration to feat directories if necessary
		# 70 hours: -T 4200 # maximum that works 
		# 48 hours: -T 2400 # also still works
		# 24 hours: -T 1200 # also still works
		# 22 hours fine for most subjects
		# 44 GB enough...? 
		# use 2400 min 50 gb as default

		 $submitdir/fsl_sub_DCCN.sh -N "GLM${GLMID}_Sub${subjectID}" -q "verylong" -T 2400 -u 50gb feat ${rootdir}/Analyses/fMRI_Scripts/FEAT_GLM${GLMID}_Subject_Scripts/feat_GLM${GLMID}_Subject_Sub${subjectID}.fsf;

done # end subject loop
