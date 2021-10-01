#!/bin/bash

# Sets regressors in .fsf scripts to empty if specified as empty in emptyregressors.txt file in specific 
# timings_regressors folder of respective GLM of respective subject.
# 
# Make executable:
# chmod a+x feat_set_empty_regressors_subject.sh
# 
# Execute:
# ./feat_set_empty_regressors_subject.sh
#
# EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
# J. Algermissen, 2018-2021.

rootdir=/project/3017042.02 # root directory--needs to be adapted to users' folder structure

GLMID=1
nRegressors=11 # 5
# nRegressors needs to be number of regressors created in Matlab file and contained in emptyregressors.txt 

# Create directory if it doesn't exist yet:
mkdir -p ${rootdir}/Analyses/fMRI_Scripts/FEAT_GLM${GLMID}_Subject_Scripts

# Set subject ID, loop:
for (( subject=1 ; subject<=36 ; subject++ )); do

	subjectID=`zeropad $subject 3` # subject ID with 3 digits
	echo "start subject ${subjectID}"
	
	# Load empty regressor file:
	file="${rootdir}/Log/fMRI/sub-${subjectID}/GLM${GLMID}/timings_regressors/emptyregressors.txt" # set file name
	emptyregressors=$(cat "$file") # store in variable

	# Check for each regressor whether 1:
	for (( regressorID=1 ; regressorID<=${nRegressors} ; regressorID++ )); do
		position=$(($regressorID  * 2 - 2)) # given the commas, we need positions 0, 2, 4, ...
		status=${emptyregressors:$position:1} # from this position, only 1 character
		# echo "Regressor No. ${regressorID} has value $status"
		if (($status == 1)); then

			# Set regressor status from 3 to 10 (empty)
			echo "Regressor No. $regressorID is empty, set to status empty (10)"
			sed -i "s/set fmri(shape$regressorID) 3/set fmri(shape$regressorID) 10/g" ${rootdir}/Analyses/fMRI_Scripts/FEAT_GLM${GLMID}_Subject_Scripts/feat_GLM${GLMID}_Subject_Sub${subjectID}.fsf
		fi

	done # end regressor loop

done # end subject loop
