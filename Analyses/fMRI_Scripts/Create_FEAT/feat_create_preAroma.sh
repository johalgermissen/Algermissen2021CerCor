#!/bin/bash

# Automatically creates .fsf files per subject per block for preprocessing based on template feat_template_preAroma.fsf.
# 
# Make executable:
# chmod a+x feat_create_preAroma.sh
#
# Submit as job to cluster:
# qsub -N "feat_create_preAroma" -l walltime=0:02:00,mem=1mb feat_create_preAroma.sh
#
# EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
# J. Algermissen, 2018-2021.

rootdir=/project/3017042.02 # root directory--needs to be adapted to users' folder structure

# Set subject ID, loop over subjects:
for (( subject=1 ; subject<=36 ; subject++ )); do

	subjectID=`zeropad $subject 3` # subject ID with 3 digits

	# Set block ID, loop over blocks:
	for (( block=1 ; block<=6 ; block++ ));	do

		blockID=`zeropad $block 1` # block ID with 1 digit

		# Retrieve number of volumes for this block: use fslinfo, first five rows, last row of these, print number of volumes to variable npts
		npts=`fslinfo ${rootdir}/Log/fMRI/sub-${subjectID}/raw/sub-${subjectID}_block${blockID}.nii | head -n 5 | tail -n 1 | awk '{print $2}'`
		
		# Specify template, number of time points, number of blocks, output variable
		cat ${rootdir}/Analyses/fMRI_Scripts/FEAT_Templates/feat_template_preAroma.fsf | sed s/TIMEPOINTS/$npts/g | sed s/SUBJECT/$subjectID/g | sed s/BLOCK/$blockID/g > ${rootdir}/Analyses/fMRI_Scripts/FEAT_preAroma_Scripts/feat_Sub${subjectID}_Block${blockID}_preAroma.fsf

	done # end of block loop

done # end of subject loop
