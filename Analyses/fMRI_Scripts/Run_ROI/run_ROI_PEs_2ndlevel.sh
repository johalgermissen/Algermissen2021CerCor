#!/bin/bash

# Loop over subjects, extract summary statistics of parameter estimates within given mask using featquery.
# Steps involved:
# 1) Run invert_highres2standard_warp_subject.sh first to invert FNIRT registration matrix.
# 2) Brings given mask from standard to functional space for each subject.
# 3) Uses featquery to get summary statistics from parameter estimates within given mask for each subject.
# 
# Make executable:
# chmod a+x run_ROI_PEs_2ndlevel.sh # make executable
#
# Execute:
# ./run_ROI_PEs_2ndlevel.sh # run
# takes about 7-8 minutes per subject, so 4-5 h in total
# qsub -N "HippocampusValence" -l walltime=5:00:00,mem=10gb /project/3017042.02/Analyses/fMRI_Scripts/Run_ROI/run_ROI_PEs_2ndlevel.sh
#
# EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
# J. Algermissen, 2018-2021.

# Follow Andrew Jahn's video under: https://www.youtube.com/watch?v=SNVt7smHLm8
# Perform on single-subject level, thus loop over subjects
# Extract mean parameter estimates and analyse in R or Matlab or whatever

# Featquery options:
# -p convert PE/COPE values into % 
# -s create time-series plots
# -b popup results in browser when finished (rather disable)

# Settings:

rootdir=/project/3017042.02 # root directory--needs to be adapted to users' own directory structure

fsldir=/opt/fsl/6.0.0/bin # FSL's directory--needs to be adapted to users' own directory structure

# GLM to extract from:
GLMIDextract=1

# a) Mask from 2nd-level GLM:
GLMIDmask=1

# Select mask for which to extract mean PEs: 
maskID=1${GLMIDmask}MedialCaudateValence
#maskID=1${GLMIDmask}LeftPutamenValence
#maskID=1${GLMIDmask}JBLIncongruency

## File of mask:
maskfile=${rootdir}/Log/fMRI/fMRI_Masks/masksTAfTCueLocked/${maskID}.nii.gz # functional masks from cue-locked GLMs

# set subject ID, loop
for (( subject=1 ; subject<=36; subject++ )); do

	# Extract for each subject; mind: registration to MNI space already applied at first level; so no inverse registration needed
	subjectID=`zeropad $subject 3` # subject ID with 3 digits
	echo "Subject ${subjectID}: Start"

	# 1) Invert FNIRT registration matrix per subject:
	# run invert_highres2standard_warp_subject.sh first !!!!!

	# 2) Bring mask from standard into functional space:
	echo "Subject ${subjectID}: Bring mask from standard to functinoal space"
	$fsldir/applywarp -i ${maskfile} -r ${rootdir}/Log/fMRI/sub-${subjectID}/GLM${GLMIDextract}/GLM${GLMIDextract}_sub${subjectID}.feat/example_func -o ${rootdir}/Log/fMRI/sub-${subjectID}/GLM${GLMIDextract}/GLM${GLMIDextract}_sub${subjectID}.feat/${maskID}_func -w ${rootdir}/Log/fMRI/sub-${subjectID}/GLM${GLMIDextract}/GLM${GLMIDextract}_sub${subjectID}.feat/reg/standard2highres_warp --postmat=${rootdir}/Log/fMRI/sub-${subjectID}/GLM${GLMIDextract}/GLM${GLMIDextract}_sub${subjectID}.feat/reg/highres2example_func.mat

	# 3) Apply functional-space mask to PEs (first 4 PEs: pe1 pe2 pe3 pe4):
	echo "Subject ${subjectID}: Extract PEs"

	$fsldir/featquery 1 ${rootdir}/Log/fMRI/sub-${subjectID}/GLM${GLMIDextract}/GLM${GLMIDextract}_sub${subjectID}.feat/ 4 stats/pe1 stats/pe2 stats/pe3 stats/pe4 featquery_${maskID} -a 9 -p ${rootdir}/Log/fMRI/sub-${subjectID}/GLM${GLMIDextract}/GLM${GLMIDextract}_sub${subjectID}.feat/${maskID}_func

	echo "Subject ${subjectID}: Finished"

done # end subject loop
