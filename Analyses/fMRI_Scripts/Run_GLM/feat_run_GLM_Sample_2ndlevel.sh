#!/bin/bash

# Runs second-level GLM one sample level based on .fsf file (see respective "FEAT_GLM${GLMID}_Sample_Scripts" folder)
#
# Make executable:
# chmod a+x feat_run_GLM_Sample_2ndlevel.sh
#
# ./feat_run_GLM_Sample_2ndlevel.sh # run
#use fsl_sub_DCCN.sh (in my home directory), run on -T, minutes, storage space, feat command, file to run (with complete path)
#
# EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
# J. Algermissen, 2018-2021.

rootdir=/project/3017042.02 # root directory--needs to be adapted to users' folder structure
submitdir=${rootdir}/Analyses/fMRI_Scripts/Run_GLM # where fsl_sub_DCCN.sh is

GLMID=1

# NOTE: Because this is second-level, you can fit all contrasts in one single command. No need to fit separate GLMs per contrast

# Regular GLMs:
$submitdir/fsl_sub_DCCN.sh -N "GLM${GLMID}_pos" -T 240 -u 10gb feat ${rootdir}/Analyses/fMRI_Scripts/FEAT_GLM${GLMID}_Sample_Scripts/feat_GLM${GLMID}_Sample_Scripts_pos.fsf;
$submitdir/fsl_sub_DCCN.sh -N "GLM${GLMID}_neg" -T 240 -u 10gb feat ${rootdir}/Analyses/fMRI_Scripts/FEAT_GLM${GLMID}_Sample_Scripts/feat_GLM${GLMID}_Sample_Scripts_neg.fsf;

# -----------------------------------------------------------------------------------
# With SVC for mask of striatum:
$submitdir/fsl_sub_DCCN.sh -N "GLM${GLMID}_pos_striatum" -T 240 -u 10gb feat ${rootdir}/Analyses/fMRI_Scripts/FEAT_GLM${GLMID}_Sample_Scripts/feat_GLM${GLMID}_Sample_Scripts_pos_striatum.fsf;
$submitdir/fsl_sub_DCCN.sh -N "GLM${GLMID}_neg_striatum" -T 240 -u 10gb feat ${rootdir}/Analyses/fMRI_Scripts/FEAT_GLM${GLMID}_Sample_Scripts/feat_GLM${GLMID}_Sample_Scripts_neg_striatum.fsf;

# -----------------------------------------------------------------------------------
# With SVC for mask of ACC_JBL:
$submitdir/fsl_sub_DCCN.sh -N "GLM${GLMID}_pos_ACC_JBL" -T 240 -u 10gb feat ${rootdir}/Analyses/fMRI_Scripts/FEAT_GLM${GLMID}_Sample_Scripts/feat_GLM${GLMID}_Sample_Scripts_pos_ACC_JBL.fsf;

# -----------------------------------------------------------------------------------
# Control with only 30 subjects (included in TAfT):
$submitdir/fsl_sub_DCCN.sh -N "GLM${GLMID}_pos" -T 240 -u 10gb feat ${rootdir}/Analyses/fMRI_Scripts/FEAT_GLM${GLMID}_Sample_Scripts/feat_GLM${GLMID}_without6_Sample_Scripts_pos.fsf;
$submitdir/fsl_sub_DCCN.sh -N "GLM${GLMID}_neg" -T 240 -u 10gb feat ${rootdir}/Analyses/fMRI_Scripts/FEAT_GLM${GLMID}_Sample_Scripts/feat_GLM${GLMID}_without6_Sample_Scripts_neg.fsf;
$submitdir/fsl_sub_DCCN.sh -N "GLM${GLMID}_pos_striatum" -T 240 -u 10gb feat ${rootdir}/Analyses/fMRI_Scripts/FEAT_GLM${GLMID}_Sample_Scripts/feat_GLM${GLMID}_without6_Sample_Scripts_pos_striatum.fsf;
$submitdir/fsl_sub_DCCN.sh -N "GLM${GLMID}_neg_striatum" -T 240 -u 10gb feat ${rootdir}/Analyses/fMRI_Scripts/FEAT_GLM${GLMID}_Sample_Scripts/feat_GLM${GLMID}_without6_Sample_Scripts_neg_striatum.fsf;
$submitdir/fsl_sub_DCCN.sh -N "GLM${GLMID}_pos_ACC_JBL" -T 240 -u 10gb feat ${rootdir}/Analyses/fMRI_Scripts/FEAT_GLM${GLMID}_Sample_Scripts/feat_GLM${GLMID}_without6_Sample_Scripts_pos_ACC_JBL.fsf;

# END
