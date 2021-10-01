#!/bin/bash

# Create ROIs as conjunctions of an anatomical mask (from Harvard-Oxford Atlas; potentially combine multiple smaller regions to one bigger one) and selected clusters above zThreshold in GLM1.
# Used for creating volume-by-volume time-series regressors for TAfT and for extract mean PE estimates from selected regions.
#
# Involves the following steps:
# 1) Extract masks from NIFTI of Harvard-Oxford Atlas provided by FSL
# 2) Binarize masks based on probabiliy pThreshold
# 3) Combine smaller masks into bigger masks (like striatum, vmPFC)
# 4) Binarize again
# 5) Binarize pThresholded z-map from GLM1 (select respective contrasts; also do for GLM with SVC)
# 6) Multiply anatomical masks from HOA with pThresholded z-maps of respective contrast
#
# Make executable:
# chmod a+x extract_ROI_HOA_GLM1_TAfT.sh # make executable
#
# Execute:
# ./extract_ROI_HOA_GLM1_TAfT.sh # just run locally, even just single jobs 
#
# EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
# J. Algermissen, 2018-2021.

# Inventory cortical (region numbers): https://neurovault.org/media/images/262/HarvardOxford-cort-maxprob-thr25-2mm.xml
# Inventory subcortical (region numbers): https://neurovault.org/media/images/262/HarvardOxford-sub-maxprob-thr25-1mm.xml

# Set directories:
rootdir=/project/3017042.02 # root directory--needs to be adapted to users' own directory structure
anatdir=$rootdir/Log/fMRI/fMRI_Masks/masksHarvardOxford
targetdir=$rootdir/Log/fMRI/fMRI_Masks/masksTAfTCueLocked

fsldir=/opt/fsl/6.0.0 # FSL's directory--needs to be adapted to users' own directory structure
fslbindir=$fsldir/bin
fsldatadir=$fsldir/data

pThreshold=10 # extract masks based on probabilities > 10%
zThreshold=3.1 # binarize thresholded z-masks at 3.1

# ---------------------------------------------------------------------------------------------------------------- #
## A) Cortex: CingulateAnteriorAction, LeftMotorHand, RightMotorHand:

echo "Extract and cortical masks"

## 1) Extract anatomical masks:
# ACC:
$fslbindir/fslroi $fsldatadir/atlases/HarvardOxford/HarvardOxford-cort-prob-2mm.nii.gz $anatdir/CingulateAnterior 28 1
# JBL (used for mean PE extraction only):
$fslbindir/fslroi $fsldatadir/atlases/HarvardOxford/HarvardOxford-cort-prob-2mm.nii.gz $anatdir/JuxtapositionalLobule 25 1
# Motor cortex:
$fslbindir/fslroi $fsldatadir/atlases/HarvardOxford/HarvardOxford-cort-prob-2mm.nii.gz $anatdir/PrecentralGyrus 6 1
$fslbindir/fslroi $fsldatadir/atlases/HarvardOxford/HarvardOxford-cort-prob-2mm.nii.gz $anatdir/PostcentralGyrus 16 1
# vmPFC:
$fslbindir/fslroi $fsldatadir/atlases/HarvardOxford/HarvardOxford-cort-prob-2mm.nii.gz $anatdir/FrontalPole 0 1
$fslbindir/fslroi $fsldatadir/atlases/HarvardOxford/HarvardOxford-cort-prob-2mm.nii.gz $anatdir/FrontalMedial 24 1
$fslbindir/fslroi $fsldatadir/atlases/HarvardOxford/HarvardOxford-cort-prob-2mm.nii.gz $anatdir/Paracingulate 27 1

## 2) Binarize:
# ACC:
$fslbindir/fslmaths $anatdir/CingulateAnterior -thr $pThreshold -bin $anatdir/CingulateAnterior_bin
# JBL:
$fslbindir/fslmaths $anatdir/JuxtapositionalLobule -thr $pThreshold -bin $anatdir/JuxtapositionalLobule_bin
# Motor cortex:
$fslbindir/fslmaths $anatdir/PrecentralGyrus -thr $pThreshold -bin $anatdir/PrecentralGyrus_bin
$fslbindir/fslmaths $anatdir/PoscentralGyrus -thr $pThreshold -bin $anatdir/PostcentralGyrus_bin
# vmPFC:
$fslbindir/fslmaths $anatdir/FrontalPole -thr $pThreshold -bin $anatdir/FrontalPole_bin
$fslbindir/fslmaths $anatdir/FrontalMedial -thr $pThreshold -bin $anatdir/FrontalMedial_bin
$fslbindir/fslmaths $anatdir/Paracingulate -thr $pThreshold -bin $anatdir/Paracingulate_bin

## 3) Combine to bigger masks:
# Motor cortex:
$fslbindir/fslmaths $anatdir/PrecentralGyrus_bin -add $anatdir/PostcentralGyrus_bin $anatdir/MotorCortex
# vMPFC:
$fslbindir/fslmaths $anatdir/FrontalPole_bin -add $anatdir/FrontalMedial_bin -add $anatdir/Paracingulate_bin $anatdir/vmPFC

## 4) Binarize again:
$fslbindir/fslmaths $anatdir/MotorCortex -thr $pThreshold -bin $anatdir/MotorCortex_bin
$fslbindir/fslmaths $anatdir/vmPFC -thr $pThreshold -bin $anatdir/vmPFC_bin

# ---------------------------------------------------------------------------------------------------------------- #
## B) Striatum:

echo "Extract and binarize striatal masks"

# !) Extract all subparts:
$fslbindir/fslroi $fsldatadir/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm.nii.gz $anatdir/AccumbensLeft 10 1
$fslbindir/fslroi $fsldatadir/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm.nii.gz $anatdir/AccumbensRight 20 1
$fslbindir/fslroi $fsldatadir/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm.nii.gz $anatdir/CaudateLeft 4 1
$fslbindir/fslroi $fsldatadir/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm.nii.gz $anatdir/CaudateRight 15 1
$fslbindir/fslroi $fsldatadir/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm.nii.gz $anatdir/PutamenLeft 5 1
$fslbindir/fslroi $fsldatadir/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm.nii.gz $anatdir/PutamenRight 16 1

# 2) Binarize:
$fslbindir/fslmaths $anatdir/AccumbensLeft -thr $pThreshold -bin $anatdir/AccumbensLeft_bin
$fslbindir/fslmaths $anatdir/AccumbensRight -thr $pThreshold -bin $anatdir/AccumbensRight_bin
$fslbindir/fslmaths $anatdir/CaudateLeft -thr $pThreshold -bin $anatdir/CaudateLeft_bin
$fslbindir/fslmaths $anatdir/CaudateRight -thr $pThreshold -bin $anatdir/CaudateRight_bin
$fslbindir/fslmaths $anatdir/PutamenLeft -thr $pThreshold -bin $anatdir/PutamenLeft_bin
$fslbindir/fslmaths $anatdir/PutamenRight -thr $pThreshold -bin $anatdir/PutamenRight_bin

# 3) Combine to entire bilateral striatum and accumbens masks:
$fslbindir/fslmaths $anatdir/AccumbensLeft_bin -add $anatdir/AccumbensRight_bin -add $anatdir/CaudateLeft_bin -add $anatdir/CaudateRight_bin -add $anatdir/PutamenLeft_bin -add $anatdir/PutamenRight_bin $anatdir/Striatum
$fslbindir/fslmaths $anatdir/AccumbensLeft_bin -add $anatdir/AccumbensRight_bin $anatdir/Accumbens

# 4) Binarize again:
$fslbindir/fslmaths $anatdir/Striatum -bin $anatdir/Striatum_bin 
$fslbindir/fslmaths $anatdir/Accumbens -bin $anatdir/Accumbens_bin

# ---------------------------------------------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------------------------------------------- #
# ---------------------------------------------------------------------------------------------------------------- #
# 5) Binarize group-level z-maps from GLM1:
echo "Binarize group-level masks at $zThreshold"

# a) Valence positive: to be multiplied with vmPFC mask:
$fslbindir/fslmaths ${rootdir}/Log/fMRI/GLM1_FEAT_Combined_pos.gfeat/cope1.feat/thresh_zstat1 -thr $zThreshold -bin ${rootdir}/Log/fMRI/GLM1_FEAT_Combined_pos.gfeat/cope1.feat/thresh_zstat1_bin

# Also for GLM with SVC with striatal mask to get posterior putamen:
$fslbindir/fslmaths ${rootdir}/Log/fMRI/GLM1_FEAT_Combined_pos_striatum.gfeat/cope1.feat/thresh_zstat1 -thr $zThreshold -bin ${rootdir}/Log/fMRI/GLM1_FEAT_Combined_pos_striatum.gfeat/cope1.feat/thresh_zstat1_bin

# b) Valence negative: to be multiplied with ACC mask:
$fslbindir/fslmaths ${rootdir}/Log/fMRI/GLM1_FEAT_Combined_neg.gfeat/cope1.feat/thresh_zstat1 -thr $zThreshold -bin ${rootdir}/Log/fMRI/GLM1_FEAT_Combined_neg.gfeat/cope1.feat/thresh_zstat1_bin

# Also for GLM with SVC with striatal mask to get medial caudate:
$fslbindir/fslmaths ${rootdir}/Log/fMRI/GLM1_FEAT_Combined_neg_striatum.gfeat/cope1.feat/thresh_zstat1 -thr $zThreshold -bin ${rootdir}/Log/fMRI/GLM1_FEAT_Combined_neg_striatum.gfeat/cope1.feat/thresh_zstat1_bin

# c) Action positive: to be multiplied with striatum and ACC masks:
$fslbindir/fslmaths ${rootdir}/Log/fMRI/GLM1_FEAT_Combined_pos.gfeat/cope2.feat/thresh_zstat1 -thr $zThreshold -bin ${rootdir}/Log/fMRI/GLM1_FEAT_Combined_pos.gfeat/cope2.feat/thresh_zstat1_bin

# d) Conflict: to be multipled with JBL (only used for mean PE extract):
$fslbindir/fslmaths ${rootdir}/Log/fMRI/GL1_FEAT_Combined_pos_ACC_JBL.gfeat/cope3.feat/thresh_zstat1 -thr $zThreshold -bin ${rootdir}/Log/fMRI/GLM1_FEAT_Combined_pos_ACC_JBL.gfeat/cope3.feat/thresh_zstat1_bin

# e) Hand left: to be multiplied with right motor cortex mask:
$fslbindir/fslmaths ${rootdir}/Log/fMRI/GLM1_FEAT_Combined_pos.gfeat/cope4.feat/thresh_zstat1 -thr $zThreshold -bin ${rootdir}/Log/fMRI/GLM1_FEAT_Combined_pos.gfeat/cope4.feat/thresh_zstat1_bin

# f) Hand right: to be multiplied with left motor cortex mask:
$fslbindir/fslmaths ${rootdir}/Log/fMRI/GLM1_FEAT_Combined_neg.gfeat/cope4.feat/thresh_zstat1 -thr $zThreshold -bin ${rootdir}/Log/fMRI/GLM1_FEAT_Combined_neg.gfeat/cope4.feat/thresh_zstat1_bin

# ---------------------------------------------------------------------------------------------------------------- #
# 6) Multiply anatomical and group-level mask:
echo "Multiply anatomical and group masks"

# ----------------------------------------------------------------------------- #
## a) Valence contrast:

# vmPFCValence (adjust manually, see at the bottom of this script):
$fslbindir/fslmaths ${rootdir}/Log/fMRI/GLM1_FEAT_Combined_pos.gfeat/cope1.feat/thresh_zstat1_bin -mul $anatdir/vmPFC_bin -bin $targetdir/GLM1vmPFCValence.nii.gz

# ACCValence:
$fslbindir/fslmaths ${rootdir}/Log/fMRI/GLM1_FEAT_Combined_neg.gfeat/cope1.feat/thresh_zstat1_bin -mul $anatdir/CingulateAnterior_bin -bin $targetdir/GLM1CingulateAnteriorValence.nii.gz

# LeftPutamenValence:
$fslbindir/fslmaths ${rootdir}/Log/fMRI/GLM1_FEAT_Combined_pos_striatum.gfeat/cope1.feat/thresh_zstat1_bin -mul $anatdir/Striatum_bin -bin $targetdir/GLM1LeftPutamenValence.nii.gz

# MedialCaudateValence:
$fslbindir/fslmaths ${rootdir}/Log/fMRI/GLM1_FEAT_Combined_neg_striatum.gfeat/cope1.feat/thresh_zstat1_bin -mul $anatdir/Striatum_bin -bin $targetdir/GLM1MedialCaudateValence.nii.gz

# ----------------------------------------------------------------------------- #
## b) Action contrast:

# StriatumAction:
$fslbindir/fslmaths ${rootdir}/Log/fMRI/GLM1_FEAT_Combined_pos.gfeat/cope2.feat/thresh_zstat1_bin -mul $anatdir/Striatum_bin -bin $targetdir/GLM1StriatumAction.nii.gz

# ACCAction:
$fslbindir/fslmaths ${rootdir}/Log/fMRI/GLM1_FEAT_Combined_pos.gfeat/cope2.feat/thresh_zstat1_bin -mul $anatdir/CingulateAnterior_bin -bin $targetdir/GLM1CingulateAnteriorAction.nii.gz

# ----------------------------------------------------------------------------- #
## c) Conflict contrast:

$fslbindir/fslmaths ${rootdir}/Log/fMRI/GLM1_FEAT_Combined_pos_ACC_JBL.gfeat/cope3.feat/thresh_zstat1_bin -mul $anatdir/JuxtapositionalLobule_bin -bin $targetdir/GLM1JBLIncongruency.nii.gz

# ----------------------------------------------------------------------------- #
## d) Left/ right hand contrast:

# RightMotorHand:
$fslbindir/fslmaths ${rootdir}/Log/fMRI/GLM1_FEAT_Combined_pos.gfeat/cope4.feat/thresh_zstat1_bin -mul $anatdir/MotorCortex_bin -bin $targetdir/GLM1RightMotorHand.nii.gz

# LeftMotorHand:
$fslbindir/fslmaths ${rootdir}/Log/fMRI/GLM1_FEAT_Combined_neg.gfeat/cope4.feat/thresh_zstat1_bin -mul $anatdir/MotorCortex_bin -bin $targetdir/GLM1LeftMotorHand.nii.gz

echo "Finished :-)"

# e) Adjust manually in FSLeyes --> GLM1vmPFCValenceMan.nii.gz
# gunzip ${rootdir}/Log/fMRI/fMRI_Masks/masksTAfTCueLocked/GLM1vmPFCValenceMan.nii.gz
# END
