#!/bin/bash

# Create an ROI as a conjunction of an anatomical mask (from Harvard-Oxford Atlas; potentially combine multiple smaller regions to one bigger one) and selected clusters above pThreshold in GLM1.
# Involves the following steps:
# 1) Extract masks from NIFTI of Harvard-Oxford Atlas provided by FSL
# 2) Binarize masks based on probabiliy pThreshold
# 3) Combine smaller masks into bigger masks (like striatum, vmPFC)
# 4) Binarize again
#
# Make executable:
# chmod a+x extract_ROI_HOA_GLM1_SVC.sh
#
# Execute:
# ./extract_ROI_HOA_GLM1_SVC.sh # just run locally 
#
# EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
# J. Algermissen, 2018-2021.

# Inventory cortical (region numbers): https://neurovault.org/media/images/262/HarvardOxford-cort-maxprob-thr25-2mm.xml
# Inventory subcortical (region numbers): https://neurovault.org/media/images/262/HarvardOxford-sub-maxprob-thr25-1mm.xml

# Set directories:
rootdir=/project/3017042.02 # root directory--needs to be adapted to users' own directory structure
targetdir=$rootdir/Log/fMRI/fMRI_Masks/masksHarvardOxford

fsldir=/opt/fsl/6.0.0 # FSL's directory--needs to be adapted to users' own directory
fslbindir=$fsldir/bin
fsldatadir=$fsldir/data

pThreshold=10 # extract masks based on probabilities > 10%

# ---------------------------------------------------------------------------------------------------------------- #
## 1) Extract masks from Harvard-Oxford Atlas:
echo "Extract masks from Harvard-Oxford Atlas"

# ACC:
$fslbindir/fslroi $fsldatadir/atlases/HarvardOxford/HarvardOxford-cort-prob-2mm.nii.gz $targetdir/CingulateAnterior 28 1
# JBL:
$fslbindir/fslroi $fsldatadir/atlases/HarvardOxford/HarvardOxford-cort-prob-2mm.nii.gz $targetdir/JuxtapositionalLobule 25 1

# Subparts of striatum:
$fslbindir/fslroi $fsldatadir/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm.nii.gz $targetdir/AccumbensLeft 10 1
$fslbindir/fslroi $fsldatadir/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm.nii.gz $targetdir/AccumbensRight 20 1
$fslbindir/fslroi $fsldatadir/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm.nii.gz $targetdir/CaudateLeft 4 1
$fslbindir/fslroi $fsldatadir/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm.nii.gz $targetdir/CaudateRight 15 1
$fslbindir/fslroi $fsldatadir/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm.nii.gz $targetdir/PutamenLeft 5 1
$fslbindir/fslroi $fsldatadir/atlases/HarvardOxford/HarvardOxford-sub-prob-2mm.nii.gz $targetdir/PutamenRight 16 1

# ---------------------------------------------------------------------------------------------------------------- #
# 2) Binarize:
echo "Binarize masks at threshold $pThreshold"

$fslbindir/fslmaths $targetdir/CingulateAnterior -thr $pThreshold -bin $targetdir/CingulateAnterior_bin
$fslbindir/fslmaths $targetdir/JuxtapositionalLobule -thr $pThreshold -bin $targetdir/JuxtapositionalLobule_bin

$fslbindir/fslmaths $targetdir/AccumbensLeft -thr $pThreshold -bin $targetdir/AccumbensLeft_bin
$fslbindir/fslmaths $targetdir/AccumbensRight -thr $pThreshold -bin $targetdir/AccumbensRight_bin
$fslbindir/fslmaths $targetdir/CaudateLeft -thr $pThreshold -bin $targetdir/CaudateLeft_bin
$fslbindir/fslmaths $targetdir/CaudateRight -thr $pThreshold -bin $targetdir/CaudateRight_bin
$fslbindir/fslmaths $targetdir/PutamenLeft -thr $pThreshold -bin $targetdir/PutamenLeft_bin
$fslbindir/fslmaths $targetdir/PutamenRight -thr $pThreshold -bin $targetdir/PutamenRight_bin

# ---------------------------------------------------------------------------------------------------------------- #
# 3) Combine ACC and JBL to entire bilateral striatum:
echo "Combine smaller masks into bigger masks"

$fslbindir/fslmaths $targetdir/CingulateAnterior_bin -add $targetdir/JuxtapositionalLobule_bin $targetdir/ACC_JBL

$fslbindir/fslmaths $targetdir/AccumbensLeft_bin -add $targetdir/AccumbensRight_bin -add $$targetdir/CaudateLeft_bin -add $targetdir/CaudateRight_bin -add $targetdir/PutamenLeft_bin -add $targetdir/PutamenRight_bin $targetdir/Striatum

# ---------------------------------------------------------------------------------------------------------------- #
# 4) Binarize again:
echo "Binarize again"
$fslbindir/fslmaths $targetdir/Striatum -bin $targedir/Striatum_bin 

# Binarize again:
$fslbindir/fslmaths $targetdir/ACC_JBL -bin $targetdir/ACC_JBL_bin

echo "Finished :-)"
# END
