#!/bin/bash

# FigureS05.sh

# Plots for Figure S05 in supplementary material.
# Mind adjusting root directory in analyze_BOLD_effect_over_time.R.
#
# EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
# J. Algermissen, 2018-2021.

# make executable:
# chmod a+x FigureS05.sh # make executable
# submit this job with:
# qsub -N "FigureS05" -l walltime=30:00:00,mem=2gb FigureS05.sh

module unload R
module load R/4.1.0

rootdir=/project/3017042.02

modID=10
echo "Start script"
R CMD BATCH $rootdir/Analyses/fMRI_Scripts/Run_ROI/analyze_BOLD_effect_over_time.R
echo "Start script"


