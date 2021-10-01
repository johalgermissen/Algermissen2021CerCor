#' analyze_BOLD_effect_over_time.R
#' 
#' Load trial-by-trial deconvolved BOLD signal in selected regions (conjunctions anatomical & functional contrast).
#' Perform mixed-effects regression with trial number and cognitive variable as predictors.
#' Make bar plots with "task block" on x-axis, BOLD signal on y-axis, and cognitive variable as color.
#' Interactive script, run step-by-step, NOT in one Go.
#' 
#' Requires that .csv files have been recoded and concatenated using 
#' EEGfMRIPav_1_Read_Initial.R
#' 
#' Calls EEGfMRIPav_2_Preprocess_Automated.R
#' 
#' Calls EEGfMRIPav_0_CreateBarPlots.R to create bar plots with 2 IVs
#'
#' Assumes that trial-by-trial deconvolved BOLD signal as be created and saved per subject under paste0(rootdir,"Log/EEG/CueLockedResults/TAfT_Betas/...)
#'
#' Mind to adjust rootdir.
#'
#' INPUTS:
#' none
#'
#' OUTPUTS:
#' no outputs, fits (and saves) models.
#'
#' EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
#' J. Algermissen, 2018-2021.

# ====================================================================================================================== #
#### Load packages ####

cat("Load packages\n")
library(ggplot2)
library(plyr)
library(stringr)

# Analyses:
library(lme4)
library(afex)

options(scipen = 20)
options(contrasts = c("contr.sum", "contr.poly"))
rm(list=ls())
# ================================================================================================================================================ #
#### 1) Load full data set ####

## Set directories:
rootdir   <- "/project/3017042.02/"
scriptdir <- paste0(rootdir,"Analyses/Behavior_Scripts")
plotdir <- paste0(rootdir,"Log/Behavior/Plots_R/")

## Load functions:
setwd(scriptdir)
source("EEGfMRIPav_1_Read_Initial.R"); # Load variables and bring to R data frame format
source("EEGfMRIPav_2_Preprocess_Automated.R"); # Recode variables
source("EEGfMRIPav_2_Cue_Position.R") # Retrieve cue position
source("EEGfMRIPav_0_CreateBarPlots.R") # plotting functions

## Apply functions:
# mydata <- EEGfMRIPav_1_Read_Initial() # Initial reading, converting to factors 
mydata <- EEGfMRIPav_2_Preprocess_Automated() # Recode variables based on uninstructed key presses etc.
mydata <- EEGfMRIPav_2_Cue_Position(mydata) # Compute cue position
# mydata$cue_pos

# setwd(moddir) # where to store models

# ================================================================================================================================================ #
#### 2a) Load outcome-locked BOLD trial-by-trial data: ####

nSub <- length(unique(mydata$PPN_n))
# nSub <- 36

## Load all variables:
allVarNames <- c("GLM5FvmPFCValenceMan","GLM5FCingulateAnteriorValence","GLM5FLeftPutamenValence","GLM5FMedialCaudateValence","GLM5FHippocampusValence","GLM5FStriatumAction","GLM5FCingulateAnteriorAction","GLM5FJBLIncongruency")
## Set directory for respective ROI:
for (iVar in 1:length(allVarNames)){
  varName <- allVarNames[iVar]
  fmridir <- paste0(rootdir,"Log/EEG/CueLockedResults/TAfT_Betas/TAfT_BOLD_",varName,"/")
  ## Load:
  mydata[,varName] <- NA # initialize
  ## Loop over subject files:
  for (iSub in 1:nSub){ # iSub <- 1
    cat(paste0("Start subject ",iSub, "\n"))
    fileName <- paste0(fmridir,"TAfT_BOLD_",varName,"_sub",str_pad(iSub,3,side = "left",pad = 0),".csv")
    tmp <- as.matrix(read.csv(fileName, header = F))
    if (length(tmp) != 640){stop(paste0("Subject ",iSub,": Input not of length 640 trials"))}
    subIdx <- which(mydata$PPN_n==iSub) 
    mydata[subIdx,varName] <- tmp
  }
}
cat("Finished :-)\n")
# names(mydata) # check out names of new variables

# ================================================================================================================================================ #
#### 2b) Select subjects: #####

invalidSubs <- c(15,25) # failed registration
validSubs <- which(!(1:nSub %in% invalidSubs))
modData <- subset(mydata, PPN_n %in% validSubs)

# ================================================================================================================================================ #
#### 3) Mixed effects model: trial number x cognitive variable of interest ####

## Trial number cross blocks: 
modData$trial_nr_all_z <- as.numeric(scale(modData$trial_nr_all_n))

## Select respective ROI (conjunction anatomical region and function contrast):
# Valence:
varName <- "GLM5FvmPFCValenceMan"
varName <- "GLM5FCingulateAnteriorValence"
varName <- "GLM5FLeftPutamenValence"
varName <- "GLM5FMedialCaudateValence"

# Action:
varName <- "GLM5FStriatumAction"
varName <- "GLM5FCingulateAnteriorAction"

# Conflict:
varName <- "GLM5FJBLIncongruency"

## Select region & contrast, copy to variable fMRIPRed:
modData$fMRIPred <- modData[,varName]
modData$fMRIPred_z <- as.numeric(scale(modData$fMRIPred))

## Select formula based on contrast of interest:

# Valence only:
formula <- "fMRIPred_z ~ stim_win_f * trial_nr_all_z + (stim_win_f * trial_nr_all_z|PPN_f)"
# Action only:
formula <- "fMRIPred_z ~ is_go_cor_f * trial_nr_all_z + (is_go_cor_f * trial_nr_all_z|PPN_f)"
# Conflict only:
formula <- "fMRIPred_z ~ conflict_f * trial_nr_all_z + (conflict_f * trial_nr_all_z|PPN_f)"

## Fit model:
mod <- lmer(formula = formula, modData,
            control = lmerControl(optCtrl = list(maxfun = 1e+9), calc.derivs = F, optimizer = c("bobyqa")))

## Inspect model output:
summary(mod)

## P-values with LRTs:
mod_LRT <- mixed(mod, modData, method = "LRT", type = "III",
                 control = lmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)

# ================================================================================================================================================ #
#### 4) Bar plots: block number (x-axis) x cognitive variable of interest (color) ####
# all saved as pngs

## 1) Loop for valence:
for (iVar in 1:length(allVarNames)){ # iVar <- 1
  selVar <- allVarNames[iVar]
  custom_barplot2(modData, xVar="block_n", yVar=selVar, zVar="stim_win_f", subVar="PPN_f",
                  xLab = "Block number", yLab = "BOLD signal (a.u.)", zLab = "Valence", selCol = c("#CC0000","#009933"),
                  isPoint = F, savePNG = T, saveEPS = F)
  # With dots:
  custom_barplot2(modData, xVar="block_n", yVar=selVar, zVar="stim_win_f", subVar="PPN_f",
                  xLab = "Block number", yLab = "BOLD signal (a.u.)", zLab = "Valence", selCol = c("#CC0000","#009933"),
                  isPoint = T, savePNG = T, saveEPS = F)
}

## 2) Loop for action:
for (iVar in 1:length(allVarNames)){
  selVar <- allVarNames[iVar]
  custom_barplot2(modData, xVar="block_n", yVar=selVar, zVar="is_go_cor_f", subVar="PPN_f",
                  xLab = "Block number", yLab = "BOLD signal (a.u.)", zLab = "Performed action", selCol = c("blue","red"),
                  isPoint = F, savePNG = T, saveEPS = F)
  # With dots:
  custom_barplot2(modData, xVar="block_n", yVar=selVar, zVar="is_go_cor_f", subVar="PPN_f",
                  xLab = "Block number", yLab = "BOLD signal (a.u.)", zLab = "Performed action", selCol = c("blue","red"),
                  isPoint = T, savePNG = T, saveEPS = F)
}

## 3) Loop for congruency:
for (iVar in 1:length(allVarNames)){
  selVar <- allVarNames[iVar]
  custom_barplot2(modData, xVar="block_n", yVar=selVar, zVar="conflict_f", subVar="PPN_f",
                  xLab = "Block number", yLab = "BOLD signal (a.u.)", zLab = "Congruency", selCol = c("#F5BD1F","#5F2689"),
                  isPoint = F, savePNG = T, saveEPS = F)
  # With dots:
  custom_barplot2(modData, xVar="block_n", yVar=selVar, zVar="conflict_f", subVar="PPN_f",
                  xLab = "Block number", yLab = "BOLD signal (a.u.)", zLab = "Congruency", selCol = c("#F5BD1F","#5F2689"),
                  isPoint = T, savePNG = T, saveEPS = F)
}

# END