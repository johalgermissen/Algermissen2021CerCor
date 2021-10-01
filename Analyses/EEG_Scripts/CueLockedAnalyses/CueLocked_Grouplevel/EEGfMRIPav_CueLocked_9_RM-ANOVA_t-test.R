#' EEGfMRIPav_CueLocked_9_RM-ANOVA_t-test.R
#' 
#' Perform RM-ANOVAs and t-tests on data exported from Matlab as .csv files.
#' Interactive script, run step-by-step, NOT in one Go.
#' 
#' Requires that .csv has been exported using 
#' EEGfMRIPav_CueLocked_8_TF_grouplevel_test.m
#'
#' Mind to adjust rootdir.
#'
#' INPUTS:
#' none
#'
#' OUTPUTS:
#' no outputs, fits and prints models.
#'
#' EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
#' J. Algermissen, 2018-2021.
#' 

# ====================================================================================================================== #
#### Load packages ####

cat("Load packages\n")
# Effect size:
library(effsize)
# RM-ANOVA:
library(ez)

options(scipen = 20)
rm(list=ls())

# ====================================================================================================================== #
#### Set options, delete existing objects: ####

cat("Penalize scientific notation\n")
options(scipen = 20) # scientific notation: penalize more heavily

# rm(list=ls()) # delete all previously loaded objects

# ====================================================================================================================== #
#### Specify directories: ####

cat("Set directories\n")
rootdir   <- "/project/3017042.02/" # root directory--needs to be adapted to users' folder structure
inputdir <- paste0(rootdir,"Log/EEG/CueLockedResults/TF_grouplevel/")

# ============================================================================================================== #
#### Load data: ####

mydata_wide <- read.csv(paste0(inputdir,
                               "TFSubCond_Congruency_stimlocked_baseCor_trend_baseTime_0_0_Go_exeAct_bothAcc_midfrontal_alpha.csv"),
                        header = F)
names(mydata_wide) <- c("G2WCor","G2ACor","NG2WCor","NG2ACor","G2WIncor","G2AIncor","NG2WIncor","NG2AIncor")
mydata_wide$PPN_n <- 1:nrow(mydata_wide) # add subject number

# ============================================================================================================== #
#### Cast into long format: ####

mydata_long <- reshape(mydata_wide, idvar = 'PPN_n', 
                       varying = c("G2WCor","G2ACor","NG2WCor","NG2ACor","G2WIncor","G2AIncor","NG2WIncor","NG2AIncor"),
                       timevar = 'condition', v.names = "score", direction = "long")

## Recode:
mydata_long$valence_f <- factor(ifelse(mydata_long$condition %in% c(1,3,5,7), 1, 2), levels = 1:2, labels = c("Win","Avoid"))
mydata_long$action_f <- factor(ifelse(mydata_long$condition %in% c(1,2,5,6), 1, 2), levels = 1:2, labels = c("Go","NoGo"))
mydata_long$accuracy_f <- factor(ifelse(mydata_long$condition %in% c(1,2,3,4), 1, 2), levels = 1:2, labels = c("Correct","Incorrect"))
mydata_long$condition_f <- factor(mydata_long$condition, levels = 1:8, labels = c("G2WCor","G2ACor","NG2WCor","NG2ACor","G2WIncor","G2AIncor","NG2WIncor","NG2AIncor"))

## Subset: 
mydata_long_cor <- subset(mydata_long, accuracy_f == "Correct")
mydata_long_incor <- subset(mydata_long, accuracy_f == "Incorrect")

# ============================================================================================================== #
#### RM-ANOVA: ####

## All conditions:
ezANOVA(data = mydata_long, dv = .(score), wid = .(PPN_n),
        within = .(valence_f, action_f, accuracy_f), type = 3, 
        detailed = TRUE)
#                          Effect DFn DFd           SSn       SSd           F                            p p<.05          ges
# 1                   (Intercept)   1  35 3287.46392738 302.02352 380.9678110 0.00000000000000000002150124     * 0.8868858394
# 2                     valence_f   1  35    3.06703176  15.49091   6.9296188 0.01253030222754369232829319     * 0.0072617746
# 3                      action_f   1  35    0.07665603  16.01653   0.1675120 0.68482646987421635920156859       0.0001827918
# 4                    accuracy_f   1  35    1.12337184  19.42457   2.0241384 0.16366905017149108170215754       0.0026720911
# 5            valence_f:action_f   1  35    3.47010764  14.51335   8.3684145 0.00652632495885600329155185     * 0.0082082989
# 6          valence_f:accuracy_f   1  35    0.07260698  14.02753   0.1811612 0.67298314965304717460981010       0.0001731382
# 7           action_f:accuracy_f   1  35    1.69127406  22.23530   2.6621902 0.11172594644274365383029846       0.0040174959
# 8 valence_f:action_f:accuracy_f   1  35    2.26775248  15.55419   5.1028926 0.03022009589888973238314129     * 0.0053795110

### Correct only:
ezANOVA(data = mydata_long_cor, dv = .(score), wid = .(PPN_n), 
        within = .(valence_f, action_f), type = 3, 
        detailed = TRUE)
#               Effect DFn DFd         SSn       SSd          F                           p p<.05         ges
# 1        (Intercept)   1  35 1583.523226 170.01579 325.989208 0.0000000000000000002586664     * 0.866666774
# 2          valence_f   1  35    1.097921  21.24479   1.808784 0.1873049314752524563409963       0.004486500
# 3           action_f   1  35    1.244029  31.63250   1.376465 0.2486269388231793775467793       0.005080517
# 4 valence_f:action_f   1  35    5.674165  20.72565   9.582127 0.0038530657843402452318593     * 0.022761039

### Incorrect only:
ezANOVA(data = mydata_long_incor, dv = .(score), wid = .(PPN_n),
        within = .(valence_f, action_f), type = 3, 
        detailed = TRUE)
#               Effect DFn DFd           SSn        SSd           F                            p p<.05          ges
# 1        (Intercept)   1  35 1705.06407333 151.432292 394.0853151 0.00000000000000000001247157     * 0.9065963495
# 2          valence_f   1  35    2.04171758   8.273650   8.6370721 0.00579776338896054958121917     * 0.0114891140
# 3           action_f   1  35    0.52390067   6.619329   2.7701484 0.10496630982134541598682631       0.0029734804
# 4 valence_f:action_f   1  35    0.06369474   9.341892   0.2386365 0.62823956572474093373870119       0.0003624562

# ============================================================================================================== #
#### Simple effects with t-tests: ####

### 1) Correct congruent vs. incongruent:  
mydata_wide$CongCor <- (mydata_wide$G2WCor + mydata_wide$NG2ACor)/2
mydata_wide$IncongCor <- (mydata_wide$G2ACor + mydata_wide$NG2WCor)/2
t.test(mydata_wide$CongCor,mydata_wide$IncongCor, paired = T) # t = -3.0955, df = 35, p-value = 0.003853
cohen.d(mydata_wide$CongCor,mydata_wide$IncongCor, paired = T) # d estimate: -0.5159168

### 2) Incorrect Win vs. Avoid:  
mydata_wide$WinIncor <- (mydata_wide$G2WIncor + mydata_wide$NG2WIncor)/2
mydata_wide$AvoidIncor <- (mydata_wide$G2AIncor + mydata_wide$NG2AIncor)/2
t.test(mydata_wide$WinIncor,mydata_wide$AvoidIncor, paired = T) # t = 2.7601, df = 33, p-value = 0.009358
cohen.d(mydata_wide$WinIncor,mydata_wide$AvoidIncor, paired = T) # d estimate: 0.4733579

### 3) Incorrect congruent vs. incongruent:  
mydata_wide$CongIncor <- (mydata_wide$G2WIncor + mydata_wide$NG2AIncor)/2
mydata_wide$IncongIncor <- (mydata_wide$G2AIncor + mydata_wide$NG2WIncor)/2
t.test(mydata_wide$CongIncor,mydata_wide$IncongIncor, paired = T) # t = -0.4885, df = 35, p-value = 0.6282
cohen.d(mydata_wide$CongIncor,mydata_wide$IncongIncor, paired = T) # d estimate: 0.08141738

# END