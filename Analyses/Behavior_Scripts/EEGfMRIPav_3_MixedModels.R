#' EEGfMRIPav_3_MixedModels()
#' 
#' Perform mixed-effects linear and logistic regression reported in behavioral results section.
#' Interactive script, run step-by-step, NOT in one Go.
#' Models reported in behavior are marked with "GOES INTO"...
#' 
#' Requires that .csv files have been recoded and concatenated using 
#' EEGfMRIPav_1_Read_Initial.R
#' 
#' Calls EEGfMRIPav_2_Preprocess_Automated.R
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
library(car) # Version ‘3.0.10’ # for Anova()

library(lme4) # Version ‘1.1.26’ # for lmer and glmer
library(afex) # Version ‘0.28.1’ # for mixed
library(emmeans) # Version ‘1.5.5.1’ # for lsmeans and glht
library(effects) # Version ‘4.2.0’ # for effects plots

# ====================================================================================================================== #
#### Set options, delete existing objects: ####

cat("Penalize scientific notation\n")
options(scipen = 20) # scientific notation: penalize more heavily

cat("Set factors to sum-to-zero coding\n")
options(contrasts = c("contr.sum", "contr.poly")) # default factor coding schemes

# rm(list=ls()) # delete all previously loaded objects

# ====================================================================================================================== #
#### Specify directories: ####

cat("Set directories\n")
rootdir   <- "/project/3017042.02/" # root directory--needs to be adapted to users' folder structure

scriptdir <- paste0(rootdir,"Analyses/Behavior_Scripts/")
inputdir  <- paste0(rootdir,"Log/Behavior/Behavior_concatenated/")
moddir <- paste0(rootdir,"Log/Behavior/Behavior_Models/")

# ====================================================================================================================== #
#### Load full data set ####

## Load data:
cat("Load data\n")
source(paste0(scriptdir,"EEGfMRIPav_2_Preprocess_Automated.R")); # Load function
mydata <- EEGfMRIPav_2_Preprocess_Automated() # Apply function

# ====================================================================================================================== #
#### Initial checks for descriptives/ how many responses got recoded: ####

#### 0a) Check invalid responses: ##### 
sum(is.na(mydata$outcome_n)) # 43
mean(is.na(mydata$outcome_n)) # 0.001866319
table(is.na(mydata$outcome_n), mydata$PPN_f) # max. 14 for subject 30
#         1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27  28
# FALSE 640 640 637 640 640 640 640 634 640 639 639 640 640 640 640 639 640 640 640 640 640 640 640 640 636 637 640 640
# TRUE    0   0   3   0   0   0   0   6   0   1   1   0   0   0   0   1   0   0   0   0   0   0   0   0   4   3   0   0
# 
#        29  30  31  32  33  34  35  36
# FALSE 635 626 637 640 640 640 638 640
# TRUE    5  14   3   0   0   0   2   0

# ---------------------------------------------------------------------------------- #
#### 0b) Check too fast RTs: ##### 

sum(mydata$RT < 0.2, na.rm = T) # 12
mean(mydata$RT < 0.2, na.rm = T) # 0.0008007474
table(mydata$RT < 0.2, mydata$PPN_f) # max 5 for sub 030
#         1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27  28
# FALSE 359 419 392 445 452 332 380 505 237 381 446 484 352 430 416 369 470 453 303 308 509 344 407 446 518 434 345 434
# TRUE    0   0   0   0   0   0   0   0   0   0   1   0   0   0   0   0   0   0   0   0   0   0   0   0   1   2   0   0
# 
#        29  30  31  32  33  34  35  36
# FALSE 332 514 408 451 456 534 505 404
# TRUE    0   5   0   0   0   0   3   0

# ---------------------------------------------------------------------------------- #
#### 0c) Check too late RTs: ##### 

sum(mydata$RT > 1.3035, na.rm = T) # 980
mean(mydata$RT > 1.3035, na.rm = T) # 0.06539437
table(mydata$RT > 1.3035, mydata$PPN_f) # 
#         1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18  19  20  21  22  23  24  25  26  27  28
# FALSE 347 402 369 433 442 321 374 454 230 376 433 442 344 291 390 366 442 441 279 285 480 332 345 428 463 404 330 425
# TRUE   12  17  23  12  10  11   6  51   7   5  14  42   8 139  26   3  28  12  24  23  29  12  62  18  56  32  15   9
# 
#        29  30  31  32  33  34  35  36
# FALSE 306 483 405 365 397 513 471 398
# TRUE   26  36   3  86  59  21  37   6
tmp <- as.matrix(table(mydata$RT > 1.3035, mydata$PPN_f))
max(tmp[2,]) # 139

# ---------------------------------------------------------------------------------- #
#### 0d) Summary statistics per subject: ####

round(tapply(mydata$is_go_cor_n,mydata$PPN_f,mean,na.rm=T),3)
round(tapply(mydata$ACC_n,mydata$PPN_f,mean,na.rm=T),3)
round(tapply(mydata$RT,mydata$PPN_f,mean,na.rm=T),3)

# ---------------------------------------------------------------------------------- #
#### 0e) Compare valid vs. invalid subjects: ####

mydata$exclude_f <- factor(ifelse(mydata$PPN_n %in% c(1,11,15,19,21,25,26),"exclude","include"))
round(tapply(mydata$is_go_cor_n,mydata$exclude_f,mean,na.rm=T),3)
round(tapply(mydata$ACC_n,mydata$exclude_f,mean,na.rm=T),3)
round(tapply(mydata$RT,mydata$exclude_f,mean,na.rm=T),3)

# ====================================================================================================================== #
#### Prepare sub-data sets: ####

## 1) Only Session 1
mydata_ses1 <- subset(mydata, Session_n ==1)

## 2) Only session 2
mydata_ses2 <- subset(mydata, Session_n ==2)

## 3) Select only trials with correct Go responses (actual response):
mydata_cor_go <- subset(mydata, stim_go_f=="Go" & (is_go_cor_n==0 | is_go_cor_n==1 & ACC_cor_n ==1)) # not incorrect Go

## 4) Select only trials with incorrect Go responses (actual response):
mydata_incor_go <- subset(mydata, stim_go_f=="Go" & (is_go_cor_n==0 | is_go_cor_n==1 & ACC_cor_n ==0)) # not correct Go

## 5) Select Go cues (to test differences in accuracy between valence):
mydata_go <- subset(mydata, stim_go_f=="Go") # not incorrect Go

# ====================================================================================================================== #
#### 1a) Proportion Go responses corrected (is_go_cor_n; wrong key presses as Go) as function of valence and required action ####

## Mixed-effects logistic regression:
mod <- glmer(is_go_cor_n ~ stim_go_f*stim_win_f + (stim_go_f*stim_win_f|PPN_f), data = mydata, family = binomial())
saveRDS(mod,paste0(moddir,"mod_pGo_cor_valence_action.rds"))
mod <- readRDS(paste0(moddir,"mod_pGo_cor_valence_action.rds"))
summary(mod)
# Fixed effects:
#                        Estimate Std. Error z value          Pr(>|z|)    
# (Intercept)             0.62123    0.08964   6.930 0.000000000004207 ***
# stim_go_f1              0.81505    0.11343   7.185 0.000000000000671 ***
# stim_win_f1            -0.42304    0.07327  -5.774 0.000000007756351 ***
# stim_go_f1:stim_win_f1 -0.03025    0.06825  -0.443             0.658 

## P-values with chisquare tests:
Anova(mod,type="III")
#                        Chisq Df         Pr(>Chisq)    
# (Intercept)          48.0255  1 0.0000000000042072 ***
# stim_go_f            51.6268  1 0.0000000000006712 ***
# stim_win_f           33.3353  1 0.0000000077563506 ***
# stim_go_f:stim_win_f  0.1965  1             0.6576

## P-values with LRTs:
mod_LRT <- mixed(mod, mydata, family = binomial(), method = "LRT", type = "III",
                 control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
saveRDS(mod_LRT,paste0(moddir,"LRT_mod_pGo_cor_valence_action.rds"))
mod_LRT <- readRDS(paste0(moddir,"LRT_mod_pGo_cor_valence_action.rds"))
anova(mod_LRT)
#                      Df   Chisq Chi Df    Pr(>Chisq)    
# stim_go_f            13 32.0075      1 0.00000001536 ***
# stim_win_f           13 23.6951      1 0.00000112871 ***
# stim_go_f:stim_win_f 13  0.1957      1        0.6583    
# --> GOES INTO PAPER 

## Plots:
plot(effect("stim_go_f",mod))
plot(effect("stim_win_f",mod))
plot(effect("stim_go_f:stim_win_f",mod))

# ====================================================================================================================== #
#### Pairwise bonferroni-corrected t-tests: ####

lsmeans(mod, pairwise ~ stim_go_f * stim_win_f, adj = "bonferroni")
# stim_go_f stim_win_f lsmean    SE  df asymp.LCL asymp.UCL
# Go        Avoid       0.983 0.117 Inf     0.754     1.212
# NoGo      Avoid      -0.587 0.156 Inf    -0.892    -0.281
# Go        Win         1.890 0.198 Inf     1.501     2.278
# NoGo      Win         0.199 0.215 Inf    -0.223     0.621

# contrast                estimate        SE df z.ratio p.value
# Go,Avoid - NoGo,Avoid    1.570 0.220 Inf  7.136  <.0001 
# Go,Avoid - Go,Win       -0.907 0.163 Inf -5.564  <.0001 
# Go,Avoid - NoGo,Win      0.784 0.264 Inf  2.964  0.0182 
# NoGo,Avoid - Go,Win     -2.476 0.276 Inf -8.986  <.0001 
# NoGo,Avoid - NoGo,Win   -0.786 0.232 Inf -3.391  0.0042 
# Go,Win - NoGo,Win        1.691 0.303 Inf  5.579  <.0001 

# ====================================================================================================================== #
#### On subset of people used in EEG-fMRI analyses (N = 30) ####

## Select subsect of data;
subsetdata <- subset(mydata, !(PPN_n %in% c(1,11,15,19,21,25))) # cue-locked (N = 30)
length(unique(subsetdata$PPN_n)) # 29
sort(unique(subsetdata$PPN_n))

## Mixed-effects logistic regression:
mod <- glmer(is_go_cor_n ~ stim_go_f*stim_win_f + (stim_go_f*stim_win_f|PPN_f), data = subsetdata, family = binomial())
saveRDS(mod,paste0(moddir,"mod_pGo_cor_valence_action_withoutInvalid.rds"))
mod <- readRDS(paste0(moddir,"mod_pGo_cor_valence_action_withoutInvalid.rds"))
summary(mod)
#                         Estimate Std. Error z value       Pr(>|z|)    
# (Intercept)             0.611881   0.096141   6.364 0.000000000196 ***
# stim_go_f1              0.836242   0.130959   6.386 0.000000000171 ***
# stim_win_f1            -0.432223   0.084599  -5.109 0.000000323708 ***
# stim_go_f1:stim_win_f1 -0.007576   0.066813  -0.113           0.91 

## P-values with LRTs:
mod_LRT <- mixed(mod, subsetdata, family = binomial(), method = "LRT", type = "III",
                 control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
saveRDS(mod_LRT,paste0(moddir,"LRT_mod_pGo_cor_valence_action_withoutInvalid.rds"))
mod_LRT <- readRDS(paste0(moddir,"LRT_mod_pGo_cor_valence_action_withoutInvalid.rds"))
anova(mod_LRT)
#                      Df   Chisq Chi Df   Pr(>Chisq)    
# stim_go_f            13 25.7728      1 0.0000003841 ***
# stim_win_f           13 18.8726      1 0.0000139743 ***
# stim_go_f:stim_win_f 13  0.0129      1       0.9097   
# --> GOES INTO SUPPLEMENTARY MATERIAL

# ====================================================================================================================== #
#### Combine required action and cue valence to one factor: ####

## Combine cue valence and required action into one single factor:
mydata$stim_go_win_f <- factor(paste(mydata$stim_go_o,mydata$stim_win_o,sep = "2"))
contrasts(mydata$stim_go_win_f)

## Mixed-effects logistic regression:
mod <- glmer(is_go_cor_n ~ stim_go_win_f + (stim_go_win_f|PPN_f), data = mydata, family = binomial())

## P-values with chi-square tests:
Anova(mod,type = "III")
#                Chisq Df            Pr(>Chisq)    
# (Intercept)   81.104  1 < 0.00000000000000022 ***
# stim_go_win_f 83.390  3 < 0.00000000000000022 ***

## Plot:
plot(effect("stim_go_win_f",mod))

## Follow-ups with lsmeans:
lsmeans(mod, pairwise ~ stim_go_win_f)
lsmeans(mod, pairwise ~ stim_go_win_f, adj = "bonferroni")
summary(glht(mod, linfct = mcp(stim_go_win_f = "Tukey")))
summary(effect("stim_go_win_f", mod, se=TRUE))

# ====================================================================================================================== #
#### Only session 1: ####

## Mixed-effects logistic regression:
mod <- glmer(is_go_cor_n ~ stim_go_f*stim_win_f + (stim_go_f*stim_win_f|PPN_f), data = mydata_ses1, family = binomial(),
             control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))

## P-values with chisquare tests:
Anova(mod,type="III")
#                        Chisq Df      Pr(>Chisq)    
# (Intercept)          38.1634  1 0.0000000006506 ***
# stim_go_f            34.0179  1 0.0000000054607 ***
# stim_win_f           38.2073  1 0.0000000006361 ***
# stim_go_f:stim_win_f  0.0218  1          0.8825    

## P-values with LRTs:
mod_LRT <- mixed(mod, mydata_ses1, family = binomial(), method = "LRT", type = "III",
                 control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#                      Df   Chisq Chi Df   Pr(>Chisq)    
# stim_go_f            13 24.0957      1 0.0000009166 ***
# stim_win_f           13 26.3662      1 0.0000002824 ***
# stim_go_f:stim_win_f 13  0.0219      1       0.8823 

# ====================================================================================================================== #
#### Only session 2: ####

## Mixed-effects logistic regression:
mod <- glmer(is_go_cor_n ~ stim_go_f*stim_win_f + (stim_go_f*stim_win_f|PPN_f), data = mydata_ses2, family = binomial(),
             control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))

## P-values with chisquare tests:
Anova(mod,type="III")
#                        Chisq Df          Pr(>Chisq)    
# (Intercept)          21.5783  1 0.00000339669069406 ***
# stim_go_f            55.4024  1 0.00000000000009822 ***
# stim_win_f           10.7786  1            0.001027 ** 
# stim_go_f:stim_win_f  1.2228  1            0.268816 

## P-values with LRTs:
mod_LRT <- mixed(mod, mydata_ses2, family = binomial(), method = "LRT", type = "III",
                 control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#                      Df   Chisq Chi Df     Pr(>Chisq)    
# stim_go_f            13 33.6248      1 0.000000006683 ***
# stim_win_f           13  9.4224      1       0.002144 ** 
# stim_go_f:stim_win_f 13  1.1991      1       0.273509 

## Plots:
plot(effect("stim_go_f",mod))
plot(effect("stim_win_f",mod))

# ====================================================================================================================== #
#### Only correct Go vs. NoGo for Go cues: ####

## Mixed-effects logistic regression:
mod <- glmer(is_go_cor_n ~ stim_win_f + (stim_win_f|PPN_f), data = mydata_cor_go, family = binomial(),
             control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
summary(mod)
#             Estimate Std. Error z value Pr(>|z|)  
# (Intercept)   0.2010     0.1796   1.119   0.2630  
# stim_win_f1  -0.2323     0.1188  -1.956   0.0505 

## P-values with chisquare tests:
Anova(mod,type="III")
#              Chisq Df  Pr(>Chisq)    
# (Intercept) 1.2527  1    0.26303  
# stim_win_f  3.8252  1    0.05049 

## P-values with LRTs:
mod_LRT <- mixed(mod, mydata_cor_go, family = binomial(), method = "LRT", type = "III",
                 control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#            Df  Chisq Chi Df Pr(>Chisq)    
# stim_win_f  4 3.5197      1    0.06064 .

## Plot:
plot(effect("stim_win_f",mod))

# ====================================================================================================================== #
#### Only incorrect Go vs. NoGo for Go cues: ####

## Mixed-effects logistic regression:
mod <- glmer(is_go_cor_n ~ stim_win_f + (stim_win_f|PPN_f), data = mydata_incor_go, family = binomial(),
             control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
summary(mod)
#             Estimate Std. Error z value           Pr(>|z|)    
# (Intercept)  1.01273    0.13187   7.680 0.0000000000000159 ***
# stim_win_f1 -0.50397    0.08538  -5.902 0.0000000035806244 ***
  
## P-values with chisquare tests:
Anova(mod,type="III")
#               Chisq Df          Pr(>Chisq)    
# (Intercept) 58.981  1 0.00000000000001592 ***
# stim_win_f  34.839  1 0.00000000358062443 ***

## P-values with LRTs:
mod_LRT <- mixed(mod, mydata_incor_go, family = binomial(), method = "LRT", type = "III",
                 control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#            Df  Chisq Chi Df Pr(>Chisq)    
# stim_win_f  4 24.786      1 0.0000006406 ***

# ====================================================================================================================== #
#### 1b) Proportion Go responses uncorrected (is_go_n) as function of valence and required action ####

# is_go_n is raw Go/NoGo variable read from raw data, not changed
table(mydata$is_go_n==mydata$is_go_cor_n) # identical with corrected Go variable

## Mixed effecs logistic regression:
mod <- glmer(is_go_n ~ stim_go_f*stim_win_f + (stim_go_f*stim_win_f|PPN_f), data = mydata, family = binomial())
saveRDS(mod,paste0(moddir,"mod_pGo_valence_action.rds"))
mod <- readRDS(paste0(moddir,"mod_pGo_valence_action.rds"))
summary(mod)
# Fixed effects:
#                        Estimate Std. Error z value          Pr(>|z|)    
# (Intercept)             0.62123    0.08964   6.930 0.000000000004207 ***
# stim_go_f1              0.81505    0.11343   7.185 0.000000000000671 ***
# stim_win_f1            -0.42304    0.07327  -5.774 0.000000007756351 ***
# stim_go_f1:stim_win_f1 -0.03025    0.06825  -0.443             0.658  

## P-values with chisquare tests:
Anova(mod,type="III")
#                         Chisq Df         Pr(>Chisq)    
# (Intercept)          48.0255  1 0.0000000000042072 ***
# stim_go_f            51.6268  1 0.0000000000006712 ***
# stim_win_f           33.3353  1 0.0000000077563506 ***
# stim_go_f:stim_win_f  0.1965  1             0.6576 

## P-values with LRTs:
mod_LRT <- mixed(mod, mydata, family = binomial(), method = "LRT", type = "III",
                 control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
saveRDS(mod_LRT,paste0(moddir,"LRT_mod_pGo_valence_action.rds"))
mod_LRT <- readRDS(paste0(moddir,"LRT_mod_pGo_valence_action.rds"))
anova(mod_LRT)
#                      Df   Chisq Chi Df    Pr(>Chisq)    
# stim_go_f            13 32.0075      1 0.00000001536 ***
# stim_win_f           13 23.6951      1 0.00000112871 ***
# stim_go_f:stim_win_f 13  0.1957      1        0.6583

## Plots:
plot(effect("stim_go_f",mod))
plot(effect("stim_win_f",mod))
plot(effect("stim_go_f:stim_win_f",mod))

# ====================================================================================================================== #
#### Robustness checks: ####

# ======================================================= #
### Only session 1:

## Mixed-effects logistic regression:
mod <- glmer(is_go_n ~ stim_go_f*stim_win_f + (stim_go_f*stim_win_f|PPN_f), data = mydata_ses1, family = binomial())

## P-values with chisquare tests:
Anova(mod,type="III")
#                        Chisq Df      Pr(>Chisq)    
# (Intercept)          38.1730  1 0.0000000006474 ***
# stim_go_f            34.0456  1 0.0000000053836 ***
# stim_win_f           38.2142  1 0.0000000006339 ***
# stim_go_f:stim_win_f  0.0218  1          0.8825    

## P-values with LRTs:
mod_LRT <- mixed(mod, mydata_ses1, family = binomial(), method = "LRT", type = "III", control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#                      Df   Chisq Chi Df   Pr(>Chisq)    
# stim_go_f            13 24.0957      1 0.0000009166 ***
# stim_win_f           13 26.3662      1 0.0000002824 ***
# stim_go_f:stim_win_f 13  0.0219      1       0.8823

# ======================================================= #
### Only session 2:

# Mixed-effects logistic regression:
mod <- glmer(is_go_n ~ stim_go_f*stim_win_f + (stim_go_f*stim_win_f|PPN_f), data = mydata_ses2, family = binomial(),
             control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))

## P-values with chisquare tests:
Anova(mod,type="III")
#                        Chisq Df          Pr(>Chisq)    
# (Intercept)          21.5783  1 0.00000339669069406 ***
# stim_go_f            55.4024  1 0.00000000000009822 ***
# stim_win_f           10.7786  1            0.001027 ** 
# stim_go_f:stim_win_f  1.2228  1            0.268816

## P-values with LRTs:
mod_LRT <- mixed(mod, mydata_ses2, family = binomial(), method = "LRT", type = "III", control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#                      Df   Chisq Chi Df     Pr(>Chisq)    
# stim_go_f            13 33.6248      1 0.000000006683 ***
# stim_win_f           13  9.4224      1       0.002144 ** 
# stim_go_f:stim_win_f 13  1.1991      1       0.273509 

# ======================================================= #
### Only correct Go:

### Mixed-effects logistic regression:
mod <- glmer(is_go_n ~ stim_win_f + (stim_win_f||PPN_f), data = mydata_cor_go, family = binomial(),
             control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))

## P-values with chisquare tests:
Anova(mod,type="III")
#             Chisq Df Pr(>Chisq)  
# (Intercept) 1.253  1    0.26298  
# stim_win_f  3.826  1    0.05046 .

## P-values with LRTs:
mod_LRT <- mixed(mod, mydata_cor_go, family = binomial(), method = "LRT", type = "III", control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#            Df  Chisq Chi Df Pr(>Chisq)    
# stim_win_f  5 3.5197      1    0.06064 .

# ======================================================= #
### Only incorrect Go:

## Mixed-effects logistic regression:
mod <- glmer(is_go_n ~ stim_win_f + (stim_win_f|PPN_f), data = mydata_incor_go, family = binomial())

## P-values with chisquare tests:
Anova(mod,type="III")
#              Chisq Df          Pr(>Chisq)    
# (Intercept) 58.892  1 0.00000000000001666 ***
# stim_win_f  34.813  1 0.00000000362863708 ***

## P-values with LRTs:
mod_LRT <- mixed(mod, mydata_incor_go, family = binomial(), method = "LRT", type = "III", control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#            Df  Chisq Chi Df   Pr(>Chisq)    
# stim_win_f  4 24.786      1 0.0000006406 ***

# ============================================================================================================================= #
# ============================================================================================================================= #
# ============================================================================================================================= #
# ============================================================================================================================= #
#### 2a) Accuracy corrected (ACC_cor_n; wrong keys as right side) as function of valence and required action ####

## Mixed-effects logistic regression:
mod <- glmer(ACC_cor_n ~ stim_go_f*stim_win_f + (stim_go_f*stim_win_f|PPN_f), data = mydata, family = binomial())
saveRDS(mod,paste0(moddir,"mod_ACC_cor_valence_action.rds"))
mod <- readRDS(paste0(moddir,"mod_ACC_cor_valence_action.rds"))
summary(mod)
# Fixed effects:
#                        Estimate Std. Error z value             Pr(>|z|)    
# (Intercept)            -0.89492    0.05802 -15.425 < 0.0000000000000002 ***
# stim_go_f1              0.32026    0.07159   4.474           0.00000769 ***
# stim_win_f1            -0.12745    0.06016  -2.119               0.0341 *  
# stim_go_f1:stim_win_f1  0.35416    0.07468   4.742           0.00000211 ***

## P-values with chisquare tests:
Anova(mod, type = "III")
#                         Chisq Df            Pr(>Chisq)    
# (Intercept)          237.9384  1 < 0.00000000000000022 ***
# stim_go_f             20.0142  1           0.000007687 ***
# stim_win_f             4.4885  1               0.03412 *  
# stim_go_f:stim_win_f  22.4885  1           0.000002114 ***

## P-values with LRTs:
mod_LRT <- mixed(mod, mydata, family = binomial(), method = "LRT", type = "III",
                 control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
saveRDS(mod_LRT,paste0(moddir,"LRT_mod_ACC_cor_valence_action.rds"))
mod_LRT <- readRDS(paste0(moddir,"LRT_mod_ACC_cor_valence_action.rds"))
anova(mod_LRT)
#                      Df   Chisq Chi Df Pr(>Chisq)    
# stim_go_f            13 15.8626      1 0.00006811 ***
# stim_win_f           13  4.2217      1    0.03991 *  
# stim_go_f:stim_win_f 13 17.4819      1 0.00002901 ***

## Plots:
plot(effect("stim_win_f",mod))
plot(effect("stim_go_f",mod))
plot(effect("stim_go_f:stim_win_f",mod))

# ============================================================================================================================= #
#### Followups: ####

lsmeans(mod, pairwise ~ stim_go_f * stim_win_f, adj = "bonferroni")
# $lsmeans
# stim_go_f stim_win_f lsmean     SE  df asymp.LCL asymp.UCL
# Go        Avoid      -0.348 0.0456 Inf    -0.437    -0.259
# NoGo      Avoid      -1.697 0.1333 Inf    -1.958    -1.436
# Go        Win        -0.801 0.1317 Inf    -1.059    -0.543
# NoGo      Win        -0.734 0.1832 Inf    -1.093    -0.375

# $ contrasts
# contrast                 estimate        SE df z.ratio p.value
# Go,Avoid - NoGo,Avoid -0.18653000 0.1580305 NA  -1.180  1.0000
# Go,Avoid - Go,Win     -0.20182566 0.1628117 NA  -1.240  1.0000
# Go,Avoid - NoGo,Win    0.60568315 0.2457346 NA   2.465  0.0823
# NoGo,Avoid - Go,Win   -0.01529567 0.2213846 NA  -0.069  1.0000
# NoGo,Avoid - NoGo,Win  0.79221315 0.2517314 NA   3.147  0.0099
# Go,Win - NoGo,Win      0.80750881 0.2896179 NA   2.788  0.0318

# ============================================================================================================================= #
#### Combine factors to one factor ####

## Combine cue valence and required action into one factor:
mydata$stim_go_win_f <- factor(paste(mydata$stim_go_f,mydata$stim_win_f,sep = "2"))

## Mixed-effects logistic regression:
mod <- glmer(ACC_cor_n ~ stim_go_win_f + (stim_go_win_f|PPN_f), data = mydata, family = binomial())

## P-values with chisquare tests:
Anova(mod,type = "III")
#                 Chisq Df            Pr(>Chisq)    
# (Intercept)   237.880  1 < 0.00000000000000022 ***
# stim_go_win_f  81.066  3 < 0.00000000000000022 ***

## Plot:
plot(effect("stim_go_win_f",mod))

### Follow-ups:
lsmeans(mod, pairwise ~ stim_go_win_f)
lsmeans(mod, pairwise ~ stim_go_win_f, adj = "bonferroni")
summary(glht(mod, linfct = mcp(stim_go_win_f = "Tukey")))
summary(effect("stim_go_win_f", mod, se=TRUE))

# ============================================================================================================================= #
#### Only session 1: ####

## Mixed-effects logistic regression:
mod <- glmer(ACC_cor_n ~ stim_go_f*stim_win_f + (stim_go_f*stim_win_f|PPN_f), data = mydata_ses1, family = binomial())

## P-values with chisquare tests:
Anova(mod, type = "III")
#                         Chisq Df            Pr(>Chisq)    
# (Intercept)          110.2896  1 < 0.00000000000000022 ***
# stim_go_f              2.2086  1                0.1372    
# stim_win_f             0.4567  1                0.4992    
# stim_go_f:stim_win_f  22.9865  1           0.000001631 ***

## P-values with LRTs:
mod_LRT <- mixed(mod, mydata_ses1, family = binomial(), method = "LRT", type = "III",
                 control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#                      Df   Chisq Chi Df Pr(>Chisq)    
# stim_go_f            13  2.1244      1     0.1450    
# stim_win_f           13  0.4498      1     0.5024    
# stim_go_f:stim_win_f 13 17.8837      1 0.00002348 ***

# ============================================================================================================================= #
#### Only session 2; ####

## Mixed-effects logistic regression:
mod <- glmer(ACC_cor_n ~ stim_go_f*stim_win_f + (stim_go_f*stim_win_f|PPN_f), data = mydata_ses2, family = binomial())

## P-values with chisquare tests:
Anova(mod, type = "III")
#                        Chisq Df            Pr(>Chisq)    
# (Intercept)          96.0718  1 < 0.00000000000000022 ***
# stim_go_f            23.7557  1           0.000001094 ***
# stim_win_f            0.3309  1               0.56515    
# stim_go_f:stim_win_f  6.2434  1               0.01247 * 

## P-values with LRTs:
mod_LRT <- mixed(mod, mydata_ses2, family = binomial(), method = "LRT", type = "III",
                 control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#                      Df   Chisq Chi Df Pr(>Chisq)    
# stim_go_f            13 17.8544      1 0.00002385 ***
# stim_win_f           13  0.3261      1    0.56797    
# stim_go_f:stim_win_f 13  5.7245      1    0.01673 *  

## Plot:
plot(effect("stim_go_f:stim_win_f",mod))

# ============================================================================================================================= #
#### Only correct Go responses to Go cues: ####

## Mixed-effects logistic regression:
mod <- glmer(ACC_cor_n ~ stim_win_f + (stim_win_f|PPN_f), data = mydata_cor_go, family = binomial(),
             control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))

## P-values with chisquare tests:
Anova(mod, type = "III")
#               Chisq Df           Pr(>Chisq)    
# (Intercept) 94.5731  1 < 0.0000000000000002 ***
# stim_win_f   4.4767  1              0.03436 *  

## P-values with LRTs:
mod_LRT <- mixed(mod, mydata_cor_go, family = binomial(), method = "LRT", type = "III",
                 control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#            Df  Chisq Chi Df Pr(>Chisq)  
# stim_win_f  4 4.3488      1    0.03704 *

# ============================================================================================================================= #
## Only incorrect Go: response is always wrong, so constant
# ============================================================================================================================= #
#### Accuracy as function of valence for Go cues: ####

## Mixed-effects logistic regression:
mod <- glmer(ACC_cor_n ~ stim_win_f + (stim_win_f|PPN_f), data = mydata_go, family = binomial())
summary(mod)
#             Estimate Std. Error z value             Pr(>|z|)    
# (Intercept) -0.57400    0.06746  -8.508 < 0.0000000000000002 ***
# stim_win_f1  0.22598    0.07124   3.172              0.00151 ** 
  
## P-values with chisquare tests:
Anova(mod, type = "III")
#               Chisq Df            Pr(>Chisq)    
# (Intercept) 74.2404  1 < 0.00000000000000022 ***
# stim_win_f   8.3931  1              0.003766 ** 

## P-values with LRTs:
mod_LRT <- mixed(mod, mydata_go, family = binomial(), method = "LRT", type = "III",
                 control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#            Df  Chisq Chi Df Pr(>Chisq)
# stim_win_f  4 7.6921      1   0.005546 **

# ============================================================================================================================= #
#### 2b) Accuracy uncorrected (ACC_n) as function of valence and required action ####

table(mydata$ACC_n==mydata$ACC_cor_n)
# FALSE  TRUE 
# 13859  9181 

## Mixed-effects logistic regression:
mod <- glmer(ACC_n ~ stim_go_f*stim_win_f + (stim_go_f*stim_win_f|PPN_f), data = mydata, family = binomial())
saveRDS(mod,paste0(moddir,"mod_ACC_valence_action.rds"))
mod <- readRDS(paste0(moddir,"mod_ACC_valence_action.rds"))
summary(mod)
# Fixed effects:
#                        Estimate Std. Error z value Pr(>|z|)    
# (Intercept)             0.17737    0.11380   1.559 0.119071    
# stim_go_f1             -0.01662    0.08319  -0.200 0.841615    
# stim_win_f1             0.12967    0.07029   1.845 0.065086 .  
# stim_go_f1:stim_win_f1 -0.26357    0.07020  -3.754 0.000174 ***

## P-values with chisquare tests:
Anova(mod,type="III")
#                       Chisq Df Pr(>Chisq)    
# (Intercept)           2.4295  1  0.1190710    
# stim_go_f             0.0399  1  0.8416151    
# stim_win_f            3.4028  1  0.0650856 .  
# stim_go_f:stim_win_f 14.0952  1  0.0001738 ***

## P-values with LRTs:
mod_LRT <- mixed(mod, mydata, family = binomial(), method = "LRT", type = "III",
                 control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
saveRDS(mod_LRT,paste0(moddir,"LRT_mod_ACC_valence_action.rds"))
mod_LRT <- readRDS(paste0(moddir,"LRT_mod_ACC_valence_action.rds"))
anova(mod_LRT)
#                      Df   Chisq Chi Df Pr(>Chisq)    
# stim_go_f            13  0.0399      1  0.8417335    
# stim_win_f           13  3.2623      1  0.0708913 .  
# stim_go_f:stim_win_f 13 11.9272      1  0.0005532 ***

## Plots:
plot(effect("stim_go_f",mod))
plot(effect("stim_win_f",mod))
plot(effect("stim_go_f:stim_win_f",mod))
#- -> Higher accuracy for congruent trials (win-go and avoid-NoGo), lower for incongruent ones

# ============================================================================================================================= #
#### Robustness checks: ####

# ======================================================= #
### Only session 1:

## Mixed-effects logistic regression:
mod <- glmer(ACC_n ~ stim_go_f*stim_win_f + (stim_go_f*stim_win_f|PPN_f), data = mydata_ses1, family = binomial())
Anova(mod,type = "III")
#                        Chisq Df Pr(>Chisq)    
# (Intercept)           0.1406  1  0.7076672    
# stim_go_f             0.0182  1  0.8928313    
# stim_win_f            4.8088  1  0.0283151 *  
# stim_go_f:stim_win_f 13.8969  1  0.0001931 ***

## P-values with LRTs:
mod_LRT <- mixed(mod, mydata_ses1, family = binomial(), method = "LRT", type = "III",
                 control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#                      Df   Chisq Chi Df Pr(>Chisq)    
# stim_go_f            13  0.0182      1  0.8928069    
# stim_win_f           13  4.5350      1  0.0332086 *  
# stim_go_f:stim_win_f 13 11.7732      1  0.0006009 ***

# ======================================================= #
### Only session 2:

## Mixed-effects logistic regression:
mod <- glmer(ACC_n ~ stim_go_f*stim_win_f + (stim_go_f*stim_win_f|PPN_f), data = mydata_ses2, family = binomial())

## P-values with chisquare tests:
Anova(mod,type="III")
#                       Chisq Df Pr(>Chisq)   
# (Intercept)          8.6208  1   0.003323 **
# stim_go_f            0.0630  1   0.801835   
# stim_win_f           0.0693  1   0.792301   
# stim_go_f:stim_win_f 5.7832  1   0.016180 * 

## P-values with LRTs:
mod_LRT <- mixed(mod, mydata_ses2, family = binomial(), method = "LRT", type = "III",
                 control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#                      Df  Chisq Chi Df Pr(>Chisq)  
# stim_go_f            13 0.0627      1    0.80227  
# stim_win_f           13 0.0688      1    0.79311  
# stim_go_f:stim_win_f 13 5.3286      1    0.02098 *

# ======================================================= #
### Only correct Go:

## Mixed-effects logistic regression:
mod <- glmer(ACC_n ~ stim_win_f + (stim_win_f|PPN_f), data = mydata_cor_go, family = binomial(),
             control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))

## P-values with chisquare tests:
Anova(mod, type = "III")
#              Chisq Df Pr(>Chisq)  
# (Intercept) 1.2331  1    0.26680  
# stim_win_f  3.7846  1    0.05173 .

## P-values with LRTs:
mod_LRT <- mixed(mod, mydata_cor_go, family = binomial(), method = "LRT", type = "III",
                 control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#            Df  Chisq Chi Df Pr(>Chisq)    
# stim_win_f  4 3.4842      1    0.06196 .

# ======================================================= #
## Only incorrect Go: response always incorrect, thus constant

# ====================================================================================================================== #
### 2c) Accuracy corrected (wrong keys as right side) as a function of congruency ####

## Mixed-effects logistic regression:
mod <- glmer(ACC_cor_n ~ stim_cong_f + (stim_cong_f|PPN_f), data = mydata, family = binomial())
saveRDS(mod,paste0(moddir,"mod_ACC_cor_cong.rds"))
mod <- readRDS(paste0(moddir,"mod_ACC_cor_cong.rds"))
summary(mod)
# Fixed effects:
#              Estimate Std. Error z value             Pr(>|z|)    
# (Intercept)  -0.81862    0.04922 -16.632 < 0.0000000000000002 ***
# stim_cong_f1 -0.32002    0.06828  -4.687           0.00000278 ***

## P-values with chisquare tests:
Anova(mod, type = "III")
#               Chisq Df            Pr(>Chisq)    
# (Intercept) 276.611  1 < 0.00000000000000022 ***
# stim_cong_f  21.965  1           0.000002777 ***

## P-values with LRTs:
mod_LRT <- mixed(mod, mydata, family = binomial(), method = "LRT", type = "III",
                 control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#             Df Chisq Chi Df Pr(>Chisq)    
# stim_cong_f  4 17.16      1 0.00003437 ***

## Plot:
plot(effect("stim_cong_f",mod))

# ====================================================================================================================== #
#### Robustness checks: ####

# ======================================================= #
### Only session 1:

## Mixed-effects logistic regression:
mod <- glmer(ACC_cor_n ~ stim_cong_f + (stim_cong_f|PPN_f), data = mydata_ses1, family = binomial())

## P-values with chisquare tests:
Anova(mod, type = "III")
#               Chisq Df            Pr(>Chisq)    
# (Intercept) 120.509  1 < 0.00000000000000022 ***
# stim_cong_f  25.475  1          0.0000004482 ***

## P-values with LRTs:
mod_LRT <- mixed(mod, mydata_ses1, family = binomial(), method = "LRT", type = "III",
                 control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#             Df  Chisq Chi Df Pr(>Chisq)   
# stim_cong_f  4 19.363      1 0.00001081 ***

# ======================================================= #
### Only session 2:

## Mixed-effects logistic regression:
mod <- glmer(ACC_cor_n ~ stim_cong_f + (stim_cong_f|PPN_f), data = mydata_ses2, family = binomial())

## P-values with chisquare tests:
Anova(mod, type = "III")
#                Chisq Df            Pr(>Chisq)    
# (Intercept) 124.4808  1 < 0.00000000000000022 ***
# stim_cong_f   6.6485  1              0.009924 ** 

## P-values with LRTs:
mod_LRT <- mixed(mod, mydata_ses2, family = binomial(), method = "LRT", type = "III",
                 control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#             Df  Chisq Chi Df Pr(>Chisq)  
# stim_cong_f  4 6.1166      1    0.01339 *

# ======================================================= #
### Only correct Go:

## Mixed-effects logistic regression:
mod <- glmer(ACC_cor_n ~ stim_cong_f + (stim_cong_f|PPN_f), data = mydata_cor_go, family = binomial())

## P-values with chisquare tests:
Anova(mod, type = "III")
#               Chisq Df           Pr(>Chisq)    
# (Intercept) 94.5707  1 < 0.0000000000000002 ***
# stim_cong_f  4.4762  1              0.03437 * 

## P-values with LRTs:
mod_LRT <- mixed(mod, mydata_cor_go, family = binomial(), method = "LRT", type = "III",
                 control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#             Df  Chisq Chi Df Pr(>Chisq)  
# stim_cong_f  4 4.3488      1    0.03704 *

# ======================================================= #
## Only incorrect Go: response is always wrong, so constant

# ====================================================================================================================== #
### 2d) Accuracy as a function of congruency ####

## Mixed-effects logistic regression:
mod <- glmer(ACC_n ~ stim_cong_f + (stim_cong_f|PPN_f), data = mydata, family = binomial())
saveRDS(mod,paste0(moddir,"mod_ACC_cong.rds"))
mod <- readRDS(paste0(moddir,"mod_ACC_cong.rds"))
summary(mod)
# Fixed effects:
#              Estimate Std. Error z value Pr(>|z|)    
# (Intercept)   0.18806    0.10569   1.779 0.075175 .  
# stim_cong_f1  0.23287    0.06081   3.830 0.000128 ***

## P-values with chisquare tests:
Anova(mod, type = "III")
#               Chisq Df Pr(>Chisq)    
# (Intercept)  3.1663  1  0.0751745 .  
# stim_cong_f 14.6659  1  0.0001283 ***
mod_LRT <- mixed(mod, mydata, family = binomial(), method = "LRT", type = "III",
                 control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#             Df Chisq Chi Df Pr(>Chisq)    
# stim_cong_f  4 12.31      1  0.0004506 ***

## Plot:
plot(effect("stim_cong_f",mod))

# ====================================================================================================================== #
#### Robustness checks: ####

# ======================================================= #
### Only session 1:

## Mixed-effects logistic regression:
mod <- glmer(ACC_n ~ stim_cong_f + (stim_cong_f|PPN_f), data = mydata_ses1, family = binomial())

## P-values with chisquare tests:
Anova(mod, type = "III")
#              Chisq Df Pr(>Chisq)    
# (Intercept)  0.029  1  0.8648152    
# stim_cong_f 13.207  1  0.0002789 ***

## P-values with LRTs:
mod_LRT <- mixed(mod, mydata_ses1, family = binomial(), method = "LRT", type = "III",
                 control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#          Df  Chisq Chi Df Pr(>Chisq)    
# stim_cong_f  4 11.257      1   0.000793 ***

# ======================================================= #
### Only session 2:

## Mixed-effects logistic regression:
mod <- glmer(ACC_n ~ stim_cong_f + (stim_cong_f|PPN_f), data = mydata_ses2, family = binomial())

## P-values with chisquare tests:
Anova(mod, type = "III")
#              Chisq Df Pr(>Chisq)   
# (Intercept) 9.2549  1   0.002349 **
# stim_cong_f 5.9760  1   0.014502 * 

## P-values with LRTs:
mod_LRT <- mixed(mod, mydata_ses2, family = binomial(), method = "LRT", type = "III",
                 control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#             Df  Chisq Chi Df Pr(>Chisq)  
# stim_cong_f  4 5.5398      1    0.01859 *

# ======================================================= #
### Only correct Go:

## Mixed-effects logistic regression:
mod <- glmer(ACC_n ~ stim_cong_f + (stim_cong_f|PPN_f), data = mydata_cor_go, family = binomial())

## P-values with chisquare tests:
Anova(mod, type = "III")
#              Chisq Df Pr(>Chisq)  
# (Intercept) 1.2331  1    0.26680  
# stim_cong_f 3.7850  1    0.05171 .

## P-values with LRTs:
mod_LRT <- mixed(mod, mydata_cor_go, family = binomial(), method = "LRT", type = "III",
                 control = glmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#             Df  Chisq Chi Df Pr(>Chisq)    
# stim_cong_f  4 3.4842      1    0.06196 
## Only incorrect Go: response is always wrong, so constant

# ==================================================================================================================== #
# ==================================================================================================================== #
# ==================================================================================================================== #
#### 3a) RTs corrected (delete too fast and too slow responses) as function of valence and required action ####

## Mixed-effects linear regression:
mod <- lmer(RT_nna ~ stim_go_f*stim_win_f + (stim_go_f*stim_win_f|PPN_f), data = mydata,na.action = na.exclude)
saveRDS(mod,paste0(moddir,"mod_RT_nna_valence_action.rds"))
mod <- readRDS(paste0(moddir,"mod_RT_nna_valence_action.rds"))
summary(mod)
# Fixed effects:
#                         Estimate Std. Error        df t value             Pr(>|t|)    
# (Intercept)             0.779358   0.012062 34.931353  64.614 < 0.0000000000000002 ***
# stim_go_f1             -0.027059   0.004513 32.152463  -5.995       0.000001080170 ***
# stim_win_f1             0.053831   0.005801 35.309023   9.280       0.000000000053 ***
# stim_go_f1:stim_win_f1 -0.004512   0.003676 33.541429  -1.227                0.228 

## P-values with chisquare tests:
Anova(mod, type = "III")
#                          Chisq Df            Pr(>Chisq)    
# (Intercept)          4175.0286  1 < 0.00000000000000022 ***
# stim_go_f              35.9445  1         0.00000000203 ***
# stim_win_f             86.1188  1 < 0.00000000000000022 ***
# stim_go_f:stim_win_f    1.5063  1                0.2197 

## P-values with LRTs:
mod_LRT <- mixed(mod, mydata, method = "LRT", type = "III", control = lmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
saveRDS(mod_LRT,paste0(moddir,"LRT_mod_RT_nna_valence_action.rds"))
mod_LRT <- readRDS(paste0(moddir,"LRT_mod_RT_nna_valence_action.rds"))
anova(mod_LRT)
#                      Df   Chisq Chi Df       Pr(>Chisq)    
# stim_go_f            14 25.6436      1 0.00000041065524 ***
# stim_win_f           14 44.5756      1 0.00000000002447 ***
# stim_go_f:stim_win_f 14  1.5135      1           0.2186 
# --> GOES INTO PAPER 

## Plots:
plot(effect("stim_win_f",mod)) # faster for Win
plot(effect("stim_go_f",mod)) # faster fo Go
plot(effect("stim_go_f:stim_win_f",mod))

# ==================================================================================================================== #
#### Follow-ups: ####

lsmeans(mod, pairwise ~ stim_go_f * stim_win_f, adj = "bonferroni")
# stim_go_f stim_win_f lsmean     SE  df asymp.LCL asymp.UCL
# Go        Avoid       0.802 0.0145 Inf     0.773     0.830
# NoGo      Avoid       0.865 0.0137 Inf     0.838     0.892
# Go        Win         0.703 0.0124 Inf     0.679     0.727
# NoGo      Win         0.748 0.0173 Inf     0.714     0.782
# 
# contrast              estimate     SE  df z.ratio p.value
# Go,Avoid - NoGo,Avoid  -0.0631 0.0121 Inf -5.207  <.0001 
# Go,Avoid - Go,Win       0.0986 0.0115 Inf  8.584  <.0001 
# Go,Avoid - NoGo,Win     0.0535 0.0147 Inf  3.644  0.0016 
# NoGo,Avoid - Go,Win     0.1618 0.0147 Inf 11.001  <.0001 
# NoGo,Avoid - NoGo,Win   0.1167 0.0157 Inf  7.451  <.0001 
# Go,Win - NoGo,Win      -0.0451 0.0111 Inf -4.049  0.0003 

# ==================================================================================================================== #
#### Merge both factors in one: ####

### Combine cue valence and required action into one factor:
mydata$stim_go_win_f <- factor(paste(mydata$stim_go_o,mydata$stim_win_o,sep = "2"))

## Mixed-effects linear regression: 
mod <- lmer(RT_nna ~ stim_go_win_f + (stim_go_win_f|PPN_f), data = mydata)

## P-values with chisquare tests:
Anova(mod,type = "III")
#                 Chisq Df            Pr(>Chisq)    
# (Intercept)   4172.69  1 < 0.00000000000000022 ***
# stim_go_win_f  128.32  3 < 0.00000000000000022 ***

## Plot:
plot(effect("stim_go_win_f",mod))
# fastest Go2Win, slowest NoGo2Avoid (few Go responses very late)

# Follow-ups:
lsmeans(mod, pairwise ~ stim_go_win_f)
lsmeans(mod, pairwise ~ stim_go_win_f, adj = "bonferroni")
summary(glht(mod, linfct = mcp(stim_go_win_f = "Tukey")))
# --> conclusion: all conditions significantly different from each other
summary(effect("stim_go_win_f", mod, se=TRUE))

# ==================================================================================================================== #
#### On subset of people used in EEG-fMRI analyses (N = 30) ####

## Select subsect of data;
subsetdata <- subset(mydata, !(PPN_n %in% c(1,11,15,19,21,25)))
length(unique(subsetdata$PPN_n)) # 29
sort(unique(subsetdata$PPN_n))

## Mixed-effects linear regression:
mod <- lmer(RT_nna ~ stim_go_f*stim_win_f + (stim_go_f*stim_win_f|PPN_f), data = subsetdata,na.action = na.exclude)
saveRDS(mod,paste0(moddir,"mod_RT_nna_valence_action_withoutInvalid.rds"))
mod <- readRDS(paste0(moddir,"mod_RT_nna_valence_action_withoutInvalid.rds"))
summary(mod)
#                         Estimate Std. Error        df t value             Pr(>|t|)    
# (Intercept)             0.776849   0.013389 28.931106  58.022 < 0.0000000000000002 ***
# stim_go_f1             -0.028010   0.005186 26.646369  -5.401        0.00001082882 ***
# stim_win_f1             0.054687   0.006445 29.195069   8.485        0.00000000226 ***
# stim_go_f1:stim_win_f1 -0.006434   0.004279 28.244948  -1.504                0.144    

## P-values with LRTs:
mod_LRT <- mixed(mod, subsetdata, method = "LRT", type = "III", control = lmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
saveRDS(mod_LRT,paste0(moddir,"LRT_mod_RT_nna_valence_withoutInvalid.rds"))
mod_LRT <- readRDS(paste0(moddir,"LRT_mod_RT_nna_valence_withoutInvalid.rds"))
anova(mod_LRT)
#                      Df   Chisq Chi Df     Pr(>Chisq)    
# stim_go_f            14 21.0184      1 0.000004548886 ***
# stim_win_f           14 37.3140      1 0.000000001006 ***
# stim_go_f:stim_win_f 14  2.2511      1         0.1335   
# --> GOES INTO SUPPLEMENTARY MATERIAL 

# ==================================================================================================================== #
#### Robustness checks: ####

# ======================================================= #
### Only session 1:

## Mixed-effects linear regression:
mod <- lmer(RT_nna ~ stim_go_f*stim_win_f + (stim_go_f*stim_win_f|PPN_f), data = mydata_ses1,na.action = na.exclude)
Anova(mod, type = "III")
#                          Chisq Df            Pr(>Chisq)    
# (Intercept)          3325.2753  1 < 0.00000000000000022 ***
# stim_go_f              23.3138  1           0.000001376 ***
# stim_win_f             83.0905  1 < 0.00000000000000022 ***
# stim_go_f:stim_win_f    0.0016  1                0.9686  

## P-values with LRTs:
mod_LRT <- mixed(mod, mydata_ses1, method = "LRT", type = "III", control = lmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#                      Df   Chisq Chi Df       Pr(>Chisq)    
# stim_go_f            14 18.7241      1 0.00001510586817 ***
# stim_win_f           14 43.4369      1 0.00000000004379 ***
# stim_go_f:stim_win_f 14  0.0011      1           0.9732 

# ======================================================= #
### Only session 2:

### Mixed-effects linear regression:
mod <- lmer(RT_nna ~ stim_go_f*stim_win_f + (stim_go_f*stim_win_f|PPN_f), data = mydata_ses2,na.action = na.exclude)
Anova(mod, type = "III")
#                         Chisq Df            Pr(>Chisq)    
# (Intercept)          4177.105  1 < 0.00000000000000022 ***
# stim_go_f              21.058  1      0.00000445478357 ***
# stim_win_f             42.229  1      0.00000000008117 ***
# stim_go_f:stim_win_f    0.782  1                0.3765  

## P-values with LRTs:
mod_LRT <- mixed(mod, mydata_ses2, method = "LRT", type = "III", control = lmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#                      Df   Chisq Chi Df   Pr(>Chisq)    
# stim_go_f            14 17.0226      1 0.0000369367 ***
# stim_win_f           14 27.7955      1 0.0000001348 ***
# stim_go_f:stim_win_f 14  0.7931      1       0.3732   

# ======================================================= #
### Only correct Go:

## Mixed-effects linear regression:
mod <- lmer(RT_nna ~ stim_win_f + (stim_win_f|PPN_f), data = mydata_cor_go,na.action = na.exclude)

## P-values with chisquare tests:
Anova(mod, type = "III")
#                Chisq Df            Pr(>Chisq)    
# (Intercept) 2831.099  1 < 0.00000000000000022 ***
# stim_win_f    17.184  1            0.00003393 ***

## P-values with LRTs:
mod_LRT <- mixed(mod, mydata_cor_go, method = "LRT", type = "III", control = lmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#            Df  Chisq Chi Df Pr(>Chisq)    
# stim_win_f  5 13.607      1  0.0002253 ***

# ======================================================= #
### Only incorrect Go:

## Mixed-effects linear regression:
mod <- lmer(RT_nna ~ stim_win_f + (stim_win_f|PPN_f), data = mydata_incor_go,na.action = na.exclude)

## P-values with chisquare tests:
Anova(mod, type = "III")
#                Chisq Df            Pr(>Chisq)    
# (Intercept) 3523.464  1 < 0.00000000000000022 ***
# stim_win_f    79.078  1 < 0.00000000000000022 ***

## P-values with LRTs:
mod_LRT <- mixed(mod, mydata_incor_go, method = "LRT", type = "III", control = lmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#            Df  Chisq Chi Df       Pr(>Chisq)    
# stim_win_f  5 42.563      1 0.00000000006844 ***

# ====================================================================================================================== #
#### 3b) RTs as function of valence and required action ####

## Mixed-effects linear regression:
mod <- lmer(RT ~ stim_go_f*stim_win_f + (stim_go_f*stim_win_f|PPN_f), data = mydata,na.action = na.exclude)
saveRDS(mod,paste0(moddir,"mod_RT_valence_action.rds"))
mod <- readRDS(paste0(moddir,"mod_RT_valence_action.rds"))
summary(mod)
# Fixed effects:
#                         Estimate Std. Error        df t value             Pr(>|t|)    
# (Intercept)             0.833753   0.017723 35.061120  47.045 < 0.0000000000000002 ***
# stim_go_f1             -0.029534   0.005934 31.778296  -4.977       0.000021585498 ***
# stim_win_f1             0.066877   0.007605 35.469282   8.794       0.000000000195 ***
# stim_go_f1:stim_win_f1 -0.007129   0.004333 33.372593  -1.645                0.109  

## P-values with chisquare tests:
Anova(mod, type = "III")
#                          Chisq Df            Pr(>Chisq)    
# (Intercept)          2213.2150  1 < 0.00000000000000022 ***
# stim_go_f              24.7734  1          0.0000006448 ***
# stim_win_f             77.3413  1 < 0.00000000000000022 ***
# stim_go_f:stim_win_f    2.7071  1                0.0999 .  

## P-values with LRTs:
mod_LRT <- mixed(mod, mydata, method = "LRT", type = "III", control = lmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
saveRDS(mod_LRT,paste0(moddir,"LRT_mod_RT_valence_action.rds"))
mod_LRT <- readRDS(paste0(moddir,"LRT_mod_RT_valence_action.rds"))
anova(mod_LRT)
#                      Df   Chisq Chi Df      Pr(>Chisq)    
# stim_go_f            14 19.5119      1 0.0000099976761 ***
# stim_win_f           14 41.7316      1 0.0000000001047 ***
# stim_go_f:stim_win_f 14  2.6756      1          0.1019    
se.fixef(mod)
# --> note: SEs look reasonable, although NoGo-RTs are based on (only a few) errors

## Plots:
plot(effect("stim_win_f",mod)) # faster for Win
plot(effect("stim_go_f",mod)) # faster fo Go
plot(effect("stim_go_f:stim_win_f",mod))
# --> faster for win, faster for Go

# ====================================================================================================================== #
#### Robustness checks: ####

# ======================================================= #
### Only session 1:

## Mixed-effects linear regression:
mod <- lmer(RT ~ stim_go_f*stim_win_f + (stim_go_f*stim_win_f|PPN_f), data = mydata_ses1,na.action = na.exclude)

## P-values with chisquare tests:
Anova(mod, type = "III")
#                          Chisq Df            Pr(>Chisq)    
# (Intercept)          1331.3319  1 < 0.00000000000000022 ***
# stim_go_f              11.2670  1              0.000789 ***
# stim_win_f             68.4404  1 < 0.00000000000000022 ***
# stim_go_f:stim_win_f    1.1377  1              0.286132 

## P-values with LRTs:
mod_LRT <- mixed(mod, mydata_ses1, method = "LRT", type = "III", control = lmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#                      Df   Chisq Chi Df      Pr(>Chisq)    
# stim_go_f            14 10.1898      1        0.001412 ** 
# stim_win_f           14 38.4938      1 0.0000000005493 ***
# stim_go_f:stim_win_f 14  1.1385      1        0.285967  

# ======================================================= #
### Only session 2:

## Mixed-effects linear regression:
mod <- lmer(RT ~ stim_go_f*stim_win_f + (stim_go_f*stim_win_f|PPN_f), data = mydata_ses2,na.action = na.exclude)

## P-values with chisquare tests:
Anova(mod, type = "III")
#                          Chisq Df            Pr(>Chisq)    
# (Intercept)          3376.5443  1 < 0.00000000000000022 ***
# stim_go_f              16.8022  1       0.0000414859143 ***
# stim_win_f             40.7835  1       0.0000000001701 ***
# stim_go_f:stim_win_f    0.9979  1                0.3178 

## P-values with LRTs:
mod_LRT <- mixed(mod, mydata_ses2, method = "LRT", type = "III", control = lmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#                      Df  Chisq Chi Df   Pr(>Chisq)    
# stim_go_f            14 14.459      1    0.0001432 ***
# stim_win_f           14 26.956      1 0.0000002082 ***
# stim_go_f:stim_win_f 14  1.004      1    0.3163338 

# ======================================================= #
### Only correct Go:

## Mixed-effects linear regression:
mod <- lmer(RT ~ stim_win_f + (stim_win_f|PPN_f), data = mydata_cor_go,na.action = na.exclude)

## P-values with chisquare tests:
Anova(mod, type = "III")
#                Chisq Df          Pr(>Chisq)    
# (Intercept) 936.2159  1 <0.0000000000000002 ***
# stim_win_f    1.9901  1              0.1583

## P-values with LRTs:
mod_LRT <- mixed(mod, mydata_cor_go, method = "LRT", type = "III", control = lmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#            Df  Chisq Chi Df Pr(>Chisq)    
# stim_win_f  5 1.9096      1      0.167

# ======================================================= #
### Only incorrect Go:

## Mixed-effects linear regression:
mod <- lmer(RT ~ stim_win_f + (stim_win_f|PPN_f), data = mydata_incor_go,na.action = na.exclude)

## P-values with chisquare tests:
Anova(mod, type = "III")
#                Chisq Df            Pr(>Chisq)    
# (Intercept) 1453.120  1 < 0.00000000000000022 ***
# stim_win_f    75.469  1 < 0.00000000000000022 ***

## P-values with LRTs:
mod_LRT <- mixed(mod, mydata_incor_go, method = "LRT", type = "III", control = lmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#            Df  Chisq Chi Df      Pr(>Chisq)    
# stim_win_f  5 41.211      1 0.0000000001367 **

# ====================================================================================================================== #
#### 3c) RTs (corrected) as function of correct Go/ incorrect Go/ NoGo ####

## Create variable reflecting correct Go/ incorrect Go/ NoGo
mydata$stim_go_acc <- factor(ifelse(mydata$stim_go_f == "NoGo", "NoGo",
                                    ifelse(mydata$stim_go_f == "Go" & mydata$ACC_cor_n == 1,"Correct Go",
                                           ifelse(mydata$stim_go_f == "Go" & mydata$ACC_cor_n == 0,"Incorrect Go",
                                                  NA))))
sum(is.na(mydata$stim_go_acc)) # No NAs
table(mydata$stim_go_acc)
# correct Go incorrect Go         NoGo 
#       4303         7217        11520 

## Plot:
ggplot(mydata,aes(stim_go_acc, RT_nna)) + stat_summary(fun.y = mean, geom = "bar", position = "dodge") + stat_summary(fun.data=mean_cl_normal, geom="errorbar", position=position_dodge(width=0.9), width=.1) +
  labs(x="Response and accuracy", y="Reaction time") + 
  ggtitle("Reaction time as a function action and accuracy") +  
  theme_classic() + # GG_pGo_valence_action_ACC
  theme(axis.text=element_text(size=15),axis.title=element_text(size=20), title = element_text(size=15), legend.text = element_text(size=15))

## Mixed-effects linear regression:
mod <- lmer(RT ~ stim_go_acc + (stim_go_acc|PPN_f), data = mydata,na.action = na.exclude)
summary(mod)

## P-values with chisquare tests:
Anova(mod, type = "III")
#                Chisq Df            Pr(>Chisq)    
# (Intercept) 1950.923  1 < 0.00000000000000022 ***
# stim_go_acc   25.061  2           0.000003615 ***

## P-values with LRTs:
mod_LRT <- mixed(mod, mydata, method = "LRT", type = "III", control = lmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#             Df  Chisq Chi Df  Pr(>Chisq)    
# stim_go_acc  8 19.441      2 0.00006004 ***

## Plot:
plot(effect("stim_go_acc",mod))
# --> correct Go faster than incorrect Go and NoGo, but those not different from each other

# ====================================================================================================================== #
#### 3d) RTs for left vs. right Go responses: ####

max(mydata$RT_nna,na.rm=T) # 1.3035 --> maximum otherwise NaN
tapply(mydata$RT,mydata$response_f,mean,na.rm=T)
#      Left      NoGo     Right 
# 0.8280673       NaN 0.8063441
tapply(mydata$RT_nna,mydata$response_f,mean,na.rm=T)
#      Left      NoGo     Right 
# 0.7729851       NaN 0.7447934

## With RT:
mod <- lmer(RT ~ response_f + (response_f|PPN_f), data = mydata)
summary(mod)
#              Estimate Std. Error        df t value            Pr(>|t|)    
# (Intercept)  0.820933   0.018297 34.989217  44.867 <0.0000000000000002 ***
# response_f1  0.012120   0.006499 34.725716   1.865              0.0707 . 
mod_LRT <- mixed(mod, mydata, method = "LRT", type = "III",
                 control = lmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#            Df Chisq Chi Df Pr(>Chisq)  
# response_f  5 3.417      1    0.06453 .

## With RT_nna:
mod <- lmer(RT_nna ~ response_f + (response_f|PPN_f), data = mydata)
summary(mod)
#              Estimate Std. Error        df t value            Pr(>|t|)    
# (Intercept)  0.767165   0.012269 34.974231  62.529 <0.0000000000000002 ***
# response_f1  0.013982   0.005225 34.683454   2.676              0.0113 * 
mod_LRT <- mixed(mod, mydata, method = "LRT", type = "III",
                 control = lmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#            Df  Chisq Chi Df Pr(>Chisq)   
# response_f  5 6.7089      1   0.009593 **
# --> GOES INTO SUPPLEMENTARY MATERIAL

# ====================================================================================================================== #
#### 3e) RTs for correct vs incorrect Go: ####

tapply(mydata$RT,mydata$ACC_f,mean,na.rm=T)
#   Correct Incorrect 
# 0.7840502 0.8412819 

mod <- lmer(RT ~ ACC_f + (ACC_f|PPN_f), mydata)
summary(mod)
#              Estimate Std. Error        df t value             Pr(>|t|)    
# (Intercept)  0.819099   0.018473 35.044157  44.340 < 0.0000000000000002 ***
# ACC_f1      -0.022596   0.006523 34.690957  -3.464              0.00143 ** 
plot(effect("ACC_f",mod))

## P-values with LRTs:
mod_LRT <- mixed(mod, mydata, method = "LRT", type = "III",
                 control = lmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#       Df  Chisq Chi Df Pr(>Chisq)   
# ACC_f  5 10.624      1   0.001116 **

# ====================================================================================================================== #
#### 3f) RTs for correct vs Go but incorrect vs. Go instead of NoGo: ####
# - correct Go
# - incorrect Go, should have made another Go
# - incorrect Go, should have made NoGo

mydata$ACC_types_n <- ifelse(mydata$stim_go_n==1 & mydata$is_go_cor_n==1 & mydata$ACC_cor_n==1,1,
                             ifelse(mydata$stim_go_n==1 & mydata$is_go_cor_n==1 & mydata$ACC_cor_n==0,2,
                                    ifelse(mydata$stim_go_n==0 & mydata$is_go_cor_n==1 & mydata$ACC_cor_n==0,3,
                                           NA)))
table(mydata$ACC_types_n)
mydata$ACC_types_f <- factor(mydata$ACC_types_n, levels = c(1,2,3), labels = c("correct","wrong Go","Go for NoGo"))
tapply(mydata$RT,mydata$ACC_types_f,mean,na.rm = T)
#   correct    wrong Go Go for NoGo 
# 0.7273191   0.7497505   0.8169105 

## Model:
mod <- lmer(RT ~ ACC_types_f + (ACC_types_f|PPN_f), mydata)
summary(mod)
#               Estimate Std. Error        df t value             Pr(>|t|)    
# (Intercept)   0.780063   0.011330 35.107525  68.852 < 0.0000000000000002 ***
# ACC_types_f1 -0.037870   0.007483 34.594992  -5.061            0.0000137 ***
# ACC_types_f2 -0.024860   0.007222 34.721538  -3.442              0.00152 ** 
plot(effect("ACC_types_f",mod))

## P-values with LRTs:
mod_LRT <- mixed(mod, mydata, method = "LRT", type = "III",
                 control = lmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)

# -------------------------------------------------------------- #:
### Correct Go vs. other Go response:

mydata2 <- droplevels(subset(mydata, ACC_types_f %in% c("correct","wrong Go")))
levels(mydata2$ACC_types_f)
mod <- lmer(RT ~ ACC_types_f + (ACC_types_f|PPN_f), mydata2)
summary(mod)
#               Estimate Std. Error        df t value            Pr(>|t|)    
# (Intercept)   0.748778   0.012273 35.134039  61.008 <0.0000000000000002 ***
# ACC_types_f1 -0.006449   0.005416 34.506260  -1.191               0.242  

## P-values with LRTs:
mod_LRT <- mixed(mod, mydata2, method = "LRT", type = "III",
                 control = lmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#             Df  Chisq Chi Df Pr(>Chisq)
# ACC_types_f  5 1.4274      1     0.2322

# -------------------------------------------------------------- #:
### Other Go response vs. Go for NoGo:

mydata3 <- droplevels(subset(mydata, ACC_types_f %in% c("wrong Go","Go for NoGo")))
levels(mydata3$ACC_types_f)
mod <- lmer(RT ~ ACC_types_f + (ACC_types_f|PPN_f), mydata3)
summary(mod)
#               Estimate Std. Error        df t value             Pr(>|t|)    
# (Intercept)   0.799028   0.011750 35.346664  68.001 < 0.0000000000000002 ***
# ACC_types_f1 -0.043743   0.007847 33.500123  -5.574           0.00000323 ***

## P-values with LRTs:
mod_LRT <- mixed(mod, mydata3, method = "LRT", type = "III",
                 control = lmerControl(optCtrl = list(maxfun = 1e+9), optimizer = c("bobyqa")))
anova(mod_LRT)
#             Df  Chisq Chi Df  Pr(>Chisq)    
# ACC_types_f  5 23.501      1 0.000001249 ***

# END
