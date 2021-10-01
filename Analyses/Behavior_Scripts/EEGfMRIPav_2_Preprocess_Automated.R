# EEGfMRIPav_2_Preprocess_Automated

EEGfMRIPav_2_Preprocess_Automated <- function(){
  #' EEGfMRIPav_2_Preprocess_Automated()
  #' 
  #' Reads data from EEGfMRIPav_all.csv (which contains already read-in and recoded behavioral variables),
  #' creates (refreshes) ordered factors,
  #' creates several more convenience variables.
  #' Essentialy does the pre-processing steps (recoding, deleting data) necessary to do other analyses.
  #'
  #' Requires that .csv files have been recoded and concatenated using 
  #' EEGfMRIPav_1_Read_Initial.R
  #'
  #' Mind to adjust rootdir.
  #'
  #' INPUTS:
  #' none
  #'
  #' OUTPUTS:
  #' @return data frame with original variables from EEGfMRIPav_all.csv plus additional convenience variables.
  #'
  #' EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
  #' J. Algermissen, 2018-2021.
  
  # ================================================================================================================================================ #
  #### Set directories: ####
  
  rootdir <- "/project/3017042.02/" # root directory--needs to be adapted to users' folder structure
  inputdir <- paste0(rootdir,"Log/Behavior/Behavior_concatenated/")

  # ================================================================================================================================================ #
  #### Load data: ####
  
  cat("Retrieve EEGfMRIPav_all.csv\n")
  data <- read.csv(paste0(inputdir,"EEGfMRIPav_all.csv"),header=T) # already recoded etc.
  
  # ================================================================================================================================================ #
  #### Process stimuli: ####
  
  ## 1) Stimulus number:
  cat("Absolute stimulus number (1-16)\n")
  data$stim_total_n <- ifelse(data$session_n == 1, data$stim_n, data$stim_n+8)
  
  ## 2) Ordered factors:
  cat("Stimulus factors into ordered factors\n")
  data$stim_win_o <- ordered(data$stim_win_f, levels = c("Win","Avoid"))
  data$stim_go_o <- ordered(data$stim_go_f, levels = c("Go","NoGo"))
  
  ## 3) Recode factors of stimulus type:
  cat("4 different stimulus types as separate variables\n")
  data$stim_type_f <- factor(paste(data$stim_go_f,data$stim_win_f,sep = "2"))
  data$stim_type_o <- factor(data$stim_type_f, levels=c("Go2Win", "Go2Avoid", "NoGo2Win","NoGo2Avoid"), ordered=TRUE)
  
  # ================================================================================================================================================ #
  #### Process left/ right/ NoGo responses: ####
  
  ## 4) Create dummy indicating whether RT was too fast:
  cat("Too fast responses\n")
  data$is_too_fast <- ifelse(!(is.na(data$RT) & data$RT < 0.2),1,0)
  
  ## 5) Create dummy indicating whether RT was too slow:
  cat("Too late responses\n")
  data$is_too_late <- ifelse(!(is.na(data$RT) & data$RT > 1.3035),1,0)
  
  ## 6) Alternative A: Set NoGos (no RT) or too long RTs (> 1.3035) to NoGo, otherwise Go:
  cat("Responses: Set too late responses to NoGo (_n and _f)\n")
  data$response_cor_n <- ifelse(is.na(data$RT) | data$RT > 1.3035,0,data$response_n)
  data$response_cor_f <- factor(ifelse(is.na(data$RT) | data$RT > 1.3035,"NoGo",data$response_f),labels=c("Left","NoGo","Right"))
  
  ## 7) Alternative B: Delete responses based on RT either NA or too early RT or too late RT (set to NA):
  cat("Responses: Delete too fast and too late responses to NA (_f and _n)\n")
  data$response_nna_n <- ifelse(!(is.na(data$RT)) & (data$RT < 0.2 | data$RT > 1.3035), NA, data$response_n)
  data$response_nna_f <- factor(ifelse(!(is.na(data$RT)) & (data$RT < 0.2 | data$RT > 1.3035), NA, data$response_f),labels=c("Left","NoGo","Right"))
  
  ## 8) Response index:
  cat("Response as index 1, 2, or 3\n")
  data$response_idx <- ifelse(data$response_cor_f == "Left", 1,
                                ifelse(data$response_cor_f == "Right", 2, 3))
  
  # ================================================================================================================================================ #
  #### Process Go/NoGo responses: ####
  
  ## 9) Only count actual non-responses (0) as NoGo, count everything else as Go:
  cat("Responses: new corrected Go/NoGo response variable\n")
  data$is_go_cor_n <- ifelse(data$response_cor_n==0,0,1)
  data$is_go_cor_f <- factor(ifelse(data$response_cor_n==0,"NoGo","Go"))
  
  ## 10) Conflict between performed action and bias-implied action:
  data$conflict_n <- ifelse(data$is_go_cor_n==data$stim_win_n,1,0)
  data$conflict_f <- factor(data$conflict_n, levels = c(1,0), labels = c("congruent","incongruent"))
  
  # ================================================================================================================================================ #
  #### Process RTs: ####
  
  ## 11) Delete reaction times that are > 1.3035
  cat("RTs: too late responses to NA\n")
  data$RT_nna <- ifelse(data$RT > 1.3035, NA,data$RT)
  
  # ================================================================================================================================================ #
  #### Process accuracy: ####
  
  ## 11) Alternative A: Compare required action to recoded Left/Right/NoGo responses (any left or any right button press):
  cat("Accuracy: new accuracy variable based on updated Left/Right/NoGo response variable\n")
  data$ACC_cor_n <- ifelse(data$req_action_f==data$response_cor_f,1,0)
  data$ACC_cor_f <- factor(data$ACC_cor_n, levels = c(1,0), labels = c("correct","incorrect"))
  
  ## 12) Alternative B: Delete accuracy for button presses that are technically invalid:
  cat("Accuracy: set uninstructed key presses to NA\n")
  data$ACC_nna_n <- ifelse(data$is_uninstructed_key==1,NA,data$ACC_n)
  
  ## 13) Compute accuracy only based on Go/NoGo
  cat("Accuracy: only based on Go/NoGo distinction\n")
  data$goACC_n <- ifelse(data$stim_go_n==data$is_go_cor_n,1,0)
  
  # ================================================================================================================================================ #
  #### Process outcomes: ####
  
  ## 14) Exact reward:
  cat("Exact outcome: interpret neutral outcomes\n")
  data$outcome_all_f <- NA
  data$outcome_all_f[data$stim_win_f=="Win" & data$outcome_f=="Reward"] <- "Reward"
  data$outcome_all_f[data$stim_win_f=="Win" & data$outcome_f=="Neutral"] <- "No Reward"
  data$outcome_all_f[data$stim_win_f=="Avoid" & data$outcome_f=="Neutral"] <- "No Punishment"
  data$outcome_all_f[data$stim_win_f=="Avoid" & data$outcome_f=="Punishment"] <- "Punishment"
  data$outcome_all_f <- factor(data$outcome_all_f)

  data$outcome_all_o <- ordered(data$outcome_all_f, levels = 
                                    c("Reward","No Reward","No Punishment","Punishment"))

  ## 15) Outcome valence:
  cat("Valence of outcomes\n")
  data$outcome_valence_f <- factor(ifelse(data$outcome_all_f %in% c("Reward","No Punishment"),"Positive","Negative"))
  data$outcome_valence_o <- ordered(data$outcome_valence_f, levels = c("Positive","Negative"))
  
  ## 16) Outcome salience:
  cat("Salience of outcomes\n")
  data$outcome_salience_f <- factor(ifelse(data$outcome_all_f %in% c("Reward","Punishment"),"Salient","Neutral"))
  data$outcome_salience_o <- ordered(data$outcome_salience_f, levels = c("Salient","Neutral"))

  ## c) Combine action and outcome:
  cat("Combine actions and outcomes\n")
  data$is_go_outcome_all_f <- factor(paste(data$is_go_cor_f,data$outcome_all_f,sep = " "))
  data$is_go_outcome_all_o <- ordered(data$is_go_outcome_all_f, levels = 
                                       c("Go Reward","Go No Reward","Go No Punishment","Go Punishment",
                                         "NoGo Reward","NoGo No Reward","NoGo No Punishment","NoGo Punishment"))
  # ================================================================================================================================================ #
  #### Return: ####
  
  # Return data:
  cat("Return data :-)\n")
  return(data)

}
# END