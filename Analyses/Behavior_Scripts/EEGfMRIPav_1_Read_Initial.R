# EEGfMRIPav_1_Read_Initial

EEGfMRIPav_1_Read_Initial <- function(){
  #' EEGfMRIPav_1_Read_Initial()
  #' 
  #' Reads to .csv files obtained by executing EEGfMRIPav_extract_rawdata.m,
  #' give names to columns, recode required action and response, convert to factors, 
  #' save as EEGfMRIPav_all.csv
  #'
  #' Mind to adjust rootdir and targetdir.
  #'
  #' INPUTS:
  #' none
  #'
  #' OUTPUTS:
  #' @return saves all data of all subjects under EEGfMRIPav_all.csv in targetdir
  #'
  #' EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
  #' J. Algermissen, 2018-2021.

  # ================================================================================================================================================ 
  #### 1) Load and combine all dataframes ####
  
  ## Store current directory:
  currentdir <- getwd()

  ## Set directories:
  rootdir <- "/project/3017042.02/" # root directory--needs to be adapted to users' folder structure
  inputdir <- paste0(rootdir,"Log/Behavior/Data_beh_csv")
  targetdir <- paste0(rootdir,"Log/Behavior/Behavior_concatenated/")

  ## Load all data files and concatenate them 
  mydata <- do.call("rbind", lapply( list.files(inputdir,full=TRUE), read.csv, header=F))
  # should be 36*2*320 = 23040 entries of 10 variables

  # ================================================================================================================================================ 
  #### 2) Name variables; ####
  
  cat("Set variable names\n")
  names(mydata) <- c("PPN_n", "session_n", "trial_nr_n","stim_n","req_action_n","RT","ACC_n","response_n","outcome_n","is_go_n")

  # ================================================================================================================================================ 
  #### 3) Convert subject/ session/ trial variables to factors; ####
  
  cat("Participant, session, and trial number from numeric to factor\n")
  mydata$PPN_f <- factor(mydata$PPN_n)
  mydata$session_f <- factor(mydata$session_n)
  mydata$trial_nr_f <- factor(mydata$trial_nr_n)
  mydata$trial_nr_all_n <- (mydata$session_n -1 )*320 + mydata$trial_nr_n # trial number across sessions
  mydata$PPN_Ses_Trial_f <- factor(paste(mydata$PPN_f,mydata$session_f,mydata$trial_nr_n,sep="_"))
  
  cat("Create block variable\n")
  # if (!(mydata[iRow,"trial_nr_n"]) %in% c(110,220,320,430,540,640)){
    
  mydata$block_n <- ifelse(mydata$session_n==1 & mydata$trial_nr_n %in% 1:110,1,
                           ifelse(mydata$session_n==1 & mydata$trial_nr_n %in% 111:220,2,
                                  ifelse(mydata$session_n==1 & mydata$trial_nr_n %in% 221:320,3,
                                         ifelse(mydata$session_n==2 & mydata$trial_nr_n %in% 1:110,4,
                                                ifelse(mydata$session_n==2 &mydata$trial_nr_n %in% 111:220,5,
                                                       ifelse(mydata$session_n==2 & mydata$trial_nr_n %in% 221:320,6,NA))))))
  table(mydata$block_n)
  mydata$block_f <- factor(mydata$block_n)
  
  # ================================================================================================================================================ 
  #### 4) Convert stimulus features to factors: ####
  
  cat("Numeric stimulus variables to factors\n")
  mydata$stim_win_f <- factor(ifelse(mydata$stim_n %in% c(1,2,5,6),"Win","Avoid"))
  mydata$stim_win_n <- ifelse(mydata$stim_n %in% c(1,2,5,6),1,0)
  mydata$stim_go_f <- factor(ifelse(mydata$stim_n %in% 1:4,"Go","NoGo"))
  mydata$stim_go_n <- ifelse(mydata$stim_n %in% 1:4,1,0)
  mydata$stim_cong_f <- factor(ifelse(mydata$stim_win_n == mydata$stim_go_n,"congruent","incongruent"))
  mydata$stim_cong_n <- ifelse(mydata$stim_win_n == mydata$stim_go_n,1,0)

  ## Ordinal factors:
  cat("Standard factors to ordered factors\n")
  mydata$stim_win_o <- ordered(mydata$stim_win_f, levels = c("Win","Avoid"))
  mydata$stim_go_o <- ordered(mydata$stim_go_f, levels = c("Go","NoGo"))
  newdata <- mydata[order(mydata$PPN_n,mydata$session_n,mydata$trial_nr_n),]
  # check: might have to redo every time data is loaded, 
  # use EEGfMRIPav_2_Preprocess_Automated.R
  
  # ================================================================================================================================================ 
  #### 5) Recode required actions (look up in pavParams.m: 101 = left, 97 = right): ####
  
  cat("Required action from numeric to factor\n")
  mydata$req_action_f <- factor(ifelse(mydata$req_action_n == 101,"Left",
                                       ifelse(mydata$req_action_n == 97,"Right","NoGo")))
  
  # ================================================================================================================================================ 
  #### 6) Recode recorded response, convert into factor: ####

  cat("Uninstructed keys\n")
  mydata$is_uninstructed_key <- ifelse(mydata$response_n %in% c(0,97,101),0,1)
  
  cat("Recode instructed keys to respective correct key of same response side\n")
  mydata$response_n <- ifelse(mydata$response_n %in% c(101:104,69:72),101,
                              ifelse(mydata$response_n %in% c(97:100,65:68),97,
                                     ifelse(mydata$response_n == 0,0,
                                            0)))

  cat("Performed response from numeric to factor\n")
  mydata$response_f <- factor(ifelse(mydata$response_n==101,"Left",
                                     ifelse(mydata$response_n==97,"Right",
                                            ifelse(mydata$response_n==0,"NoGo",NA))))
  
  # ================================================================================================================================================ 
  #### 7) Convert accuracy to factor: ####
  
  cat("Accuracy from numeric to factor\n")
  mydata$ACC_f <- factor(ifelse(mydata$ACC_n == 1,"Correct","Incorrect"))
  
  # ================================================================================================================================================ 
  #### 8) Convert outcome to factor: ####
  
  cat("Outcome from numeric to factor\n")
  mydata$outcome_f <- factor(ifelse(mydata$outcome_n==1,"Reward",
                                    ifelse(mydata$outcome_n==0,"Neutral",
                                           "Punishment")))
  
  # ==========================================================================================
  #### 9) Reorder rows: ####
  
  cat("Reorder rows based on participant number, session number, trial number\n")
  newdata <- mydata[order(mydata$PPN_n,mydata$session_n,mydata$trial_nr_n),]
  
  # ==========================================================================================
  #### 10) Save: ####
  
  ## Save concatenated data:
  cat("Save data as EEGfMRIPav_all.csv in ",targetdir,"\n")
  write.csv(newdata,paste0(targetdir,"EEGfMRIPav_all.csv"),row.names = F)
  
}
# END