#' preprocess_PEs_cueLocked.R
#' 
#' Load .txt file output from featquery, extract relevant summary statistic,
#' rearrange conditions, save for all subjects in one .csv file. 
#' 
#' Requires that Run_ROI/run_ROI_PEs_2ndlevel.sh has been run and .txt files 
#' are available under e.g. sub-001/GLM1/featquery_GLM1MedialCaudateValence
#'
#' Mind to adjust rootdir.
#'
#' INPUTS:
#' none
#'
#' OUTPUTS:
#' no outputs, saves .txt file to Log/fMRI/fMRI_ROIs/ROIPlots.
#'
#' EEG/fMRI STUDY, DONDERS INSTITUTE, NIJMEGEN.
#' J. Algermissen, 2018-2021.

# ========================================================================================================================================
#### General settings: ####

GLMID <- "1" # GLM where PEs were extracted from
GLMIDmask <- GLMID # GLM used for creating mask

nPE <- 4 # number of PEs to extract
PEselection <- c(1,3,2,4) # use for selecting and reordering PEs

# Valid subjects:
invalidSubs <- c(15,25) # subjects to exclude (no PEs available)
validSubs <- which(!(1:36 %in% invalidSubs))
nSub <- length(validSubs)

# ======================================================================================================================================== #
#### Create/ set directories: ####

rootdir <- "/project/3017042.02/"
plotdir <- paste0(rootdir,"Log/fMRI/fMRI_ROIs/ROIPlots/")

# Diretory dependent on GLMID:
outputdir <- paste0(rootdir,"Log/fMRI/fMRI_ROIs/GLM",GLMID,"_ROIData/")
if (!dir.exists(outputdir)) {dir.create(outputdir)}

# ======================================================================================================================================== #
#### Create wide data frame: ####

create_wide <- function(GLMID,ROIID){
  
  cat(paste0("Extract ROI ", ROIID, " from GLM",GLMID,": ",nSub, " subjects, ",nPE," PEs\n"))

  # 1) Initialize empty data frame:
  ROI_raw <- data.frame(matrix(0,nSub,nPE),stringsAsFactors=FALSE) 
  
  # 2) Loop through subjects and load:
  iFill <- 0; for (iSub in validSubs){ # iSub = 1;
    
    iFill <- iFill + 1; # counter for rows to fill, always increment by 1
    print(sprintf("Load subject %02d into row %02d",iSub,iFill))
    
    # Determine directory:
    datadirectory <- sprintf("/project/3017042.02/Log/fMRI/sub-%03d/GLM%s/GLM%s_sub%03d.feat/",iSub,GLMID,GLMID,iSub)
    
    report <- read.delim(paste0(datadirectory,"featquery_GLM",GLMIDmask,ROIID,"/report.txt"),sep = " ", header = F) # load
    names(report) <- c("X","PEname","nVoxels","min","10%","mean","median","90%","max","stddev","vox1","vox2","vox3","mm1","mm2","mm3")
    for (iCope in 1:nPE){
      ROI_raw[iFill,iCope] <- report$mean[iCope]
    }
  }
  
  # 3) Reorder PEs:
  cat(paste0("Select PEs ",
             paste(as.character(sort(PEselection)),collapse = " "),
             " in order ",
             paste(as.character(PEselection),collapse = " "),
             "\n"))
  ROI_wide <- ROI_raw[,PEselection]
  names(ROI_wide) <- paste0("PE",1:length(PEselection)) # add names
  
  # 4) Create subject ID:
  ROI_wide$PPN_n <- validSubs # add subject number
  ROI_wide$PPN_f <- as.factor(ROI_wide$PPN_n) # as factor
  ROI_wide$group <- 1 # fake x-axis for plots
  
  return(ROI_wide)
}

# ======================================================================================================================================== #
#### Loop over ROIs to use:

for (ROIID in c("MedialCaudateValence","LeftPutamenValence","JBLIncongruency")){
  
  # Load and extract:
  ROI_wide <- create_wide(GLMID,ROIID)
  
  # Write as .txt file for Matlab: 
  output <- t(as.matrix(ROI_wide[,c("PE1","PE2","PE3","PE4")]))
  write.table(output,paste0(outputdir,"GLM",GLMID,ROIID,".txt"),row.names = F, col.names = F)
  
}

cat("Finished :-)\n")

# END