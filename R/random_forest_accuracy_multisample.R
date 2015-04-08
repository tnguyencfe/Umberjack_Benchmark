# Uses random forest feature selection to determine best features that affect dnds estimate accuracy

library(caret)
library(plyr)
library(doMC)  # for parallele rfe

DNDS_FILENAME <- "/home/thuy/gitrepo/Umberjack_Benchmark/simulations/out/collate_all.small.csv"

dnds <- read.table(DNDS_FILENAME, header=TRUE, sep=",", na.strings=c("", "None"))
dim(dnds)
summary(dnds)
head(dnds)
dnds$ErrBaseRate.Act <- dnds$ErrBase.Act/dnds$Reads.Act  # per-window-codonsite-read nucleotide errors
dnds$AmbigPadBaseRate.Act <- dnds$AmbigPadBase/dnds$Reads.Act  # per-window-codonsite-read nucleotide N's or gaps
dnds$UnambigCodonRate.Act <- dnds$UnambigCodons.Act/ dnds$Reads.Act  # per-window-codonsite fraction of unambiguous codons
dnds$Subst.Act <- dnds$N.Act + dnds$S.Act
dnds$Subst.Exp <- dnds$N.Exp + dnds$S.Exp
dnds$Cov.Act <- dnds$Reads.Act/dnds$PopSize.Act
dnds <- subset(dnds, select=-c(PopSize.Act, ErrBase.Act, AmbigPadBase.Act, UnambigCodons.Act))
dim(dnds)
summary(dnds)
head(dnds)

# Average across all codon sites in a window
per_window_ave <- ddply(.data=dnds, .variables=c("File", "Window_Start"), 
                        .fun=function(x) {                            
                          data.frame(Window_Conserve.Act=mean(x$ConserveTrueBase.Act, na.rm=TRUE),
                                     Window_Entropy.Act=mean(x$EntropyTrueBase.Act, na.rm=TRUE),
                                     Window_UnambigCodonRate.Act=mean(x$UnambigCodonRate.Act, na.rm=TRUE),
                                     Window_ErrBaseRate.Act=mean(x$ErrBaseRate.Act, na.rm=TRUE))
                        })
dnds <- merge(x=dnds, y=per_window_ave, by=c("File", "Window_Start"), all=TRUE, sort=TRUE)
dim(dnds)
summary(dnds)
head(dnds)

# Now remove window-codonsites where there is no dnds and no dn-ds information because of insufficient window sequences
dnds <- dnds[!is.na(dnds$dN_minus_dS.Exp) & !is.na(dnds$N.Act), ]
dim(dnds)
summary(dnds)
head(dnds)

# Ignore dN/dS == 0 for numerical stability
dnds$LOD_dNdS <- log(dnds$dNdS.Act) - log(dnds$dNdS.Exp)
dnds$LOD_dNdS[dnds$dNdS.Exp == 0 | dnds$dNdS.Act == 0] <- NA
dnds$Dist_dn_minus_dS <- dnds$dN_minus_dS.Act - dnds$dN_minus_dS.Exp
dnds$CrapLOD <- as.factor(abs(dnds$LOD_dNdS) > 1)
dnds$CrapDist <- as.factor(abs(dnds$Dist_dn_minus_dS) > 1)
summary(dnds)
dim(dnds)
head(dnds)

NUM_RESP_NAMES <- c("LOD_dNdS", "Dist_dn_minus_dS")
CAT_RESP_NAMES <- c("CrapLOD", "CrapDist")
COVAR_NAMES <- colnames(dnds[sapply(dnds,is.numeric)])[!colnames(dnds[sapply(dnds,is.numeric)]) %in% NUM_RESP_NAMES]
LM_COVAR_NAMES <- COVAR_NAMES[!(COVAR_NAMES %in% c("dNdS.Act", "dN_minus_dS.Act", "dN_minus_dS.Exp", "Window_End"))]


#  TODO:  remove this once we are confident about the function
respname <- "CrapLOD"
feats <- LM_COVAR_NAMES


# TODO:  do PCA to determine the total features we should keep 
# based on the amount of independent signals within them
# Or Parital Least Squares Analysis which finds the area under ROC for number of orthogonal components kept

FOLDS <- 5
SEED <- 7
PROCS <- 2


rf_feat_sel_class_cv <- function(dnds, respname, feats) {
  print(feats)
  
  # Remove samples where response is NA
  cleandnds <- dnds[!is.na(dnds[, respname]), ]
  dim(cleandnds)
  
  # Take bags of samples with replacement
  # Rows = number of training samples per bag, columns = number of bags aka number of partitions
  partitionMat <-createResample(y=cleandnds[, respname], times=FOLDS, list=FALSE)  # ???? Nope, this seems to select same sample twice for same bag
  print(dim(partitionMat))
  print(summary(partitionMat))
  totalUniqSamplesPerBag <- sapply(1:ncol(partitionMat), function(x) {length(unique(partitionMat[, x]))})
  print(totalUniqSamplesPerBag)
  
  for (bag in 1:FOLDS) {
    set.seed(SEED)
    bagSamples <- partitionMat[, bag]
    bag_dnds <- cleandnds[rownames(dnds) %in% bagSamples, ]    
    outbag_dnds <- cleandnds[!(rownames(dnds) %in% bagSamples), ]
    
    trainFormula <- as.formula(paste0(respname, "~", paste0(feats, collapse=" + ")))
    print(trainFormula)
    
#     # Method 1: Cross Validation, then Average Importance Across Each Partition    
#     
#     
#     # Cross Validation Function
#     control <- trainControl(method="cv",number=FOLDS)
#     
#     fit <-train(trainFormula, data=bag_dnds, method="rf",
#                     trControl=control,
#                     prox=TRUE,allowParallel=TRUE)
#     print(fit)
#     
#     imp <- print(varImp(fit))
    
   
    
  }

  
  
  
  #
  
}

# Method 2: Random Forest, Recurisve Feature Elmination  (Backwards Selection)
rf_feat_sel_class_rfe <- function(dnds, respname, feats) {
  
  print(paste0("Response=", respname))
  print(paste0("Features=", paste0(feats, collapse=", ")))
  
  # Remove samples where response is NA
  cleandnds <- dnds[!is.na(dnds[, respname]), ]
  dim(cleandnds)
  summary(cleandnds)
  
  # Random Forest Selection Function
  # TODO:  should we rerank features after they are removed???
  # This will auto resample and create bags for us
  control <- rfeControl(functions=rfFuncs, method="boot", number=FOLDS)
  
  ## Note: if the underlying model also uses foreach, the
  ## number of cores specified above will double (along with
  ## the memory requirements)
  registerDoMC(cores = PROCS)
  
  # TODO:  this can be parallelized across each bag
  # This takes 6 minutes even for A bag only has 5000 samples.  You can get timing through results$timing
  results <- rfe(x=cleandnds[, feats], y=cleandnds[, respname], sizes=c(1:length(feats)),
                 rfeControl=control, 
                 metric="Accuracy"  # This is for classification.  use another metric for regression
  )
  print(results)
  # list the chosen features
  predictors(results)  # results$optVariables also does the same
  # plot the results
  plot(results, type=c("g", "o"))
  # per-variable importance
  varImp(results)
  
  # the time it took to finish
  print(results$times)
}


# Remove highly correlated features so that speedglm doesn't crash due to bug where it doesn't update size of response var 
# in linear sys of equations after it removes eigenvectors == 0 from data
#
# Find the correlation of all the possible variables that can go into linear modelling
# lm_var_cor <- cor(cleandnds[, LM_COVAR_NAMES], method="spearman")
# summary(lm_var_cor)
# dim(lm_var_cor)
# # find attributes that are highly corrected
# highlyCorrelated <- findCorrelation(lm_var_cor, cutoff=0.95, verbose=FALSE)
# length(highlyCorrelated)
# # print indexes of highly correlated attributes
# print(highlyCorrelated)
