# Uses random forest feature selection to determine best features that affect dnds estimate accuracy

library(knitr)
library(caret)
library(plyr)
library(doMC)  # for parallele rfe

source("./load_all_sim_dnds.R")

dnds <- get_all_sim_dnds()
dim(dnds)
summary(dnds)
head(dnds)


NUM_RESP_NAMES <- c("LOD_dNdS", "Dist_dn_minus_dS", "Ratio_dNdS", "AbsLOD_dNdS", "AbsDist_dn_minus_dS")
CAT_RESP_NAMES <- c("CrapLOD", "CrapDist")
COVAR_NAMES <- colnames(dnds[sapply(dnds,is.numeric)])[!colnames(dnds[sapply(dnds,is.numeric)]) %in% NUM_RESP_NAMES]
CAT_COVAR_NAMES <- c("IsLowSubst.Act")
LM_COVAR_NAMES <- c(CAT_COVAR_NAMES, 
                    COVAR_NAMES[!(COVAR_NAMES %in% c("dNdS.Act", "dN_minus_dS.Act", "dN_minus_dS.Exp", 
                                                     # In separate analysis, conservation and entropy are highly correlated.
                                                     # When we use speedglm, it bugs out after it removes highly correlated variables.
                                                     # So we do it for them.
                                                     "ConserveTrueBase.Act", "ConserveTrueBase.Exp", "Window_Conserve.Act",
                                                     "UnambigCodonRate.Act", "Window_UnambigCodonRate.Act",
                                                     # These are highly correlated with N, S
                                                     "Subst.Act", "Subst.Exp",
                                                     "EN.Exp", "ES.Exp", "EN.Act", "ES.Act",
                                                     "Window_Start", "Window_End", "CodonSite", "Reads.Act"
                    )
                    )])

#  TODO:  remove this once we are confident about the function
respname <- "CrapLOD"
feats <- c(COVAR_NAMES, CAT_COVAR_NAMES)


# TODO:  do PCA to determine the total features we should keep 
# based on the amount of independent signals within them
# Or Parital Least Squares Analysis which finds the area under ROC for number of orthogonal components kept

FOLDS <- 3
SEED <- 7
PROCS <- 3



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
  print(predictors(results))  # results$optVariables also does the same)
  # plot the results
  plot(results, type=c("g", "o"))
  # per-variable importance
  print(varImp(results))
  
  # the time it took to finish
  print(results$times)
}


# Method 2: Random Forest, Recurisve Feature Elmination  (Backwards Selection)  for continuous response
rf_feat_sel_cont_rfe <- function(dnds, respname, feats) {
  
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
                 metric="RMSE"  # This is for continuous response
  )
  print(results)
  # list the chosen features
  print(predictors(results))  # results$optVariables also does the same
  # plot the results
  fig <- plot(results, type=c("g", "o"))
  print(fig)
  # per-variable importance
  print(varImp(results))
  
  # the time it took to finish
  print(results$times)
}


rf_feat_sel_class_rfe(dnds=dnds, respname="CrapLOD", feats=feats)

# Random forest recursive backwards feature selection for continous response 
#rf_feat_sel_cont_rfe(dnds=dnds, respname="AbsLOD_dNdS", feats=feats)
#rf_feat_sel_cont_rfe(dnds=dnds, respname="AbsDist_dn_minus_dS", feats=feats)











# 
# rf_feat_sel_class_cv <- function(dnds, respname, feats) {
#   print(feats)
#   
#   # Remove samples where response is NA
#   cleandnds <- dnds[!is.na(dnds[, respname]), ]
#   dim(cleandnds)
#   
#   # Take bags of samples with replacement
#   # Rows = number of training samples per bag, columns = number of bags aka number of partitions
#   partitionMat <-createResample(y=cleandnds[, respname], times=FOLDS, list=FALSE)  # ???? Nope, this seems to select same sample twice for same bag
#   print(dim(partitionMat))
#   print(summary(partitionMat))
#   totalUniqSamplesPerBag <- sapply(1:ncol(partitionMat), function(x) {length(unique(partitionMat[, x]))})
#   print(totalUniqSamplesPerBag)
#   
#   for (bag in 1:FOLDS) {
#     set.seed(SEED)
#     bagSamples <- partitionMat[, bag]
#     bag_dnds <- cleandnds[rownames(dnds) %in% bagSamples, ]    
#     outbag_dnds <- cleandnds[!(rownames(dnds) %in% bagSamples), ]
#     
#     trainFormula <- as.formula(paste0(respname, "~", paste0(feats, collapse=" + ")))
#     print(trainFormula)
#     
#     # Method 1: Cross Validation, then Average Importance Across Each Partition    
#     
#     
#     # Cross Validation Function
#     control <- trainControl(method="cv",number=FOLDS)
#     
#     fit <-train(trainFormula, data=bag_dnds, method="rf",
#                 trControl=control,
#                 prox=TRUE,allowParallel=TRUE)
#     print(fit)
#     
#     imp <- print(varImp(fit)) 
#   }  
# }



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
