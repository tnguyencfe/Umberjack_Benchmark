library(plyr)


PSEUDOCOUNT <- 1e-7

get_all_sim_dnds <- function() {
  DNDS_FILENAME <- "/home/thuy/gitrepo/Umberjack_Benchmark/simulations/out/collate_all.small.csv"
  #DNDS_FILENAME <- "/home/thuy/gitrepo/Umberjack_Benchmark/simulations/out/collate_all.med.csv"
  
  dnds <- read.table(DNDS_FILENAME, header=TRUE, sep=",", na.strings=c("", "None"))
  dim(dnds)
  summary(dnds)
  head(dnds)
  #dnds$ErrBaseRate.Act <- dnds$ErrBase.Act/dnds$Reads.Act  # per-window-codonsite-read nucleotide errors
  #dnds$AmbigPadBaseRate.Act <- dnds$AmbigPadBase/dnds$Reads.Act  # per-window-codonsite-read nucleotide N's or gaps
  #dnds$UnambigCodonRate.Act <- dnds$UnambigCodons.Act/ dnds$Reads.Act  # per-window-codonsite fraction of unambiguous codons
  dnds$Subst.Act <- dnds$N.Act + dnds$S.Act
  dnds$Subst.Exp <- dnds$N.Exp + dnds$S.Exp
  #dnds$Cov.Act <- dnds$Reads.Act/dnds$PopSize.Act
  #dnds$IsLowSubst.Act <- as.factor(dnds$N.Act < 1 | dnds$S.Act < 1)
  #dnds <- subset(dnds, select=-c(PopSize.Act, ErrBase.Act, AmbigPadBase.Act, UnambigCodons.Act))
  dim(dnds)
  summary(dnds)
  
  
  
  # Average across all codon sites in a window
  per_window_ave <- ddply(.data=dnds, .variables=c("File", "Window_Start"), 
                          .fun=function(x) {                            
#                             data.frame(Window_Conserve.Act=mean(x$ConserveTrueBase.Act, na.rm=TRUE),
#                                        Window_Entropy.Act=mean(x$EntropyTrueBase.Act, na.rm=TRUE),
#                                        Window_UnambigCodonRate.Act=mean(x$UnambigCodonRate.Act, na.rm=TRUE),
#                                        Window_ErrBaseRate.Act=mean(x$ErrBaseRate.Act, na.rm=TRUE))
                            data.frame(Window_Entropy.Act=mean(x$EntropyTrueBase.Act, na.rm=TRUE),
                                       Window_UnambigCodonRate.Act=mean(x$UnambigCodonRate.Act, na.rm=TRUE),
                                       Window_ErrBaseRate.Act=mean(x$ErrBaseRate.Act, na.rm=TRUE))
                          })
  dnds <- merge(x=dnds, y=per_window_ave, by=c("File", "Window_Start"), all=TRUE, sort=TRUE)
  dim(dnds)
  summary(dnds)
  
  
  
  # Now remove window-codonsites where there is no dnds and no dn-ds information because of insufficient window sequences
  # NOOOO:  don't remove because this forms part of the crappy prediction
  #dnds <- dnds[!is.na(dnds$dN_minus_dS.Exp) & !is.na(dnds$N.Act), ]
  dnds <- dnds[!is.na(dnds$dNdS.Exp), ]
  dim(dnds)
  summary(dnds)
  
  
  # Ignore dN/dS == 0 for numerical stability
  dnds$LOD_dNdS <- log(dnds$dNdS.Act + PSEUDOCOUNT) - log(dnds$dNdS.Exp + PSEUDOCOUNT)
  dnds$AbsLOD_dNdS <- abs(dnds$LOD_dNdS)

  dnds$CrapLOD <- FALSE
  dnds$CrapLOD <- as.factor(dnds$AbsLOD_dNdS >= 1)
  dnds$CrapLOD[is.na(dnds$dNdS.Act)] <- TRUE

  dnds$wrongSelect <- FALSE
  dnds$wrongSelect <- as.factor((dnds$dNdS.Act < 1 & dnds$dNdS.Exp > 1) | (dnds$dNdS.Act < 1 & dnds$dNdS.Exp > 1))
  dnds$wrongSelect[is.na(dnds$dNdS.Act)] <- TRUE
  summary(dnds)
  dim(dnds)
  head(dnds)
  
  return (dnds)
}



# Random Forest, Recursive Feature Elmination  (Backwards Selection)  for Classification
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
  
 
  # TODO:  this can be parallelized across each bag
  # This takes 6 minutes even for A bag only has 5000 samples.  You can get timing through results$timing
  rfe_class_results <- rfe(x=cleandnds[, feats], y=cleandnds[, respname], sizes=c(1:length(feats)),
                           rfeControl=control, 
                           metric="Accuracy"  # This is for classification.  use another metric for regression
  )
  
  # save the classifier feature selection results (including its fit) to file 
  # so that we can load the environment variable back again later.
  save(rfe_class_results, file="rfe_class_results.RData")
  
  print(rfe_class_results)
  # list the chosen features
  print(predictors(rfe_class_results))  # results$optVariables also does the same)
  # plot the results
  plot(rfe_class_results, type=c("g", "o"))
  # per-variable importance
  print(varImp(rfe_class_results))
  
  # the time it took to finish
  print(rfe_class_results$times)
  return (rfe_class_results)
}


# Random Forest, Recursive Feature Elmination  (Backwards Selection)  for continuous response
rf_feat_sel_cont_rfe <- function(dnds, respname, feats) {
  
  print(paste0("Response=", respname))
  print(paste0("Features=", paste0(feats, collapse=", ")))
  
  # Remove samples where response is NA
  cleandnds <- dnds[!is.na(dnds[, respname]), ]
  print("Cleaned dnds dimensions:")
  print(dim(cleandnds))
  print("Cleaned dnds summary:")
  print(summary(cleandnds))
  
  
  new_rfFuncs <- rfFuncs
  new_rfFuncs$fit <- 
    function (x, y, first, last, ...) 
    {
      library(randomForest)
      randomForest(x, y, importance = first, ntree=501, ...)
    }
  # Random Forest Selection Function
  # This will auto resample and create bags for us.
  # Uses default number of trees = 500.  
  # TODO:  But we want to be able to break ties, so use odd number.
  control <- rfeControl(functions=new_rfFuncs, method="boot", number=FOLDS, 
                        #rerank=TRUE,  # rerank features after eliminate
                        #saveDetails=TRUE,   # save predictions and variable importances from selection process
                        #returnResamp="all",   # save all resampling summary metrics
                        verbose=TRUE
                        )
    
  # TODO:  this can be parallelized across each bag
  # This takes 6 minutes even for A bag only has 5000 samples.  You can get timing through results$timing
  rfe_cont_results <- rfe(x=cleandnds[, feats], y=cleandnds[, respname], 
                          #sizes=c(1:length(feats)),
                          sizes=c(1:2),
                          rfeControl=control, 
                          metric="RMSE"  # This is for continuous response
                          )
  
  # Save the rfe_cont_results environment object to file.
  save(rfe_cont_results, file="rfe_cont_results.RData")
  
  
  print(rfe_cont_results)
  # list the chosen features
  print(paste0("Predictors\n", predictors(rfe_cont_results)))  # results$optVariables also does the same
  # plot the results
  fig <- plot(rfe_cont_results, type=c("g", "o"))
  print(fig)
  # per-variable importance
  print("Variable Importance")
  print(varImp(rfe_cont_results))
  
  # the time it took to finish
  print("Timing")
  print(rfe_cont_results$times)
  
  return (rfe_cont_results)
}


# Cross Validation for random forest
rf_cv <- function(dnds, respname, feats) {
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
    
    # Method 1: Cross Validation, then Average Importance Across Each Partition    
    
    
    # Cross Validation Function
    control <- trainControl(method="cv",number=FOLDS)
    
    fit <-train(trainFormula, data=bag_dnds, method="rf",
                trControl=control,
                prox=TRUE,allowParallel=TRUE)
    print(fit)
    
    imp <- print(varImp(fit)) 
  }  
}

# Train random forest
rf_train <- function(dnds, respname, feats) {
  print(feats)
  
  # Remove samples where response is NA
  cleandnds <- dnds[!is.na(dnds[, respname]), ]
  dim(cleandnds)
  
  trainFormula <- as.formula(paste0(respname, "~", paste0(feats, collapse=" + ")))
  print(trainFormula)
  
  rf_output=randomForest(formula=trainFormula, data=cleandnds, importance = TRUE, ntree = 1001, proximity=TRUE, keep.forest=TRUE,
                         na.action=na.exclude)
  
  save(rf_output, file="RF_model")
  load("RF_model")
  
  rf_importances=importance(rf_output, scale=FALSE)
  print(rf_importances)
  print(rf_output$confusion)
  
  # i'th elelemtn is error rate for all trees up to the ith tree
  overall_error=rf_output$err.rate[length(rf_output$err.rate[,"OOB"]), "OOB"]*100
  print(overall_error)
  overall_accuracy=100-overall_error
  print(overall_accuracy)
  return (rf_output)
}