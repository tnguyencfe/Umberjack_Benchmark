library(plyr)


PSEUDOCOUNT <- 1e-7

get_all_sim_dnds <- function() {
  #DNDS_FILENAME <- "../simulations/out/collate_all.small.csv"
  DNDS_FILENAME <- "../simulations/out/collate_all.med.csv"
  
  dnds <- read.table(DNDS_FILENAME, header=TRUE, sep=",", na.strings=c("", "None"))
  dim(dnds)
  summary(dnds)
  head(dnds)
  #dnds$ErrBaseRate.Act <- dnds$ErrBase.Act/dnds$Reads.Act  # per-window-codonsite-read nucleotide errors
  #dnds$AmbigPadBaseRate.Act <- dnds$AmbigPadBase/dnds$Reads.Act  # per-window-codonsite-read nucleotide N's or gaps
  #dnds$UnambigCodonRate.Act <- dnds$UnambigCodons.Act/ dnds$Reads.Act  # per-window-codonsite fraction of unambiguous codons
  dnds$Subst.Act <- dnds$N.Act + dnds$S.Act
  dnds$Subst.Exp <- dnds$N.Exp + dnds$S.Exp
  dnds$Cov.Act <- dnds$Reads.Act/dnds$PopSize.Act
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
  
  
  
  # Do not remove window-codonsites where there is no dnds and no dn-ds information because of insufficient window sequences.
  # But do remove window codon sites in which there are no true dn/ds because of zero synonymous substitutions.
  dnds <- dnds[!is.na(dnds$dNdS.Exp), ]
  dim(dnds)
  summary(dnds)
  
  
  # Add a PSEUDOCOUNT so that dN/dS == 0 does not cause numerical instability
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
  
  # Remove samples where response, features is NA
  cleandnds_diversify <- na.omit(dnds[, c(feats, respname)])
  print(dim(cleandnds_diversify))
  print(summary(cleandnds_diversify))
  
  # Override default fit function from rfFuncs to use 501 trees instead of default 500 
  # so that we can break ties
  new_rfFuncs <- rfFuncs
  new_rfFuncs$fit <- 
    function (x, y, first, last, ...) 
    {
      library(randomForest)
      # keep track of which samples are in which bag in which trees
      randomForest(x, y, importance = first, ntree=501, keep.inbag=TRUE, ...)
    }
  # Random Forest Selection Function
  # This will auto resample and create bags for us.  
  # Only rank features on first iteration when all features are used, 
  #   so that we get more accurate depiction of how much
  #   feature matters when all features considered.
  control <- rfeControl(functions=new_rfFuncs, method="boot", number=FOLDS, 
                        #rerank=TRUE,  # rerank features after eliminate
                        saveDetails=TRUE,   # save predictions and variable importances from selection process
                        returnResamp="all",   # save all resampling summary metrics
                        verbose=TRUE
  )
  
 
  # This is automatically parallelized across each bag when you use doMC or doMPI library.
  # This takes 6 minutes even for A bag only has 5000 samples.  You can get timing through results$timing
  rfe_class_results <- rfe(x=cleandnds_diversify[, feats], y=cleandnds_diversify[, respname], 
                           sizes=c(1:length(feats)),
                           rfeControl=control, 
                           metric="Accuracy"  # This is for classification.  use another metric for regression
  )
  
  
  
  # save the classifier feature selection results (including its fit) to file 
  # so that we can load the environment variable back again later.
  save(rfe_class_results, file="rfe_class_results.RData")
  
  # save the training data
  save(cleandnds_diversify, file="cleandnds_diversify.RData")
  
  print(rfe_class_results)
  # list the chosen features
  print(predictors(rfe_class_results))  # results$optVariables also does the same)
  
  # plot the results
  fig <- ggplot(rfe_class_results, metric = rfe_class_results$metric[1], output = "layered")
  ggsave(filename="RandomForestClassifyDiversifyFeatureElbow.pdf", plot=fig, device=pdf)
  
  
  # print confusion matrix
  print(rfe_class_results$fit$confusion)
  
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
  
  # Remove samples where response, features is NA
  cleandnds <- na.omit(dnds[, c(respname, feats)])
  print("Cleaned dnds dimensions:")
  print(dim(cleandnds))
  print("Cleaned dnds summary:")
  print(summary(cleandnds))
  
  # Override default number of trees from 500 to 501 to break ties
  new_rfFuncs <- rfFuncs
  new_rfFuncs$fit <- 
    function (x, y, first, last, ...) 
    {
      library(randomForest)
      randomForest(x, y, importance = first, ntree=501, ...)
    }
  # Random Forest Selection Function
  # This will auto resample and create bags for us.
  control <- rfeControl(functions=new_rfFuncs, method="boot", number=FOLDS, 
                        #rerank=TRUE,  # rerank features after eliminate
                        saveDetails=TRUE,   # save predictions and variable importances from selection process
                        returnResamp="all",   # save all resampling summary metrics
                        verbose=TRUE
                        )
    
  # Automatically parallelized when you use library doMC or doMPI
  # This takes 6 minutes even for A bag only has 5000 samples.  You can get timing through results$timing
  rfe_cont_results <- rfe(x=cleandnds[, feats], y=cleandnds[, respname], 
                          sizes=c(1:length(feats)),
                          rfeControl=control, 
                          metric="RMSE"  # This is for continuous response
                          )
  
  # Save the rfe_cont_results environment object to file.
  save(rfe_cont_results, file="rfe_cont_results.RData")
  # Save the training dataset to feil
  save(cleandnds, file="cleandnds.RData")
  
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


# Does all the work for regression
do_predict_cont <- function() {
  
  dnds <- get_all_sim_dnds()
  dim(dnds)
  summary(dnds)
  head(dnds)
  object_size(dnds)
  mem_used()
  
  NUM_RESP_NAMES <- c("LOD_dNdS", "Dist_dn_minus_dS", "AbsLOD_dNdS", "AbsDist_dn_minus_dS")
  CAT_RESP_NAMES <- c("CrapLOD", "CrapDist", "wrongSelect")
  COVAR_NAMES <- colnames(dnds[sapply(dnds,is.numeric)])[!colnames(dnds[sapply(dnds,is.numeric)]) %in% NUM_RESP_NAMES]
  CAT_COVAR_NAMES <-  c() #c("IsLowSubst.Act")
  LM_COVAR_NAMES <- c(CAT_COVAR_NAMES, 
                      COVAR_NAMES[!(COVAR_NAMES %in% c("dNdS.Act", "dN_minus_dS.Act", "dN_minus_dS.Exp", 
                                                       # In separate analysis, conservation and entropy are highly correlated.
                                                       # When we use speedglm, it bugs out after it removes highly correlated variables.
                                                       # So we do it for them.
                                                       "ConserveTrueBase.Act", "ConserveTrueBase.Exp", "Window_Conserve.Act",                                                     
                                                       # These are highly correlated with N, S
                                                       #"Subst.Act", "Subst.Exp",
                                                       "EN.Exp", "ES.Exp", "EN.Act", "ES.Act",
                                                       "Window_Start", "Window_End", "CodonSite",
                                                       "N.Exp", "S.Exp", "dNdS.Exp",
                                                       "EntropyTrueBase.Exp"
                                                       #"Reads.Act"
                      )
                      )])
  
  feats <- c(LM_COVAR_NAMES)
  
  
  print("About to do random forest regression")
  rfe_cont_results <- rf_feat_sel_cont_rfe(dnds=dnds, respname="LOD_dNdS", feats=feats)
  
  # Get the predictions for all of the simulation data  
  lod_dnds_dat <- dnds[rowSums(is.na(dnds[, c("LOD_dNdS", feats)])) == 0,]  
  lod_dnds_dat$pred <- predict(rfe_cont_results$fit, lod_dnds_dat[, c(feats)])
  lod_dnds_dat$residual <- lod_dnds_dat$LOD_dNdS - lod_dnds_dat$pred
  
  # Get the MSE for all of the simulation data predictions
  mse <- mean((lod_dnds_dat$residual)^2)
  print("MSE for all simulation data")
  print(mse)
  
  # Get the RSquared for all of the simulation data predictions
  r2 <- rSquared(y=lod_dnds_dat$LOD_dNdS, resid=lod_dnds_dat$residual)
  print("RSquared for all simulation data")
  print(r2)
  
  
  # Save the predictions to file
  write.table(lod_dnds_dat, file="umberjack_accuracy_predict.csv", sep=",", row.names=FALSE)
  
  # Plot the random forest regression fit
  fig <- ggplot(lod_dnds_dat, aes(x=LOD_dNdS, y=pred)) + 
    geom_point(alpha=0.5, shape=1) +
    geom_smooth(method="lm") + 
    geom_abline(color="red") +
    xlab("\n Ln (Umberjack / True dNdS)") + 
    ylab("RF Predicted Ln (Umberjack / True dNdS) \n")
  ggtitle(paste("RandomForest Regression in R r^2=", r2, sep=""))
  
  ggsave(filename="RandomForestRegressionRsq.pdf", plot=fig, device=pdf)
  
  
  print(paste0("memused = ", mem_used()))
}



# Does all the work for finding Umberjack accuracy for classifying sites as Diversifying
do_predict_class_diversify <- function() {
  
  dnds <- get_all_sim_dnds()
  dim(dnds)
  summary(dnds)
  head(dnds)
  object_size(dnds)
  mem_used()
  
  NUM_RESP_NAMES <- c("LOD_dNdS", "Dist_dn_minus_dS", "AbsLOD_dNdS", "AbsDist_dn_minus_dS")
  CAT_RESP_NAMES <- c("CrapLOD", "CrapDist", "wrongSelect")
  COVAR_NAMES <- colnames(dnds[sapply(dnds,is.numeric)])[!colnames(dnds[sapply(dnds,is.numeric)]) %in% NUM_RESP_NAMES]
  CAT_COVAR_NAMES <-  c() #c("IsLowSubst.Act")
  LM_COVAR_NAMES <- c(CAT_COVAR_NAMES, 
                      COVAR_NAMES[!(COVAR_NAMES %in% c("dNdS.Act", "dN_minus_dS.Act", "dN_minus_dS.Exp",                                                        
                                                       "ConserveTrueBase.Act", "ConserveTrueBase.Exp", "Window_Conserve.Act",
                                                       "EN.Exp", "ES.Exp", "EN.Act", "ES.Act",
                                                       "Window_Start", "Window_End", "CodonSite",
                                                       "N.Exp", "S.Exp", "dNdS.Exp",
                                                       "EntropyTrueBase.Exp"
                                                       
                      )
                      )])
  
  feats <- c(LM_COVAR_NAMES)
  
  
  print("About to do random forest feature selection to determine what affects accuracy of umberjack predictions of diversifying sites")
  rfe_class_results <- rf_feat_sel_class_rfe(dnds=dnds, respname="wrongSelect", feats=feats)
  
  # Get the predictions for all of the simulation data
  wrongselect_dnds_dat <- dnds[rowSums(is.na(dnds[, c("wrongSelect", feats)])) == 0, ]  
  summary(wrongselect_dnds_dat)
  wrongselect_dnds_dat$pred <- predict(rfe_class_results$fit, wrongselect_dnds_dat[, c(feats)])
  
  # Make confusion matrix
  confuse <- with(wrongselect_dnds_dat, table(wrongSelect, pred))
  print(confuse)
  
  # Get accuracy for all of the simulation data predictions  (biased - should use OOB instead)
  accuracy <- sum(wrongselect_dnds_dat$wrongSelect == wrongselect_dnds_dat$pred, na.rm=TRUE)/sum(!is.na(wrongselect_dnds_dat$wrongSelect) & !is.na(wrongselect_dnds_dat$pred))
  print(accuracy)
  
    
  # Save the predictions to file
  write.table(wrongselect_dnds_dat, file="umberjack_diversify_accuracy_predict.csv", sep=",", row.names=FALSE)
  
  
  print(paste0("memused = ", mem_used()))
}