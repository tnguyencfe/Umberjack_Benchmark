library(plyr)
library(pryr)
library(caret)
library(randomForest)
library(foreach)

PSEUDOCOUNT <- 1e-7



NUM_RESP_NAMES <- c("LOD_dNdS", "Dist_dn_minus_dS", "AbsLOD_dNdS", "AbsDist_dn_minus_dS", "SqDist_dn_minus_dS")
CAT_RESP_NAMES <- c("CrapLOD", "CrapDist")

# All numeric varaibles
NUM_NAMES <- c("Window_Start",
               "Window_End",
               "CodonSite",
               "Is_Break",
               "BreakRatio.Act",
               "Reads.Act",
               "UnambigCodonRate.Act",
               "AADepth.Act",
               "PopSize.Act",
               "ConserveCodon.Act",
               "EntropyCodon.Act",
               "UnknownPerCodon.Act",
               "ErrPerCodon.Act",
               "N.Act",
               "S.Act",
               "EN.Act",
               "ES.Act",
               "dNdS.Act",
               "dN_minus_dS.Act",
               "TreeLen.Act",
               "TreeLenPerRead.Act",
               "TreeDepth.Act",
               "TreeDist.Act",
               "TreeDistPerRead.Act",
               "Polytomy.Act",
               "PolytomyPerRead.Act",
               "P_SameCodonFreq.Act",
               "ResolvedPerSub.Act",
               "ConserveCodon.Exp",
               "EntropyCodon.Exp",
               "N.Exp",
               "S.Exp",
               "EN.Exp",
               "ES.Exp",
               "dNdS.Exp",
               "dN_minus_dS.Exp")

# Numeric variables that might affect Umberjack accuracy
COVAR_NAMES <- NUM_NAMES[!NUM_NAMES %in% 
                           c(NUM_RESP_NAMES,
                             "dNdS.Act", "dNdS.Exp", "dN_minus_dS.Act", "dN_minus_dS.Exp", 
                             # In separate analysis, conservation and entropy are highly correlated.
                             # When we use speedglm, it bugs out after it removes highly correlated variables.
                             # So we do it for them.
                             "ConserveCodon.Act", "ConserveCodon.Exp", "Window_Conserve.Act",
                             #"UnambigCodonRate.Act", 
                             #"Window_UnambigCodonRate.Act",
                             # These are highly correlated with N, S
                             "Subst.Act", "Subst.Exp",
                             "EN.Exp", "ES.Exp", "EN.Act", "ES.Act",  "N.Exp", "S.Exp",
                             "Window_Start", "Window_End", "CodonSite", "Reads.Act", "PopSize.Act", "Is_Break"
                           )
                         ]

# These variables apply to the entire window not just a window-codon site
WINDOW_COVAR_NAMES <- c("BreakRatio.Act", "Window_Breaks", "TreeLen.Act", "TreeDepth.Act", "TreeDist.Act", "TreeDistPerRead.Act", 
                        "Cov.Act", "Polytomy.Act", "PolytomyPerRead.Act", "TreeLenPerRead.Act", "Reads.Act", "WinP_SameCodonFreq.Act",
                        "Window_Entropy.Act", "Window_UnambigCodonRate.Act", "Window_ErrPerCodon.Act", "Window_Subst.Act")

# These variables apply only to specific window-codon site
WINDOW_SITE_COVAR_NAMES <- COVAR_NAMES[!COVAR_NAMES %in% WINDOW_COVAR_NAMES]

# categorical variables
CAT_COVAR_NAMES <- c("IsLowSubst.Act")

# variables used in linear regression
LM_COVAR_NAMES <- c(CAT_COVAR_NAMES, COVAR_NAMES)

REAL_LM_COVAR_NAMES <- c("Reads.Act",
                         "UnambigCodonRate.Act",
                         "UnknownPerCodon.Act",
                         "AADepth.Act",
                         "ConserveCodon.Act",
                         "EntropyCodon.Act",
                         "N.Act",
                         "S.Act",
                         "TreeLen.Act",
                         "TreeDepth.Act",
                         "Polytomy.Act",
                         "ResolvedPerSub.Act",                                                  
                         "Window_Entropy.Act", "Window_UnambigCodonRate.Act", "Window_Subst.Act"
                         )



nice <- function(name) {
  if (name == "IsLowSubst.Act") {
    return ("Only Ambig Window Phylogney Subs")
  } else if (name == "BreakRatio.Act") {
    return ("Metric of window breakpoints & \ncloseness to centre of window")
  } else if (name == "Window_Breaks") {
    return ("Recombination Breaks in Window")
  } else if (name == "Genome_Breaks") {
    return ("Recombination Breaks in Genome")
  } else if (name == "UnambigCodonRate.Act") {
    return ("Window-Site Unambiguous Codons Per Read")
  } else if (name == "AADepth.Act") {
    return ("Window-Site Unambiguous AA Per Read")
  } else if (name == "EntropyCodon.Act") {
    return ("Window-Site Codon Entropy")
  } else if (name == "UnknownPerCodon.Act") {
    return ("Window-Site Unknown Bases Per Codon")
  } else if (name == "ErrPerCodon.Act") {
    return ("Window-Site Erroneously Sequenced Bases Per Codon")
  } else if (name == "N.Act") {
    return ("Window-Site Nonsyn subs")
  } else if (name == "S.Act") {
    return ("Window-Site Syn subs")
  } else if (name == "TreeLen.Act") {
    return ("Window Tree Length (Subs/Site)")
  } else if (name == "TreeDepth.Act") {
    return ("Window Tree Depth (Subs/Site)")
  } else if (name == "TreeDistPerRead.Act") {
    return ("Ave Robinson Foulds (Window Tree, Expected Tree) \n Weighted By Genome Width Covered By Expected Tree\nNormalized by Total Reads")
  } else if (name == "TreeDist.Act") {
    return ("Ave Robinson Foulds (Window Tree, Expected Tree) \n Weighted By Genome Width Covered By Expected Tree")
  } else if (name == "EntropyCodon.Exp") {
    return ("True Population Site Codon Entropy")
  } else if (name == "WinP_SameCodonFreq.Act") {
    return ("Window Ave log10 P(window codon distro matches full popn)")
  } else if (name == "P_SameCodonFreq.Act") {
    return ("log10 P(window codon distro matches full popn)")
  } else if (name == "ResolvedPerSub.Act ") {
    return ("Fraction of Substitutions from Ambiguous Codons")
  }
  else if (name == "N.Exp") {
    return ("True Population Site Nonsyn Subs")
  } else if (name == "S.Exp") {
    return ("True Population Site Syn Subs")
  } else if (name == "Cov.Act") {
    return ("Window Read Coverage Per Individual")
  } else if (name == "Window_Entropy.Act") {
    return ("Codon Entropy Ave Across Window")
  } else if (name == "Window_UnambigCodonRate.Act") {
    return ("Unambig Codons Per Read Ave Across Window")
  } else if (name == "Window_ErrPerCodon.Act") {
    return ("Erroneously Sequenced Bases Per Codon Ave Across Window")
  } else if (name == "Window_Subst.Act") {
    return ("Phylogeny-Site Substitutions Ave Across Window")
  } else if (name == "LOD_dNdS") {
    return ("log2(inferred window site dn/ds) - log2(expected site dn/ds)")
  } else if (name == "Dist_dn_minus_dS") {
    return ("(inferred window site dn-ds) - (expected site dn-ds)")
  } else if (name == "AbsLOD_dNdS") {
    return ("|log2(inferred window site dn/ds) - log2(expected site dn/ds)|")
  } else if (name == "AbsDist_dn_minus_dS") {
    return ("|(inferred window site dn-ds) - (expected site dn-ds)|")
  } else if (name == "SqDist_dn_minus_dS") {
    return ("[(inferred window site dn-ds) - (expected site dn-ds)]^2")
  } else if (name == "WinAbsLOD_dNdS") {
    return ("|log2(inferred window site dn/ds) - log2(expected site dn/ds)| \n Ave Across Window")
  } else if (name == "WinAbsDist_dn_minus_dS") {
    return ("|(inferred window site dn-ds) - (expected site dn-ds)|\nAve Across Window")
  } else if (name == "WinSqDist_dn_minus_dS") {
    return ("[(inferred window site dn-ds) - (expected site dn-ds)] ^2\n Ave Across Window")
  } else {
    return (name)
  }
}
  
  
get_all_sim_dnds <- function(dnds_filename=NULL) {
   
  DNDS_FILENAME <- "../simulations/out/collate_all.treedist.csv"
    
  if (is.null(dnds_filename))
  {
    dnds_filename <- DNDS_FILENAME
  }
  
  print (paste0("Using Training Data ", dnds_filename))
  
  # Window_Start,Window_End,CodonSite,File,Is_Break,Reads.Act,UnambigCodonRate.Act,AADepth.Act,
  #PopSize.Act,ConserveCodon.Act,EntropyCodon.Act,UnknownPerCodon.Act,ErrPerCodon.Act,N.Act,S.Act,EN.Act,ES.Act,
  #dNdS.Act,dN_minus_dS.Act,TreeLen.Act,TreeDepth.Act,TreeDistPerRead.Act,ConserveCodon.Exp,EntropyCodon.Exp,
  #N.Exp,S.Exp,EN.Exp,ES.Exp,dNdS.Exp,dN_minus_dS.Exp
  dnds <- read.table(dnds_filename, header=TRUE, sep=",", na.strings=c("", "None"))
  dnds$Subst.Act <- dnds$N.Act + dnds$S.Act
  dnds$Subst.Exp <- dnds$N.Exp + dnds$S.Exp
  dnds$Cov.Act <- dnds$Reads.Act/dnds$PopSize.Act  
  dnds$IsLowSubst.Act <- as.factor((dnds$N.Act > 0 & dnds$N.Act < 1) | (dnds$S.Act > 0 & dnds$S.Act < 1))
  dnds$PolytomyPerRead.Act <- dnds$Polytomy.Act / dnds$Reads.Act
  dnds$TreeLenPerRead.Act <- dnds$TreeLen.Act / dnds$Reads.Act
  
  if (all(0 <= dnds$P_SameCodonFreq.Act, na.rm=TRUE) & all(dnds$P_SameCodonFreq.Act <= 1, na.rm=TRUE))  # we want log10 probabilities
  {
    dnds$P_SameCodonFreq.Act <- log10(dnds$P_SameCodonFreq.Act)
    dnds$P_SameCodonFreq.Act[is.infinite(dnds$P_SameCodonFreq.Act)] <- NA
  }
  
  if (!"TreeDist.Act"  %in% colnames(dnds)) {
    dnds$TreeDist.Act <- dnds$TreeDistPerRead.Act * dnds$Reads.Act
  }
  
  # Average across all codon sites in a window
  per_window_ave <- ddply(.data=dnds, .variables=c("File", "Window_Start"), 
                          .fun=function(x) {                            
                            data.frame(Window_Entropy.Act=mean(x$EntropyCodon.Act, na.rm=TRUE),
                                       Window_UnambigCodonRate.Act=mean(x$UnambigCodonRate.Act, na.rm=TRUE),
                                       Window_ErrPerCodon.Act=mean(x$ErrPerCodon.Act, na.rm=TRUE),
                                       Window_Subst.Act=mean(x$Subst.Act, na.rm=TRUE),
                                       WinP_SameCodonFreq.Act=mean(x$P_SameCodonFreq.Act, na.rm=TRUE)
                                       )
                          })
  
  
  dnds <- merge(x=dnds, y=per_window_ave, by=c("File", "Window_Start"), all=TRUE, sort=TRUE)
  
  # Check that if a codon site in a dataset is a breakpoint, it's reported as a breakpoint in all the windows for that dataset
  uniq_breaks_per_codon <- aggregate(Is_Break ~ File + CodonSite, data=dnds, 
                                FUN=function(x) {return (length(unique(x)))})
  if (sum(uniq_breaks_per_codon$Is_Break > 1) > 0) {
    print(head(uniq_breaks_per_codon[uniq_breaks_per_codon$Is_Break > 1,]))
    stop("Codon Site in dataset is reported as both breakpoint and non-breakpoint")
  }
  
  # Count breakpoints per window
  breaks_per_win <- aggregate(Is_Break ~ File + Window_Start, data=dnds, FUN=sum)
  colnames(breaks_per_win)[grep("Is_Break", colnames(breaks_per_win))] <- "Window_Breaks"
  
  dnds <- merge(x=dnds, y=breaks_per_win, all=TRUE, sort=TRUE)
  
  # Count breakpoints per genome
  breaks_per_site <- aggregate(Is_Break ~ File + CodonSite, data=dnds, FUN=function(x){unique(x)})
  breaks_per_genome <- aggregate(Is_Break ~ File, data=breaks_per_site, FUN=sum)
  colnames(breaks_per_genome)[grep("Is_Break", colnames(breaks_per_genome))] <- "Genome_Breaks"

  dnds <- merge(x=dnds, y=breaks_per_genome, all=TRUE, sort=TRUE)
  
#   # Do not remove window-codonsites where there is no dnds and no dn-ds information because of insufficient window sequences.
#   # But do remove window codon sites in which there are no true dn/ds because of zero synonymous substitutions.
#   dnds <- dnds[!is.na(dnds$dNdS.Exp), ]
#   dim(dnds)
#   summary(dnds)
  
  
  # Add a PSEUDOCOUNT so that dN/dS == 0 does not cause numerical instability
  dnds$LOD_dNdS <- log2(dnds$dNdS.Act + PSEUDOCOUNT) - log2(dnds$dNdS.Exp + PSEUDOCOUNT)
  dnds$AbsLOD_dNdS <- abs(dnds$LOD_dNdS)


  dnds$CrapLOD <- FALSE
  dnds$CrapLOD <- as.factor(abs(dnds$LOD_dNdS) >= 1)
  dnds$CrapLOD[is.na(dnds$dNdS.Act)] <- TRUE

  dnds$Dist_dn_minus_dS <- dnds$dN_minus_dS.Act - dnds$dN_minus_dS.Exp
  dnds$AbsDist_dn_minus_dS <- abs(dnds$Dist_dn_minus_dS)
  dnds$SqDist_dn_minus_dS <- dnds$Dist_dn_minus_dS ^ 2

  dnds$CrapDist <- FALSE
  dnds$CrapDist <- as.factor(abs(dnds$Dist_dn_minus_dS) >= 1)
  dnds$CrapDist[is.na(dnds$CrapDist)] <- TRUE
  
  dnds$wrongSelect <- FALSE
  dnds$wrongSelect <- as.factor((dnds$dN_minus_dS.Act < 0 & dnds$dN_minus_dS.Exp > 0) | (dnds$dN_minus_dS.Act > 0 & dnds$dN_minus_dS.Exp < 0))
  dnds$wrongSelect[is.na(dnds$wrongSelect)] <- TRUE
  
  return (dnds)
}


get_window_sim_dnds <- function(dnds) {
  # When we compare umberjack dnds against variables that affect entire windows (as opposed to window-codon sites),
  # we need to use umberjack dnds averaged across window to avoid excess noise when plotting
  
  if (is.null(dnds)) {
    dnds <- get_all_sim_dnds()
  }
  window_means <- aggregate(cbind(AbsLOD_dNdS, AbsDist_dn_minus_dS, SqDist_dn_minus_dS)  ~ File + Window_Start,
                      data=dnds,
                      FUN=function(y) {mean(y, na.rm=TRUE)})
  colnames(window_means)[grep("AbsLOD_dNdS", colnames(window_means))] <- "WinAbsLOD_dNdS"
  colnames(window_means)[grep("AbsDist_dn_minus_dS", colnames(window_means))] <- "WinAbsDist_dn_minus_dS"
  colnames(window_means)[grep("SqDist_dn_minus_dS", colnames(window_means))] <- "WinSqDist_dn_minus_dS"
  
  
  
  vars_formula <- as.formula(paste0(
    "cbind(",
    paste0(WINDOW_COVAR_NAMES, collapse=", "),
    ") ~ File + Window_Start"))
  window_vars <- aggregate(vars_formula, data=dnds, FUN=function(x) {unique(x)})
  
  window <- merge(x=window_means, y=window_vars, all=TRUE)
  
  return (window)
}

# Random Forest, Recursive Feature Elmination  (Backwards Selection)  for Classification
rf_feat_sel_class_rfe <- function(dnds, respname, feats, xfold=1) {
  
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
  control <- rfeControl(functions=new_rfFuncs, method="boot", number=xfold, 
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
rf_feat_sel_cont_rfe <- function(dnds, respname, feats, folds, trees_per_rf, cores_per_rf) {
  
  print(paste0("Response=", respname))
  print(paste0("Features=", paste0(feats, collapse=", ")))
  
  # Remove samples where response, features is NA
  cleandnds <- na.omit(dnds[, c(respname, feats)])
  print("Cleaned dnds dimensions:")
  print(dim(cleandnds))
  print("Cleaned dnds summary:")
  print(summary(cleandnds))
  
  # Override default number of trees from 500 to 501 to break ties, Use parallelized random forest
  new_rfFuncs <- rfFuncs
  new_rfFuncs$fit <- 
    function (x, y, first, last, ...) 
    {
      library(randomForest)      
      
      # Use randomForest's builtin combine function     
      # Divide trees evenly into the amount we want to parallelize.
      # If it won't divide evenly, then the last tree will have the leftovers.
      trees_per_core <- trees_per_rf %/% cores_per_rf
      leftover_trees_per_core <- trees_per_rf - trees_per_core * (cores_per_rf-1)
      
      parallel_randomForest <- foreach(ntree=c(rep(trees_per_core, cores_per_rf-1), leftover_trees_per_core),
                                       .combine=combine, .packages='randomForest') %dopar% {
        randomForest(x, y, importance = first, ntree=ntree, ...)
      }
      return (parallel_randomForest)
      
    } 
  
  # Random Forest Selection Function
  # This will auto resample and create bags for us.
  control <- rfeControl(functions=new_rfFuncs, method="boot", number=folds, 
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
  #save(rfe_cont_results, file="rfe_cont_results.RData")
  # Save the training dataset to file
  #save(cleandnds, file="cleandnds.RData")
  
  print(rfe_cont_results)
  # list the chosen features
  print(paste0("Predictors\n", predictors(rfe_cont_results)))  # results$optVariables also does the same
  
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
do_predict_cont <- function(dnds_filename=NULL, respname, folds=5, trees_per_rf=501, cores_per_rf=1) {
  
  dnds <- get_all_sim_dnds(dnds_filename)
  dim(dnds)
  summary(dnds)
  head(dnds)
  object_size(dnds)
  print(paste0("mem used from dnds=", mem_used()))
  
  feats <- c(LM_COVAR_NAMES)
  
  
  print("About to do random forest regression")
  rfe_cont_results <- rf_feat_sel_cont_rfe(dnds=dnds, respname=respname, feats=feats, 
                                           folds=folds, trees_per_rf=trees_per_rf, cores_per_rf=cores_per_rf)
  
  print(paste0("Mem Bytes after RF=", mem_used()))
  
  # Get the predictions for all of the simulation data  
  lod_dnds_dat <- dnds[rowSums(is.na(dnds[, c("SqDist_dn_minus_dS", feats)])) == 0,]  
  lod_dnds_dat$pred <- predict(rfe_cont_results$fit, lod_dnds_dat[, c(feats)])
  lod_dnds_dat$residual <- lod_dnds_dat$SqDist_dn_minus_dS - lod_dnds_dat$pred
  
  # Get the MSE for all of the simulation data predictions
  mse <- mean((lod_dnds_dat$residual)^2)
  print("MSE for all simulation data")
  print(mse)
  
  # Get the RSquared for all of the simulation data predictions
  r2 <- rSquared(y=lod_dnds_dat$SqDist_dn_minus_dS, resid=lod_dnds_dat$residual)
  print("RSquared for all simulation data")
  print(r2)
  
  
  # Save the predictions to file
  write.table(lod_dnds_dat, file="umberjack_accuracy_predict.csv", sep=",", row.names=FALSE)
  
  # Plot the random forest regression fit
  fig <- ggplot(lod_dnds_dat, aes(x=SqDist_dn_minus_dS, y=pred)) + 
    geom_point(alpha=0.5, shape=1) +
    geom_smooth(method="lm") + 
    geom_abline(color="red") +
    xlab("\n [(Umberjack dn-ds) - (True dn-ds)]^2") + 
    ylab("RF Predicted [(Umberjack dn-ds) - (True dn-ds)]^2 \n")
  ggtitle(paste("RandomForest Regression in R r^2=", r2, sep=""))
  
  ggsave(filename="RandomForestRegressionRsq.pdf", plot=fig, device=pdf)
  
  
  print(paste0("memused after finish= ", mem_used()))
}



# Does all the work for finding Umberjack accuracy for classifying sites as Diversifying
do_predict_class_diversify <- function(train_dnds_csv=NULL, xfold=1) {
  
  dnds <- get_all_sim_dnds(dnds_filename=train_dnds_csv)
  dim(dnds)
  summary(dnds)
  head(dnds)
  object_size(dnds)
  print(paste0("dnds mem=", mem_used()))
    
  feats <- c(LM_COVAR_NAMES)
  
  
  print("About to do random forest feature selection to determine what affects accuracy of umberjack predictions of diversifying sites")
  rfe_class_results <- rf_feat_sel_class_rfe(dnds=dnds, respname="wrongSelect", feats=feats, xfold=xfold)
  
  # Get the predictions for all of the simulation data
  wrongselect_dnds_dat <- dnds[rowSums(is.na(dnds[, c("wrongSelect", feats)])) == 0, ]  
  summary(wrongselect_dnds_dat)
  wrongselect_dnds_dat$pred <- predict(rfe_class_results$fit, wrongselect_dnds_dat[, c(feats)])
  
  # Make confusion matrix
  confuse <- with(wrongselect_dnds_dat, table(wrongSelect, pred))
  print(confuse)
  
  # save env to file
  save(rfe_class_results, file="rfe_class_results.RData")
  
  # Get accuracy for all of the simulation data predictions  (biased - should use OOB instead)
  accuracy <- sum(wrongselect_dnds_dat$wrongSelect == wrongselect_dnds_dat$pred, na.rm=TRUE)/sum(!is.na(wrongselect_dnds_dat$wrongSelect) & !is.na(wrongselect_dnds_dat$pred))
  print(accuracy)
  
  # Save the predictions to file
  write.table(wrongselect_dnds_dat, file="umberjack_diversify_accuracy_predict.csv", sep=",", row.names=FALSE)
  
  print(paste0("memused = ", mem_used()))
  
}



# Does all the work for finding Umberjack accuracy for classifying real sites 
# with unknown expected values as Diversifying
do_predict_class_diversify_real <- function() {
  
  dnds <- get_all_sim_dnds()
  dim(dnds)
  summary(dnds)
  head(dnds)
  object_size(dnds)
  print(paste0("dnds mem=", mem_used()))
  
  NUM_RESP_NAMES <- c("LOD_dNdS", "Dist_dn_minus_dS", "AbsLOD_dNdS", "AbsDist_dn_minus_dS")
  CAT_RESP_NAMES <- c("CrapLOD", "CrapDist", "wrongSelect")
  COVAR_NAMES <- colnames(dnds[sapply(dnds,is.numeric)])[!colnames(dnds[sapply(dnds,is.numeric)]) %in% NUM_RESP_NAMES]
  CAT_COVAR_NAMES <-  c() #c("IsLowSubst.Act")
  LM_COVAR_NAMES <- c(CAT_COVAR_NAMES, 
                      COVAR_NAMES[!(COVAR_NAMES %in% c("dNdS.Act", "dNdS.Exp",
                                                       "dN_minus_dS.Act", "dN_minus_dS.Exp",                                                        
                                                       "ConserveTrueBase.Act", "ConserveTrueBase.Exp", 
                                                       "Window_Conserve.Act",
                                                       "Window_ErrPerCodon.Act",
                                                       "EN.Exp", "ES.Exp", "EN.Act", "ES.Act",
                                                       "Window_Start", "Window_End", "CodonSite",
                                                       "PopSize.Act",
                                                       "Cov.Act",
                                                       "ErrPerCodon.Act",
                                                       "N.Exp", "S.Exp", 
                                                       "EntropyTrueBase.Exp",
                                                       "Subst.Act", "Subst.Exp"
                                                       
                      )
                      )])
  
  feats <- c(LM_COVAR_NAMES)
  
  
  print("About to do random forest feature selection to determine what affects accuracy of umberjack predictions of diversifying real sites")
  rfe_class_results_real <- rf_feat_sel_class_rfe(dnds=dnds, respname="wrongSelect", feats=feats)
  
  # save env to file
  save(rfe_class_results_real, file="rfe_class_results_real.RData")
  
  # Get the predictions for all of the simulation data
  wrongselect_dnds_dat <- dnds[rowSums(is.na(dnds[, c("wrongSelect", feats)])) == 0, ]  
  summary(wrongselect_dnds_dat)
  wrongselect_dnds_dat$pred <- predict(rfe_class_results_real$fit, wrongselect_dnds_dat[, c(feats)])
  
  # Make confusion matrix
  confuse <- with(wrongselect_dnds_dat, table(wrongSelect, pred))
  print(confuse)
  
  # Get accuracy for all of the simulation data predictions  (biased - should use OOB instead)
  accuracy <- sum(wrongselect_dnds_dat$wrongSelect == wrongselect_dnds_dat$pred, na.rm=TRUE)/sum(!is.na(wrongselect_dnds_dat$wrongSelect) & !is.na(wrongselect_dnds_dat$pred))
  print(accuracy)
  
  # Save the predictions to file
  write.table(wrongselect_dnds_dat, file="umberjack_diversify_accuracy_predict_real.csv", sep=",", row.names=FALSE)
  
  print(paste0("memused = ", mem_used()))
  
}

# Does all the work for finding Umberjack accuracy for regression on real sites 
do_predict_cont_real <- function(dnds_filename=NULL, folds=5, trees_per_rf=501, cores_per_rf=1) {
  
  dnds <- get_all_sim_dnds(dnds_filename)
  dim(dnds)
  summary(dnds)
  head(dnds)
  object_size(dnds)
  print(paste0("mem used from dnds=", mem_used()))
  
  feats <- c(REAL_LM_COVAR_NAMES)
  
  
  print("About to do random forest regression on features available for realistic data")
  rfe_cont_results_real <- rf_feat_sel_cont_rfe(dnds=dnds, respname="SqDist_dn_minus_dS", feats=feats, 
                                                folds=folds, trees_per_rf=trees_per_rf, cores_per_rf=cores_per_rf)
  
  # Save environment var to file
  save(rfe_cont_results_real, file="rfe_cont_results_real.RData")
    
  # plot the results
  fig <- ggplot(rfe_cont_results_real, metric = rfe_cont_results_real$metric[1], output = "layered")
  ggsave(filename="RandomForestRegressionFeatureElbow.real.pdf", plot=fig, device=pdf)
    
  print(paste0("Mem Bytes after RF=", mem_used()))
  
  # Get the predictions for all of the simulation data  
  lod_dnds_dat <- dnds[rowSums(is.na(dnds[, c("SqDist_dn_minus_dS", feats)])) == 0,]  
  lod_dnds_dat$pred <- predict(rfe_cont_results_real$fit, lod_dnds_dat[, c(feats)])
  lod_dnds_dat$residual <- lod_dnds_dat$SqDist_dn_minus_dS - lod_dnds_dat$pred
  
  # Get the MSE for all of the simulation data predictions
  mse <- mean((lod_dnds_dat$residual)^2)
  print("MSE for all simulation data")
  print(mse)
  
  # Get the RSquared for all of the simulation data predictions
  r2 <- rSquared(y=lod_dnds_dat$SqDist_dn_minus_dS, resid=lod_dnds_dat$residual)
  print("RSquared for all simulation data")
  print(r2)
  
  
  # Save the predictions to file
  write.table(lod_dnds_dat, file="umberjack_accuracy_predict.real.csv", sep=",", row.names=FALSE)
  
  # Plot the random forest regression fit
  fig <- ggplot(lod_dnds_dat, aes(x=SqDist_dn_minus_dS, y=pred)) + 
    geom_point(alpha=0.5, shape=1) +
    geom_smooth(method="lm") + 
    geom_abline(color="red") +
    xlab("\n [(Umberjack dn-ds) - (True dn-ds)]^2") + 
    ylab("RF Predicted [(Umberjack dn-ds) - (True dn-ds)]^2 \n")
  ggtitle(paste("RandomForest Regression in R r^2=", r2, sep=""))
  
  ggsave(filename="RandomForestRegressionRsq.real.pdf", plot=fig, device=pdf)
  
  
  print(paste0("memused after finish= ", mem_used()))
}
