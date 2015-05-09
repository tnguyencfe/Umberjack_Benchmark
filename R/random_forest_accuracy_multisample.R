# Uses random forest feature selection to determine best features that affect dnds estimate accuracy

library(knitr)
library(caret)
library(plyr)
#library(doMC)  # for parallele rfe
library(doMPI)  # for MPI parallelism
library(randomForest)
library(pryr)


FOLDS <- 1
SEED <- 7
PROCS <- 1

slaves <- startMPIcluster(count=PROCS)
registerDoMPI(slaves)


## Note: if the underlying model also uses foreach, the
## number of cores specified above will double (along with
## the memory requirements)
#registerDoMC(cores = PROCS)


source("./load_all_sim_dnds.R")


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


# crapLOD_results <- rf_feat_sel_class_rfe(dnds=dnds[!is.na(dnds$CrapLOD) & !is.na(dnds$dNdS.Act),], respname="CrapLOD", feats=feats)
# wrongSelect_results <- rf_feat_sel_class_rfe(dnds=dnds[!is.na(dnds$wrongSelect) & !is.na(dnds$dNdS.Act),], respname="wrongSelect", feats=feats)
# 
# 
# rf_fit <- rf_train(dnds=dnds, respname="wrongSelect", feats=feats)
# 
# 
# 
# longshot_dnds <- read.table("/home/thuy/gitrepo/MutationPatterns/out_maskstopcodon_remdup/140415_M01841_0059_000000000-A64EA/collate_all.longshot.csv", sep=",", header=TRUE)
# summary(longshot_dnds)
# head(longshot_dnds)
# 
# 
# preds <- predict(wrongSelect_results$fit, dnds[!is.na(dnds$CrapLOD) & !is.na(dnds$dNdS.Act),])
# 
# preds_longshot <- predict(wrongSelect_results$fit, longshot_dnds[!is.na(longshot_dnds$dNdS.Act), 
#                                                     c("Window_Start", "File", "CodonSite", "S.Act", "N.Act", "EntropyTrueBase.Act", "UnambigCodonRate.Act", "Reads.Act")])
# sum(preds_longshot==TRUE)
# 
# length(preds_longshot)
# 
# 
# longshot_pred_dat <- longshot_dnds[!is.na(longshot_dnds$dNdS.Act), 
#                                    c("Window_Start", "File", "CodonSite", "S.Act", "N.Act", "EntropyTrueBase.Act", "UnambigCodonRate.Act", "Reads.Act")]
# 
# longshot_pred_dat$Pred <- preds_longshot
# 
# write.table(longshot_pred_dat, "longshot_pred.csv", sep=",", row.names=FALSE)
# 
# 
# head(longshot_pred_dat[longshot_pred_dat$Pred==TRUE,])


print("About to do random forest regression")
rf_feat_sel_cont_rfe <- rf_feat_sel_cont_rfe(dnds=dnds, respname="LOD_dNdS", feats=feats)
# Get the predictions for the training data?
preds <- predict(rfe_cont_results$fit, dnds[!is.na(dnds$LOD_dNdS),])
# Get the mean squared error for each tree
# mse = (sum of squared residuals)/n
rfe_cont_results$fit$mse
# r-sq for each tree = 1-mse/var(response)
rfe_cont_results$fit$mse


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


closeCluster(slaves)
mpi.quit()