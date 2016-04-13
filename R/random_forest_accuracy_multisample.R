# Uses random forest feature selection to determine best features that affect dnds estimate accuracy
# This is just an mpi wrapper for load_all_sim_dnds.R which does all the work

# Execute with optional parameters  Rscript random_forest_accuracy_multisample.R [per sample dnds csv] [cores per random forest]

library(caret)
library(plyr)
library(doMPI)  # for MPI parallelism
library(randomForest)
library(pryr)
#install.packages("pryr", repos=c("http://cran.stat.sfu.ca/"))


FOLDS <- 5
#SEED <- 7
CORES_PER_RF <- 4  # We do cross validation in parallel.  But for each cross validation, we also do the random forest trees in parallel
TREES_PER_RF <- 501  # total random forest trees to execute for each cross validation

source("./load_all_sim_dnds.R")

args <- commandArgs(TRUE)
if (length(args) >= 1) {
  dnds_filename <- args[1]
} else {
  dnds_filename <- NULL
}

if (length(args) >= 2) {
  cores_per_rf <- as.integer(args[2])
} else {
  cores_per_rf <- CORES_PER_RF
}
if (length(args) >= 3) {
  folds <- as.integer(args[3])
} else {
  folds <- FOLDS
}
if (length(args) >= 4) {
  seed <- as.integer(args[4])
} else {
  seed <- NULL
}
slaves <- startMPIcluster(count=(folds * cores_per_rf))
registerDoMPI(slaves)

do_predict_cont(dnds_filename=dnds_filename, folds=folds, trees_per_rf=TREES_PER_RF, cores_per_rf=cores_per_rf, seed=seed)

closeCluster(slaves)
mpi.quit()