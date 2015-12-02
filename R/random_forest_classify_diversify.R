# Uses random forest feature selection to determine best features that affect Umberjack accuracy in classifying
# sites as diversifying or purifying.
# This is just an mpi wrapper for load_all_sim_dnds.R which does all the work



library(plyr)
library(doMPI)  # for MPI parallelism

library(pryr)


FOLDS <- 5
SEED <- 7
PROCS <- 5

source("./load_all_sim_dnds.R")

TRAIN_DNDS_CSV <- NULL
args <- commandArgs(TRUE)
if (length(args) >= 1) {
  TRAIN_DNDS_CSV <- args[1]  
}

print (paste0("Input argument=", TRAIN_DNDS_CSV))

slaves <- startMPIcluster(count=PROCS, verbose=TRUE)
registerDoMPI(slaves)

do_predict_class_diversify(train_dnds_csv=TRAIN_DNDS_CSV, xfold=FOLDS)

closeCluster(slaves)
mpi.quit()