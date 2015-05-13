# Uses random forest feature selection to determine best features that affect Umberjack accuracy in classifying
# sites as diversifying or purifying.
# This is just an mpi wrapper for load_all_sim_dnds.R which does all the work


library(caret)
library(plyr)
library(doMPI)  # for MPI parallelism
library(randomForest)
library(pryr)


FOLDS <- 5
SEED <- 7
PROCS <- 5

source("./load_all_sim_dnds.R")

slaves <- startMPIcluster(count=PROCS)
registerDoMPI(slaves)

do_predict_class_diversify()

closeCluster(slaves)
mpi.quit()