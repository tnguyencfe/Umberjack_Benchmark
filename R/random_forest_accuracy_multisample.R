# Uses random forest feature selection to determine best features that affect dnds estimate accuracy
# This is just an mpi wrapper for load_all_sim_dnds.R which does all the work


library(knitr)
library(caret)
library(plyr)
#library(doMC)  # for parallele rfe
library(doMPI)  # for MPI parallelism
library(randomForest)
library(pryr)
library(miscTools)  # for rSquared

FOLDS <- 5
SEED <- 7
PROCS <- 5

source("./load_all_sim_dnds.R")

slaves <- startMPIcluster(count=PROCS)
registerDoMPI(slaves)

do_predict_cont()

closeCluster(slaves)
mpi.quit()