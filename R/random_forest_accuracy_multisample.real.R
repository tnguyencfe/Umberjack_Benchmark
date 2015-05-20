# Uses random forest feature selection to determine best features that affect dnds estimate accuracy
# This is just an mpi wrapper for load_all_sim_dnds.R which does all the work



library(caret)
library(plyr)
library(doMPI)  # for MPI parallelism
library(randomForest)
library(pryr)
#install.packages("pryr", repos=c("http://cran.stat.sfu.ca/"))
library(miscTools)  # for rSquared

FOLDS <- 5
SEED <- 7
PROCS <- 5

source("./load_all_sim_dnds.R")

slaves <- startMPIcluster(count=PROCS)
registerDoMPI(slaves)

do_predict_cont_real()

closeCluster(slaves)
mpi.quit()