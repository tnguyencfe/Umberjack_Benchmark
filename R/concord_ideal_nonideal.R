
source("plot_helper.R")

READLEN <- 250
BASE_SIZE <- 24
THESIS_OUTDIR <- "/home/thuy/gitrepo/MutationPatterns/tex"
  
# Find the concordance of simulated datasets, stratify by "ideal" vs "nonideal" scenario
concord <- read.table("/home/thuy/gitrepo/Umberjack_Benchmark/R/training/dataset_concord.csv", sep=",", header=TRUE)
head(concord)
summary(concord)

# simulation configs
sims <- read.table('/home/thuy/gitrepo/Umberjack_Benchmark/sim_config/sim_args.tsv', sep="\t", header=TRUE)
head(sims)
summary(sims)
dim(sims)


# what are the patterns for the datasets with concordance > 0.8?
ideal <- sims[sims$Name %in% concord$Dataset.File[concord$Concordance >= 0.8], ]
summary(ideal)
head(ideal)
dim(ideal)

aved_ideal <- read.table("site_ave_dnds_concord_by_dataset.csv", sep=",", header=TRUE)
aved_ideal <- merge(x=aved_ideal,
                    y=sims,
                    by.x=c("File"),
                    by.y="Name")
summary(aved_ideal)
head(aved_ideal)
dim(aved_ideal)

#' Only `r nrow(ideal)` / `r nrow(sims)`  =  `nrow(ideal) / nrow(sims)` simulated datasets have concordance >= 0.8
#' 
print(nrow(ideal) / nrow(sims))

#' Find adapter contamination
#' ====
ave_fraglen <- mean(ideal$FragLenAve)
std_fraglen <- mean(ideal$FragLenStd)

# Adapter contamination happens when the fragment size is < read
p_adapterContam <- pnorm(q=READLEN, mean = ave_fraglen, sd = std_fraglen, lower.tail = TRUE, log.p = FALSE)

#' The average adapter contamination in ideal sims = `r p_adapterContam`
#' 

#' Find average length of infection
#' ========
#' 
ave_d_infect <- mean(ideal$Generations)

#' Average days of infection = `r ave_d_infect`  or `r ave_d_infect/365` years


#' Average subs/site Tree Length
#' 
ideal$TreelenList <- sapply(as.character(ideal$TreeLen), function(TreeLen) {strsplit(TreeLen, ",") }, USE.NAMES=FALSE)
ideal$TreelenList <- lapply(ideal$TreelenList, as.numeric)

total_trees <- sum(unlist(lapply(ideal$TreelenList, length)))
sum_sub_persite <- sum(Reduce("+", ideal$TreelenList))
ave_sub_persite <- sum_sub_persite/total_trees

aved_ideal$TreelenList <- sapply(as.character(aved_ideal$TreeLen), function(TreeLen) {strsplit(TreeLen, ",") }, USE.NAMES=FALSE)
aved_ideal$TreelenList <- lapply(aved_ideal$TreelenList, as.numeric)

total_trees <- sum(unlist(lapply(aved_ideal$TreelenList, length)))
sum_sub_persite <- sum(Reduce("+", aved_ideal$TreelenList))
ave_sub_persite <- sum_sub_persite/total_trees

#' Average substitutions per site tree len = `r ave_sub_persite`
#' 
ave_nucsites <- mean(ideal$CodonSites) * 3
#' Average nucleotide sites per simulated population = `r ave_nucsites `

# always 1999 nodes in population
#' How many substitutions per site per node in tree?  ie average branch length  `r ave_sub_persite/1999`


ave_readq <- mean(ideal$MinQual)
#' Read Quality = `r ave_readq`
#' 



#'  We can really associate the entire dataset with concordance, because the tree lengths were randomly assigned.

# Plot ideal vs nonideal concordance
source("load_all_sim_dnds.R")
source("plot_helper.R")
dnds_filename <- "/home/thuy/gitrepo/Umberjack_Benchmark/simulations/out/collate_all.treedist.csv"


load(file="/home/thuy/gitrepo/Umberjack_Benchmark/R/lhs_regression_real_fix/rfe_cont_results_real.RData")
summary(rfe_cont_results_real)

# Use Breiman normality test for z-score
import <- varImp(rfe_cont_results_real)
import$p.value <- 2*pnorm(-abs(import$Overall))
print(import)

dnds <- get_all_sim_dnds(dnds_filename)  
dnds$Is_Ideal <- as.factor(dnds$File %in% ideal$Name)
dnds$Is_Aved_Ideal <- as.factor(dnds$File %in% aved_ideal$File[aved_ideal$concord_est >= 0.8])
summary(dnds)
head(dnds)
dim(dnds)

#' Calculate concordance 95percent CI for ideal datasets
#' 
summary(dnds[dnds$Is_Ideal == TRUE,])
head(dnds[dnds$Is_Ideal == TRUE,])
dim(dnds[dnds$Is_Ideal == TRUE,])
concord <- epi.ccc(dnds$dN_minus_dS.Act[dnds$Is_Ideal == TRUE], dnds$dN_minus_dS.Exp[dnds$Is_Ideal == TRUE])
print(concord$rho.c)

#' Predetermine ideal dnds
preideal <- merge(x=dnds, 
                  y=subset(sims, select=c(Name, Generations, FragLenAve, Cover, NumBreakpoints)),
                  by.x="File", by.y="Name")
partition_ave_treelen <- read.table("/home/thuy/gitrepo/Umberjack_Benchmark/simulations/data/collate_treelen_partition.csv", sep=",", header=TRUE)
summary(partition_ave_treelen)
head(partition_ave_treelen)
dim(partition_ave_treelen)

preideal <- merge(x=preideal, y=partition_ave_treelen,
                  by.x="File", by.y="Name")

#' What is the actual average branch length of all the trees that had concordance > 0.8?
#' 
summary(preideal[preideal$Is_Ideal == TRUE, ])
summary(preideal$AveTreeLen[preideal$Is_Ideal == TRUE]/1999)

preideal$Is_PreIdeal <- preideal$FragLenAve >= 300 & preideal$Cover >= 2 & preideal$TreeDepth.Act/preideal$Reads.Act >= (1e-4 * 5 * 365/1000)
summary(preideal)
head(preideal)
dim(preideal)

summary(preideal[preideal$Is_PreIdeal == TRUE,])
head(preideal[preideal$Is_PreIdeal == TRUE,])
dim(preideal[preideal$Is_PreIdeal == TRUE,])
concord <- epi.ccc(preideal$dN_minus_dS.Act[preideal$Is_PreIdeal == TRUE], preideal$dN_minus_dS.Exp[preideal$Is_PreIdeal == TRUE])
print(concord$rho.c)


## [1] "Using Training Data 
dnds_range <- outlier_range(c(dnds$dN_minus_dS.Act, dnds$dN_minus_dS.Exp))

#+ fig.width=10, fig.height=7
fig <- ggplot(dnds, aes(x=dN_minus_dS.Exp, y=dN_minus_dS.Act, color=Is_Ideal)) + 
  geom_abline(slope=1, color="black") + 
  geom_point(alpha=0.05, shape=1) + 
  geom_smooth(method="lm", se=FALSE, size=2) +   
  xlab("\nExpected dN-dS") + 
  ylab("Inferred dN-dS\n")  +
  theme_bw(base_size=BASE_SIZE) + 
  scale_color_discrete(name="Dataset", labels=c("Non Ideal", "Ideal")) + 
  scale_x_continuous(limits=dnds_range) + 
  scale_y_continuous(limits=dnds_range)
print(fig)

if (!file.exists(paste0(THESIS_OUTDIR, "/umberjack/scatter_expect_infer_idea_nonideal.pdf"))) {
  dir.create(paste0(THESIS_OUTDIR, "/umberjack"), showWarnings = TRUE, recursive = TRUE, mode = "0777")
  ggsave(paste0(THESIS_OUTDIR, "/umberjack/scatter_expect_infer_idea_nonideal.pdf"), fig)
}


# Now plot the random forest Rsquared
lod_dnds_dat <- read.table( '/home/thuy/gitrepo/Umberjack_Benchmark/R/lhs_regression_real_fix/umberjack_accuracy_predict.real.csv', sep=",", header=TRUE)
summary(lod_dnds_dat)
head(lod_dnds_dat)
dim(lod_dnds_dat)

mse <- mean((lod_dnds_dat$residual)^2)
print("MSE for all simulation data")
print(mse)

# Get the RSquared for all of the simulation data predictions
r2 <- rSquared(y=lod_dnds_dat$SqDist_dn_minus_dS, resid=lod_dnds_dat$residual)
print("RSquared for all simulation data")
print(r2)

# NAh, the Rsquared is legit.  But why does the line look so bad?


# http://stats.stackexchange.com/questions/7357/manually-calculated-r2-doesnt-match-up-with-randomforest-r2-for-testing
# Manually calculated rsquared using cor(yhat - y)^2
# cor_dnds <- cor(lod_dnds_dat$SqDist_dn_minus_dS, lod_dnds_dat$pred)
# 
#
# man_r2 <- 1 - sum((lod_dnds_dat$SqDist_dn_minus_dS - lod_dnds_dat$pred)^2)/sum((lod_dnds_dat$SqDist_dn_minus_dS - mean(lod_dnds_dat$SqDist_dn_minus_dS))^2)
# print(man_r2)
# 
# 
# testlm <- lm(pred ~ SqDist_dn_minus_dS, data=lod_dnds_dat)
# summary(testlm)
# 
# testdnds <- lod_dnds_dat
# testdnds$lm_fitted_vals <- testlm$fitted.values
# head(testdnds)
# summary(testdnds)
# 
# 
# fig <- ggplot(testdnds) + 
#   geom_point(aes(x=SqDist_dn_minus_dS, y=pred), alpha=0.1, shape=1, color="black") +
#   geom_line(aes(x=SqDist_dn_minus_dS, y=lm_fitted_vals), color="blue") + 
#   #geom_smooth(method="lm", se=FALSE) + 
#   geom_abline(color="red") +
#   scale_x_continuous(limits=c(0, 3)) + 
#   scale_y_continuous(limits=c(0, 3)) + 
#   xlab("\n Umberjack Inaccuracy") + 
#   ylab("RF Predicted Umberjack Inaccuracy \n") + 
#   theme_bw(BASE_SIZE)
# print(fig)

# Plot the random forest regression fit
dnds_range <- outlier_range(c(lod_dnds_dat$SqDist_dn_minus_dS, lod_dnds_dat$pred))

fig <- ggplot(lod_dnds_dat, aes(x=SqDist_dn_minus_dS, y=pred)) + 
  geom_point(alpha=0.5, shape=1) +
  geom_smooth(method="lm", se=FALSE) + 
  geom_abline(color="red") +
#   scale_x_continuous(limits=c(0, 6000)) + 
#   scale_y_continuous(limits=c(0, 6000)) + 
  xlab(expression(Delta)) + 
  ylab(expression(paste("Random Forest Predicted ", Delta))) + 
  theme_bw(BASE_SIZE)
print(fig)



# From http://stackoverflow.com/questions/5219671/it-is-possible-to-create-inset-graphs-using-ggplot2
# Make inset plot
blowfig <- ggplot(lod_dnds_dat, aes(x=SqDist_dn_minus_dS, y=pred)) + 
  geom_point(alpha=0.1, shape=1) +
  geom_smooth(method="lm", se=FALSE) + 
  geom_abline(color="red") +
  scale_x_continuous(limits=c(0, 3)) + 
  scale_y_continuous(limits=c(0, 3)) +
  theme_bw(BASE_SIZE) + 
  theme(axis.title = element_blank())
print(blowfig)

#A viewport taking up a fraction of the plot area
vp <- viewport(width = 0.4, height = 0.4, x = 0.6, y = 0.125, just=c("left", "bottom"))

#Just draw the plot twice
png(paste0(THESIS_OUTDIR, "/umberjack/RandomForestRegressionRsq_real.png"), width=560, height=560)
print(fig)
print(blowfig, vp=vp)
dev.off() 
#ggsave(filename=paste0(THESIS_OUTDIR, "/umberjack/RandomForestRegressionRsq_real.pdf"), plot=fig, device=pdf)


# Load in the full training data (not real dataset)
# Nope!  Looks like I saved th wrong .RData file  :(
#load(file="/home/thuy/gitrepo/Umberjack_Benchmark/R/lhs_regression_fix/rfe_cont_results.RData")
#summary(rfe_cont_results)

# Load in the importances manually
#importnonreal <- varImp(rfe_cont_results)
importnonreal <- read.table('/home/thuy/gitrepo/Umberjack_Benchmark/R/lhs_regression_fix/importance.csv', sep=",", header=TRUE)
# use the increase in MSE () normalized by stddev as zscore
# See http://www.statistik.uni-dortmund.de/useR-2008/slides/Strobl+Zeileis.pdf
importnonreal$p.value <- 2*pnorm(-abs(importnonreal$IncMSE))
print(importnonreal)

library(xtable)
xtable(importnonreal, 
       display=c("s", "s", "g", "E"), 
       digits=10,
       caption="Importance of Feature on Prediction of Umberjack Error.", type="latex")



#' train random Forest
#' 
#cleandnds <- na.omit(dnds)
respname <- "SqDist_dn_minus_dS"
feats <- LM_COVAR_NAMES
ids <- c("File", "Window_Start", "Window_End", "CodonSite")
cleandnds <- na.omit(dnds[, c(respname, feats, ids)])
summary(cleandnds)
dim(cleandnds)

keep_feats <- c("UnambigCodonRate.Act",
                "AADepth.Act",
                "EntropyCodon.Act",
                "ErrPerCodon.Act",
                "Subst.Act",
                "N.Act",
                "S.Act" ,
                "ES.Act",
                "TreeLenPerRead.Act",
                "TreeDistPerRead.Act",
                "P_SameCodonFreq.Act",
                "EntropyCodon.Exp")

ptm <- proc.time()
final_rf <- randomForest(x=cleandnds[, keep_feats], y=cleandnds[, respname],
             replace=TRUE,
             keep.inbag=TRUE, 
             importance=TRUE, 
             keep.forest=TRUE)
print (proc.time() - ptm)
save(final_rf, file="./lhs_regression_fix/final_rf.RData")
print(final_rf)

fitted.values <- predict(final_rf, cleandnds[, keep_feats])
plot(cleandnds[, respname],  fitted.values,col="green")

