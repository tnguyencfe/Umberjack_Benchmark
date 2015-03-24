# Uses Results from Typical Reads vs Good Reads to determine covariates for dN/dS inferrence accuracy


#+ setup, include=FALSE
library(knitr)
opts_chunk$set(warning=FALSE, message=FALSE, width=1200)


# From data from all windows, aggregates by averaging over windows
library(ggplot2)
library(reshape2)
library(epiR)
library(plyr)
library(stats)
library(psych)
library(GGally)
library(MASS)
library(caret)
library(gplots)
library(RColorBrewer)
library(nortest)
library(devtools)
library(ggbiplot)
library(directlabels)
library(gamlss)
library(logistf)

NUC_PER_CODON <- 3


# Read locations of input files from local sliding_window_tree_unit_test.config file
CONFIG_FILENAME <- "./find_covar_accuracy.config"
config<-read.table(CONFIG_FILENAME, sep="=", col.names=c("key","value"), as.is=c(1,2))

COLLATE_DNDS_FILENAME <- config[config$key=="COLLATE_DNDS_FILENAME",]$val
COLLATE_ERR_FREE_DNDS_FILENAME <- config[config$key=="COLLATE_ERR_FREE_DNDS_FILENAME",]$val
COLLATE_ERR_FREE_DOWNSAMP_DNDS_FILENAME <- config[config$key=="COLLATE_ERR_FREE_DOWNSAMP_DNDS_FILENAME",]$val
COLLATE_ERR_FREE_DOWNSAMP_MASK_DNDS_FILENAME <- config[config$key=="COLLATE_ERR_FREE_DOWNSAMP_MASK_DNDS_FILENAME",]$val


EXPECTED_DNDS_FILENAME <-  config[config$key=="EXPECTED_DNDS_FILENAME",]$val



FULL_POPN_CONSERVE_CSV <- config[config$key=="FULL_POPN_CONSERVE_CSV",]$val
ORIG_CONSERVE_CSV <- config[config$key=="ORIG_CONSERVE_CSV",]$val
ALN_CONSERVE_CSV <- config[config$key=="ALN_CONSERVE_CSV",]$val
ORIG_ERR_FREE_CONSERVE_CSV <- config[config$key=="ORIG_ERR_FREE_CONSERVE_CSV",]$val
ALN_ERR_FREE_CONSERVE_CSV <- config[config$key=="ALN_ERR_FREE_CONSERVE_CSV",]$val


SMOOTH_DIST <-  as.numeric(config[config$key=="SMOOTH_DIST",]$val)

#+ results='asis'
kable(config, format="html", caption="config")

#+
# Cols:  Window_Start,Window_End,Reads,CodonSite,CodonDepth,Conserve,Entropy,N,S,EN,ES,dN,dS,dN_minus_dS
get_collate <- function(COLLATE_DNDS_FILENAME) {
  collate_dnds <- read.table(COLLATE_DNDS_FILENAME, header=TRUE, na.strings="None", comment.char = "#", sep=",")
  collate_dnds$Subst <- collate_dnds$N + collate_dnds$S
  collate_dnds$Diverge <- 1 - collate_dnds$Conserve
  colnames(collate_dnds)[grep("^Codons$", colnames(collate_dnds), perl=TRUE)] <- "CodonDepth"
  collate_dnds$dNdS <- collate_dnds$dN/collate_dnds$dS
  collate_dnds$dNdS[collate_dnds$dS==0] <- NA
  collate_dnds$ErrPerc <- collate_dnds$Err*100/(collate_dnds$Reads)  # per-window-codon site error fraction
  collate_dnds$AmbigCodonFractn <- 1 - collate_dnds$CodonDepth/collate_dnds$Reads
  collate_dnds$Err_N_Perc <- collate_dnds$Err_N*100/collate_dnds$CodonDepth
  collate_dnds$Err_S_Perc <- collate_dnds$Err_S*100/collate_dnds$CodonDepth 
  collate_dnds$Pad_Fractn <- collate_dnds$Pad/(collate_dnds$Reads)  # per-window-codon site error fraction
  
  # Average across all codon sites in a window
  per_window_ave <- ddply(.data=collate_dnds, .variables="Window_Start", 
                          .fun=function(x) {   
                           
                            data.frame(Window_Diverge=mean(x$Diverge, na.rm=TRUE),
                                       Window_Entropy=mean(x$Entropy, na.rm=TRUE),
                                       Window_Subst=mean(x$Subst, na.rm=TRUE),
                                       Window_CodonDepth=mean(x$CodonDepth, na.rm=TRUE),
                                       Window_ErrPerc=mean(x$ErrPerc, na.rm=TRUE),
                                       Window_AmbigCodonFractn=mean(x$AmbigCodonFractn, na.rm=TRUE))
                          })
  collate_dnds <- merge(x=collate_dnds, y=per_window_ave, by="Window_Start", all=TRUE, sort=TRUE)
  
  return (collate_dnds)
}


#' **Summary Per-Window-Codon Stats**
#' 
collate_dnds <- get_collate(COLLATE_DNDS_FILENAME)
dim(collate_dnds)
head(collate_dnds)
summary(collate_dnds)

#' **Summary Err Free Per-Window-Codon Stats**
#'
collate_errfree_dnds <- get_collate(COLLATE_ERR_FREE_DNDS_FILENAME) 
dim(collate_errfree_dnds)
head(collate_errfree_dnds)
summary(collate_errfree_dnds)

#' **Summary Err Free Downsampled to same Read Depth as Typical Per-Window-Codon Stats**
#'
collate_errfree_down_dnds <- get_collate(COLLATE_ERR_FREE_DOWNSAMP_DNDS_FILENAME) 
dim(collate_errfree_down_dnds)
head(collate_errfree_down_dnds)
summary(collate_errfree_down_dnds)

#' **Summary Err Free Downsampled and N-masked same as Typical Per-Window-Codon Stats**
#'
collate_errfree_down_mask_dnds <- get_collate(COLLATE_ERR_FREE_DOWNSAMP_MASK_DNDS_FILENAME) 
dim(collate_errfree_down_mask_dnds)
head(collate_errfree_down_mask_dnds)
summary(collate_errfree_down_mask_dnds)

#'
#'  Read in Expected dN/dS
#'  ==============================
#'  

# Cols: Observed S Changes  Observed NS Changes  E[S Sites]  E[NS Sites]	Observed S. Prop.	P{S}	dS	dN	dN-dS	P{S leq. observed}	P{S geq. observed}	Scaled dN-dS	dn/ds
# Parse the expected dnds filename for it start and end nucleotide positions (1-based)
expected_dnds <- read.table(EXPECTED_DNDS_FILENAME, header=TRUE, sep="\t")  # Site  Interval	Scaling_factor	Rate_class	Omega
expected_dnds$CodonSite <- as.numeric(rownames(expected_dnds))
expected_dnds$dNdS <- expected_dnds$dN/expected_dnds$dS
expected_dnds$dNdS[expected_dnds$dS == 0] <- NA
colnames(expected_dnds)[grep("Observed.S.Changes", colnames(expected_dnds))] <- "S"
colnames(expected_dnds)[grep("Observed.NS.Changes", colnames(expected_dnds))] <- "N"
colnames(expected_dnds)[grep("E.S.Sites.", colnames(expected_dnds))] <- "ES"
colnames(expected_dnds)[grep("E.NS.Sites.", colnames(expected_dnds))] <- "EN"
colnames(expected_dnds)[grep("Scaled.dN.dS", colnames(expected_dnds))] <- "dN_minus_dS"
expected_dnds$Subst <- expected_dnds$S + expected_dnds$N


#' **Summary Expected dN/dS Stats**
#' 
dim(expected_dnds)
head(expected_dnds)
summary(expected_dnds)

# Takes a *.conserve.csv file with nucleotide stats for unsliced sequences
# Expects cols: "NucSite", "Conserve", "Entropy", "NucDepth", "CodonDepth"  
# Returns a dataframe with elements:
# - CodonDat = dataframe with values average per-codon values for Conserve, Diverge, Entropy, CodonDepth
# - NucStats = matrix summary of per-nucleotide stats
fill_nowindow_codon_conserve_dat <- function(nuc_conserve_csv) {
  # Per-nucleotide conservation, Divergence, Entropy
  conserve_nuc <- read.table(nuc_conserve_csv, header=TRUE, sep=",")
  conserve_nuc$Diverge <- 1 - conserve_nuc$Conserve
  
  NUM_NUC_SITES <- max(conserve_nuc$NucSite)
  NUM_CODON_SITES <- ceiling(NUM_NUC_SITES / 3)
  
  # Get per-codon average stats for Divergence, Entropy
  conserve_codon <- adply(.data=c(1:NUM_CODON_SITES), .margins=1,
                          .fun=function(codonsite) {
                            start_nuc_site <- (codonsite * 3) -2  # 1-based nucleotide sites
                            end_nuc_site <- codonsite * 3
                            codon_set <- conserve_nuc[conserve_nuc$NucSite >= start_nuc_site & conserve_nuc$NucSite <= end_nuc_site,]
                            data.frame(  
                              CodonDepth=conserve_nuc[start_nuc_site, "CodonDepth"],
                              Conserve=mean(codon_set$Conserve),
                              Diverge=mean(codon_set$Diverge),
                              Entropy=mean(codon_set$Entropy))
                          })
  colnames(conserve_codon)[1] <- "CodonSite"  # rename X1 column to CodonSite
  conserve_codon$CodonSite <- as.numeric(conserve_codon$CodonSite)
  returnDat <- list(CodonDat=conserve_codon,
                    NucStat=summary(conserve_nuc))
  return (returnDat)
}

full_popn_conserve_dat <- fill_nowindow_codon_conserve_dat(FULL_POPN_CONSERVE_CSV)

#' What covariates determine the distance from the actual dN/dS to expected dN/dS?
#' --------------------------------------------------------------------------------------
#' 
#'  What are the patterns when the actual dn/ds is far from the true dn/ds?
#'  Find distance from actual dn/ds to true dn/ds.
#'  What is that distance by Full Population N, S, Subst?
#'  What is the distance by Full Population Site Diversity?
#'  what is the distance by Full Population Site Entropy?
#'  What is the distance by Nucleotide Site?
#'  What is the distance by Codon Depth?
#'  What is the distance by Window Entropy?
#'  What is the distance by Window Diversity?
#' Do linear modelling to find the covariates.




collate_subset <- subset(collate_dnds, select=c(Window_Start, CodonSite, 
                                                Reads, CodonDepth, 
                                                Diverge, Entropy, 
                                                AmbigBase, ErrPerc, Err_N_Perc, Err_S_Perc, AmbigCodonFractn,
                                                N, S, Subst, EN, ES, dNdS,  dN_minus_dS,
                                                Window_Diverge, Window_Entropy, 
                                                Window_Subst, Window_CodonDepth, Window_ErrPerc, Window_AmbigCodonFractn))
collate_errfree_subset <- subset(collate_errfree_dnds, select=c(Window_Start, CodonSite, 
                                                                Reads, CodonDepth, 
                                                                Diverge, Entropy, 
                                                                AmbigBase, ErrPerc, Err_N_Perc, Err_S_Perc, AmbigCodonFractn,
                                                                N, S, Subst, EN, ES, dNdS,  dN_minus_dS,
                                                                Window_Diverge, Window_Entropy, 
                                                                Window_Subst, Window_CodonDepth, Window_ErrPerc, Window_AmbigCodonFractn))
collate_errfree_down_subset <- subset(collate_errfree_down_dnds, select=c(Window_Start, CodonSite, 
                                                                Reads, CodonDepth, 
                                                                Diverge, Entropy, 
                                                                AmbigBase, ErrPerc, Err_N_Perc, Err_S_Perc, AmbigCodonFractn,
                                                                N, S, Subst, EN, ES, dNdS,  dN_minus_dS,
                                                                Window_Diverge, Window_Entropy, 
                                                                Window_Subst, Window_CodonDepth, Window_ErrPerc, Window_AmbigCodonFractn))
collate_errfree_down_mask_subset <- subset(collate_errfree_down_mask_dnds, select=c(Window_Start, CodonSite, 
                                                                          Reads, CodonDepth, 
                                                                          Diverge, Entropy, 
                                                                          AmbigBase, ErrPerc, Err_N_Perc, Err_S_Perc, AmbigCodonFractn,
                                                                          N, S, Subst, EN, ES, dNdS,  dN_minus_dS,
                                                                          Window_Diverge, Window_Entropy, 
                                                                          Window_Subst, Window_CodonDepth, Window_ErrPerc, Window_AmbigCodonFractn))

collate_all <- rbind(data.frame(Source="Typical", collate_subset),
                     data.frame(Source="ErrFree", collate_errfree_subset),
                     data.frame(Source="ErrFreeDown", collate_errfree_down_subset),
                     data.frame(Source="ErrFreeDownMask", collate_errfree_down_mask_subset))

colnames(collate_all)[-c(1,2,3)] <- paste0(colnames(collate_all)[-c(1,2,3)], ".Act")
summary(collate_all)

exp_subset <- subset(expected_dnds, select=c(CodonSite, N, S, Subst, EN, ES, dNdS,  dN_minus_dS))
colnames(exp_subset)[-1] <- paste0(colnames(exp_subset)[-1], ".Exp")



collate_exp_dnds <- merge(x=collate_all, 
                          y=exp_subset,
                          by="CodonSite", all.x=TRUE, all.y=FALSE)

full_popn_conserve_subset <- subset(full_popn_conserve_dat$CodonDat, select=-c(CodonDepth, Conserve))
colnames(full_popn_conserve_subset)[-1] <- paste0(colnames(full_popn_conserve_subset)[-1], ".Exp")

collate_exp_dnds <- merge(x=collate_exp_dnds, y=full_popn_conserve_subset, by="CodonSite", all.x=TRUE, all.y=FALSE)

collate_exp_dnds$LOD_dNdS <- log(collate_exp_dnds$dNdS.Act + .Machine$double.eps) - log(collate_exp_dnds$dNdS.Exp + .Machine$double.eps)
collate_exp_dnds$Dist_dn_minus_dS <- abs(collate_exp_dnds$dN_minus_dS.Act - collate_exp_dnds$dN_minus_dS.Exp)
collate_exp_dnds$Dist_dNdS <- abs(collate_exp_dnds$dNdS.Act - collate_exp_dnds$dNdS.Exp)
collate_exp_dnds$Low_Dist_dNdS <- factor(collate_exp_dnds$Dist_dNdS < 0.2, exclude=NA)
collate_exp_dnds$Low_Dist_dn_minus_dS <- factor(collate_exp_dnds$Dist_dn_minus_dS < 0.2, exclude=NA)
#collate_exp_dnds$PercDist_dNdS <- collate_exp_dnds$Dist_dNdS/(collate_exp_dnds$dNdS.Exp + .Machine$double.eps)
#collate_exp_dnds$PercDist_dn_minus_dS <- collate_exp_dnds$Dist_dn_minus_dS/(collate_exp_dnds$dN_minus_dS.Exp + .Machine$double.eps)
collate_exp_dnds$LowLOD <- factor(abs(collate_exp_dnds$LOD_dNdS) < 0.5, exclude=NA)
#collate_exp_dnds$LowPercOff_dNdS <- as.factor(collate_exp_dnds$PercDist_dNdS < 0.1)
#collate_exp_dnds$LowPercOff_dN_minus_dS <- as.factor(collate_exp_dnds$PercDist_dn_minus_dS < 0.1)
summary(collate_exp_dnds)
dim(collate_exp_dnds)
head(collate_exp_dnds)

NUM_RESP_NAMES <- c("LOD_dNdS", "Dist_dn_minus_dS", "Dist_dNdS")
CAT_RESP_NAMES <- c("Low_Dist_dNdS", "Low_Dist_dn_minus_dS", "LowLOD")
COVAR_NAMES <- colnames(collate_exp_dnds[sapply(collate_exp_dnds,is.numeric)])[!colnames(collate_exp_dnds[sapply(collate_exp_dnds,is.numeric)]) %in% NUM_RESP_NAMES]



# Check normality for Response
qqnorm(collate_exp_dnds$LOD_dNdS[!is.na(collate_exp_dnds$LOD_dNdS)])
ad.test(collate_exp_dnds$LOD_dNdS[!is.na(collate_exp_dnds$LOD_dNdS)])

qqnorm(collate_exp_dnds$Dist_dn_minus_dS[!is.na(collate_exp_dnds$Dist_dn_minus_dS)])
ad.test(collate_exp_dnds$Dist_dn_minus_dS[!is.na(collate_exp_dnds$Dist_dn_minus_dS)])

qqnorm(collate_exp_dnds$Dist_dNdS[!is.na(collate_exp_dnds$Dist_dNdS)])
qqnorm(collate_exp_dnds$Dist_dNdS[!is.na(collate_exp_dnds$Dist_dNdS) & collate_exp_dnds$Dist_dNdS < 10])
ad.test(collate_exp_dnds$Dist_dNdS[!is.na(collate_exp_dnds$Dist_dNdS)])



#' Correlation Between Features
#' -----------------------------------
#' 

# Find correlation between features.  Don't bother scaling.  The correlation heatmap is the same whether we scale & centre or not.
corMat <- cor(collate_exp_dnds[sapply(collate_exp_dnds,is.numeric)], use="complete.obs", method="spearman")
#+ fig.width=15, fig.height=15
heatmap.2(corMat,  col=bluered, density.info="none", trace="none", srtCol=45, main="Feature Correlation", margins=c(12,12), 
          colsep=(1:ncol(corMat)), rowsep=(1:nrow(corMat)), sepwidth=c(0.05, 0.05), sepcolor="black")



#' Density Plots of Each Numerical Column
#' =============================================
#' 
#' 
plot_density <- function(colname) {
  fig <- ggplot(collate_exp_dnds, aes_string(x=colname)) + 
    geom_density(na.rm=TRUE, color="black") + 
    geom_density(aes_string(x=colname, color="Source"), na.rm=TRUE) + 
    ggtitle(paste0("Density plot ", colname))
  print(fig)
}

#sapply(COVAR_NAMES, plot_density)


#' Plot Response Vs Variable
#' =============================================
#' 
#' 

plot_resp_vs_var <- function(resp_colname) {
  sapply(COVAR_NAMES, 
         function(var_colname) {
          fig <- ggplot(collate_exp_dnds, aes_string(x=var_colname, y=resp_colname)) +             
            geom_point(shape=1, alpha=0.5, aes_string(color="Source"), na.rm=TRUE) + 
            ggtitle(paste0(resp_colname, " vs ", var_colname))
          print(fig)
  })
}
plot_resp_vs_var("Dist_dNdS")


#' Verify Homoscedasticity.
#' =============================================
#' 
#' 
#' **Plot Mean-Variance of Codon Site Response.  Each Window is a sample.**
#' 
#' 
plot_resp_meanvar <- function(resp_colname) {
  meanvar <- ddply(.data=collate_exp_dnds, .variables=c("CodonSite", "Source"), .fun=function(x) {
                  returndat <- data.frame(site_ave=mean(x[, resp_colname], na.rm=TRUE),
                                          site_var=var(x[, resp_colname], na.rm=TRUE))
                  return (returndat)
                   })
  print(summary(meanvar))
  fig <- ggplot(meanvar, aes(x=site_ave, y=log(site_var), color=Source)) +     
    geom_smooth(method="loess", na.rm=TRUE) + 
    geom_point(shape=1, alpha=0.5, na.rm=TRUE) +     
    ggtitle(paste0("Mean-Variance Relationship for Reponse ", resp_colname))
  print(fig)
  

  fig <- ggplot(meanvar[abs(meanvar$site_ave-mean(meanvar$site_ave, na.rm=TRUE)) < 2*sd(meanvar$site_ave, na.rm=TRUE), ]) +     
    geom_smooth(aes(x=site_ave, y=site_var), method="loess", color="black", na.rm=TRUE) +   
    geom_smooth(aes(x=site_ave, y=site_var, color=Source), method="loess", na.rm=TRUE) +   
    geom_point(aes(x=site_ave, y=site_var, color=Source), shape=1, alpha=0.5, na.rm=TRUE) + 
    ggtitle(paste0("Mean-Variance Relationship for Reponse ", resp_colname, " Ignore outliers"))
  print(fig)

  log_meanvar <- ddply(.data=collate_exp_dnds, .variables=c("CodonSite", "Source"), .fun=function(x) {
    returndat <- data.frame(site_ave=mean(log(x[, resp_colname] + 0.000001), na.rm=TRUE),
                            site_var=var(log(x[, resp_colname] + 0.000001), na.rm=TRUE))
    return (returndat)
  })
  print(summary(log_meanvar))
  fig <- ggplot(log_meanvar[abs(log_meanvar$site_ave-mean(log_meanvar$site_ave, na.rm=TRUE)) < 2*sd(log_meanvar$site_ave, na.rm=TRUE), ]) +     
    geom_smooth(aes(x=site_ave, y=site_var), method="loess", color="black", na.rm=TRUE) +   
    geom_smooth(aes(x=site_ave, y=site_var, color=Source), method="loess", na.rm=TRUE) +   
    geom_point(aes(x=site_ave, y=site_var, color=Source), shape=1, alpha=0.5, na.rm=TRUE) + 
    ggtitle(paste0("Log Mean-Variance Relationship for Reponse ", resp_colname, " Ignore outliers"))
  print(fig)
}

#sapply(NUM_RESP_NAMES, plot_resp_meanvar)


#' **Nonparametric Test for homoscedasticity.  They all fail homoscedasticity.**
#' 
test <- fligner.test(LOD_dNdS~CodonSite, data=collate_exp_dnds, na.action=na.exclude)
print(test)
test <- fligner.test(Dist_dn_minus_dS~CodonSite, data=collate_exp_dnds, na.action=na.exclude)
print(test)
test <- fligner.test(Dist_dNdS~CodonSite, data=collate_exp_dnds, na.action=na.exclude)
print(test)



# #' Impact of Principal Components
# #' =====================================
# #' 
# # Centre all columns to 0 and scale such that variance is 1 before principal component
# # Get rid of missing data since prcomp seems to choke on them
# pc_collate_exp_dnds <- na.omit(collate_exp_dnds)
# pcs <- prcomp(subset(pc_collate_exp_dnds, select=c(COVAR_NAMES)), 
#               centre=TRUE, scale.=TRUE)
# summary(pcs)
# 
# # Variances of principal components
# variances <- data.frame(variances=pcs$sdev**2, pcomp=1:length(pcs$sdev))
# 
# 
# #Plot of variances
# fig <- ggplot(variances, aes(pcomp, variances)) + 
#   geom_bar(stat="identity") + 
#   ggtitle("Variances of Each Prinicipal Component")
# print(fig)
# 
# # Cumulative Distro of the pcomp
# fig <- ggplot(variances, aes(x=pcomp)) + 
#   geom_line(aes(y = cumsum(variances)/sum(variances))) + 
#   ylab("Cumulative Distribution of Variance")
#   ggtitle("Cumulative Distribution of Variances of Each Prinicipal Component")
# print(fig)
# 
# # biplot showing projection of original variables on 1st 2 principal components
# fig <- ggbiplot(pcs, obs.scale = 1, var.scale=1, groups=pc_collate_exp_dnds$LowLOD, alpha=0.1, 
#                 ellipse=TRUE,circle=FALSE,varname.size=6)
# print(fig)


#' Why is there such a big spread in N.Act for Typical reads but not for Error Free Downsampled N.Mask Reads?
#' ================================================================
#' 
#' 

#' **Even when we get rid of the low codon depth windows, we still see the spread in N.Act for Typical  Reads**
#' 
collate_exp_dnds_Ncheck <- na.omit(collate_exp_dnds[collate_exp_dnds$Source %in% c("Typical", "ErrFreeDownMask"), ])
summary(collate_exp_dnds_Ncheck)


#' **Little difference in N.Act for inaccurate dN/dS and accurate dN/dS is very small for ErrFreeDownMask reads, but in general these reads have mostly accurate dN/dS**.
#' **Big gap in N.Act for inaccurate dN/dS and accurate dN/dS for Typical reads**.
#' 
# Plot difference between Typical, ErrFreeDownMask reads vs Full Population
plot_fullpopn_v_actual_lowdist_dnds <- function(colname) {
  fig <- ggplot(collate_exp_dnds_Ncheck, aes_string(x=paste0(colname, ".Exp"), y=paste0(colname, ".Act"), color="Low_Dist_dNdS")) + 
    geom_abline(slope=1, color="black") + 
    geom_point(shape=1, alpha=0.5) + 
    geom_smooth(method="loess") +   
    xlab(paste0("Full Population ", colname)) + 
    ylab(paste0("Actual ", colname)) + 
    ggtitle(paste0("Expected Full Population ", colname, " vs Actual Window-CodonSite ", colname)) +
    facet_wrap(~Source)
  print(fig) 
}

plot_fullpopn_v_actual_lowdist_dnds("N")
plot_fullpopn_v_actual_lowdist_dnds("S")
plot_fullpopn_v_actual_lowdist_dnds("Entropy")
plot_fullpopn_v_actual_lowdist_dnds("Diverge")




#' Scatterplot Typical vs ErrFreeDownMask reads
#' 
#' ???? Perhaps we matched the same number of N's but we didn't match the same left and right pads????
#' This would mean that the error free reads have fewer ambiguous codons.  But this shows that they have more ambiguous codons in a window.
#' 
#' ??? Perhaps the average read length in a window is shorter for error free reads because they don't align as well?
#' ??? Perhaps we need to extend the error free reads so that they are just as long as the typical reads? 
#' NO -- this won't really do anything for us.  It will just make them more accurate. and we want to know why typical are so crappy.
#' 
#' Typical reads have fewer window ambiguous codons.  Why is that?  They should have the same number of Ns and padding!
#' Typical reads sometimes have every so slightly more codon depth.  ??? they should have same number of Ns and padding?
#' Difference between typical and errfree divergence & entropy is so low, that it doesn't seem to matter
collate_typical_errfreedownmask <- merge(x=collate_dnds, y=collate_errfree_down_mask_dnds, 
                                         by=c("Window_Start", "Window_End", "CodonSite"),  all=TRUE,
                                         suffixes=c(".Typ", ".ErrFreeDownMask"))
collate_typical_errfreedownmask <- merge(x=collate_typical_errfreedownmask, 
                                         y=subset(expected_dnds, select=c(CodonSite, dNdS)),
                                         by="CodonSite")
colnames(collate_typical_errfreedownmask)[grep("^dNdS$", colnames(collate_typical_errfreedownmask), perl=TRUE)] <- "dNdS.Exp"
collate_typical_errfreedownmask$Dist_dNdS.Factor.Typ <- NA
collate_typical_errfreedownmask$Dist_dNdS.Factor.Typ[
  collate_typical_errfreedownmask$dNdS.Exp - collate_typical_errfreedownmask$dNdS.Typ >= 0.2] <- "TOOLOW"
collate_typical_errfreedownmask$Dist_dNdS.Factor.Typ[
  collate_typical_errfreedownmask$dNdS.Typ - collate_typical_errfreedownmask$dNdS.Exp >= 0.2] <- "TOOHI"
collate_typical_errfreedownmask$Dist_dNdS.Factor.Typ[
  abs(collate_typical_errfreedownmask$dNdS.Typ - collate_typical_errfreedownmask$dNdS.Exp) < 0.2] <- "GOOD"
collate_typical_errfreedownmask$Dist_dNdS.Factor.Typ <- factor(collate_typical_errfreedownmask$Dist_dNdS.Factor.Typ, exclude=NA)
summary(collate_typical_errfreedownmask)
plot_typical_errfreedownmask <- function(colname) {
  concord <- epi.ccc(collate_typical_errfreedownmask[, paste0(colname, ".ErrFreeDownMask")],
                     collate_typical_errfreedownmask[, paste0(colname, ".Typ")])
  
  fig <- ggplot(collate_typical_errfreedownmask[!is.na(collate_typical_errfreedownmask[, paste0(colname, ".ErrFreeDownMask")]) &
                                                  !is.na(collate_typical_errfreedownmask[, paste0(colname, ".Typ")]) &
                                                  !is.na(collate_typical_errfreedownmask[, "Dist_dNdS.Factor.Typ"]), ],
                aes_string(x=paste0(colname, ".ErrFreeDownMask"), y=paste0(colname, ".Typ"), color="Dist_dNdS.Factor.Typ")) + 
    geom_abline(slope=1) + 
    geom_smooth(method="lm", na.rm=TRUE) + 
    geom_point(shape=1, alpha=0.3, na.rm=TRUE) + 
    xlab(paste0("Error Free Downsampled N-Masked ", colname)) + 
    ylab(paste0("Typical ", colname)) + 
    ggtitle(paste0("Error Free Downsampled N-Masked  vs Typical ", colname, " Concordance=", concord$rho$est))
  print(fig) 
}
sapply(colnames(collate_dnds)[-c(1,2,4)], plot_typical_errfreedownmask)

#' **After N.Exp > 30, either the N.Act for Typical reads are good or they are way too low.  Why?**
#' 
 
bound_errFreeDownMask <- ddply(.data=collate_exp_dnds_Ncheck[collate_exp_dnds_Ncheck$Source=="ErrFreeDownMask",], 
                               .variables="N.Exp", .fun=function(x) {
  data.frame(Max.ErrFreeDownMask.N.Act=max(x$N.Act, na.rm=TRUE), 
             Min.ErrFreeDownMask.N.Act=min(x$N.Act, na.rm=TRUE))
  })

pred_Typical_N.Act_low <- merge(x=collate_exp_dnds_Ncheck, y=bound_errFreeDownMask, by="N.Exp", all=TRUE, sort=FALSE)
summary(pred_Typical_N.Act_low)


pred_Typical_N.Act_low$Is_CrapLow_N.Act <- factor(pred_Typical_N.Act_low$Min.ErrFreeDownMask.N.Act > pred_Typical_N.Act_low$N.Act, exclude=NA)
pred_Typical_N.Act_low <- subset(pred_Typical_N.Act_low, 
                                 subset=(Source == "Typical" & N.Act <= Max.ErrFreeDownMask.N.Act))
summary(pred_Typical_N.Act_low)
dim(pred_Typical_N.Act_low)


# Violin Plot N.Act vs Covariates for Typical Reads when N.Act lower than ErrorFree Downsampled Masked Reads
# The gap between expected N.Exp and crappy low N.Act widens with:
# Entropy.Exp - higher expected entropy
# Diverge.Exp - higher expected diversity
# dNdS.Exp - higher means crappier lower N.Act
# Subst.Exp - higher expected substitutions = crappier lower N.Act
# S.Exp - lower S.Exp = crappier lower N.Act
# N.Exp - higher N.Exp = crappier lower N.Act
# Window_AmbigCodonFracn.Act  - lower means crappier lower N.Act.  ??? Does crappier lower N.Act correspond with inaccurate dn/ds?  No.
# Window_ErrPerc.Act - higher window errors = crappier lower N.Act.  ??? Expect higher window errors to mean higher N.Act albeit inaccurate N.Act
# Window_Subst.Act - ?????more subtitutions in the window, more likely to get crappier N.Act
# Subst.Act - lower subst means crappier lower n.act
# S.Act - lower S.act means crappier lower n.act
# AmbigBase.Act - lower ambigBase
# Entropy.Act - lower entropy  = lower crappier N.Act
# Diverge.Act - lower divergence = lower crappier N.Act
sapply(COVAR_NAMES,
       function(colname) {
         fig <- ggplot(pred_Typical_N.Act_low, aes_string(x="Is_CrapLow_N.Act", y=colname, color="Is_CrapLow_N.Act")) +
           geom_point(shape=1, alpha=0.5) + 
           geom_violin() +
           stat_summary(fun.ymin=function(y) {quantile(y, probs=c(0.25))},
                        fun.y=function(y){quantile(y, probs=c(0.5))},
                        fun.ymax=function(y) {quantile(y, probs=c(0.75))}, 
                        geom='crossbar') +
           xlab("Crappy and Low N.Act?") + 
           ylab(colname) + 
           ggtitle(paste0("Violin plot Crappy N.Act vs ", colname))
         print(fig)
       })


#' Find significance of predictors for crappy extra low N.Act for Typical Reads when N.Exp > 30
#' 

# Manual Logistic regression.  What determines crappy or not crappy for Typical Reads when N.Exp > 30?

#' Compare dN_minus_dS.Exp vs dNdS.Exp as covariates:  dN_minus_dS.Exp smaller deviance, smaller AIC
gamlssfit <- gamlss(Is_CrapLow_N.Act~dN_minus_dS.Exp, data=pred_Typical_N.Act_low,family=BI)
summary(gamlssfit)
plot(gamlssfit)
EDA(residuals(gamlssfit))

gamlssfit <- gamlss(Is_CrapLow_N.Act~dNdS.Exp, data=pred_Typical_N.Act_low, family=BI)
summary(gamlssfit)
plot(gamlssfit)
EDA(residuals(gamlssfit))

#' Compare different links:  they're all the same for dN_minus_dS.Exp, but normal dN_minus_dS.Exp is best with lowest deviance
fit <-  glm(Is_CrapLow_N.Act~dN_minus_dS.Exp, data=pred_Typical_N.Act_low, family=binomial(link="logit"))
summary(fit)

fit <-  glm(Is_CrapLow_N.Act~dN_minus_dS.Exp, data=pred_Typical_N.Act_low, family=quasibinomial)
summary(fit)

fit <-  gamlss(Is_CrapLow_N.Act~dN_minus_dS.Exp, data=pred_Typical_N.Act_low, family=BB)
summary(fit)



# Entropy.Exp - higher expected entropy
# Diverge.Exp - higher expected diversity
# dNdS.Exp - higher means crappier lower N.Act
# Subst.Exp - higher expected substitutions = crappier lower N.Act
# S.Exp - lower S.Exp = crappier lower N.Act
# N.Exp - higher N.Exp = crappier lower N.Act
# Window_AmbigCodonFracn.Act  - lower means crappier lower N.Act.  ??? Does crappier lower N.Act correspond with inaccurate dn/ds?  No.
# Window_ErrPerc.Act - higher window errors = crappier lower N.Act.  ??? Expect higher window errors to mean higher N.Act albeit inaccurate N.Act
# Window_Subst.Act - ?????more subtitutions in the window, more likely to get crappier N.Act
# Subst.Act - lower subst means crappier lower n.act
# S.Act - lower S.act means crappier lower n.act
# AmbigBase.Act - lower ambigBase
# Entropy.Act - lower entropy  = lower crappier N.Act
# Diverge.Act - lower divergence = lower crappier N.Act

fit <- glm(Is_CrapLow_N.Act~dNdS.Exp+N.Exp+S.Exp+Entropy.Exp+Diverge.Exp+
             (Window_AmbigCodonFractn.Act*Reads.Act)+Window_ErrPerc.Act+Window_Subst.Act+Window_Entropy.Act+Window_Diverge.Act+
             Subst.Act+S.Act+
             AmbigBase.Act+Entropy.Act+Diverge.Act+
             ErrPerc.Act+Err_N_Perc.Act+Err_S_Perc.Act,
           data=pred_Typical_N.Act_low, family=binomial(link="logit"))
summary(fit)
#' Don't trust this model since we get separation warnings.  This is because Subst.Act+S.Act are highly linked to N.Act


#' Remove S.Act
fit <- glm(Is_CrapLow_N.Act~dNdS.Exp+N.Exp+S.Exp+Entropy.Exp+Diverge.Exp+
             (Window_AmbigCodonFractn.Act*Reads.Act)+Window_ErrPerc.Act+Window_Subst.Act+Window_Entropy.Act+Window_Diverge.Act+
             Subst.Act+
             AmbigBase.Act+Entropy.Act+Diverge.Act+
             ErrPerc.Act+Err_N_Perc.Act+Err_S_Perc.Act,
           data=pred_Typical_N.Act_low, family=binomial(link="logit"))
summary(fit)

# Readd S.Act and remove Subst.Act.  Subst.Act is better AIC.
fit <- glm(Is_CrapLow_N.Act~dNdS.Exp+N.Exp+S.Exp+Entropy.Exp+Diverge.Exp+
             (Window_AmbigCodonFractn.Act*Reads.Act)+Window_ErrPerc.Act+Window_Subst.Act+Window_Entropy.Act+Window_Diverge.Act+
             S.Act+
             AmbigBase.Act+Entropy.Act+Diverge.Act+
             ErrPerc.Act+Err_N_Perc.Act+Err_S_Perc.Act,
           data=pred_Typical_N.Act_low, family=binomial(link="logit"))
summary(fit)

# Remove S.Act again.  Keep Subst.Act
fit <- glm(Is_CrapLow_N.Act~dNdS.Exp+N.Exp+S.Exp+Entropy.Exp+Diverge.Exp+
             (Window_AmbigCodonFractn.Act*Reads.Act)+Window_ErrPerc.Act+Window_Subst.Act+Window_Entropy.Act+Window_Diverge.Act+
             Subst.Act+
             AmbigBase.Act+Entropy.Act+Diverge.Act+
             ErrPerc.Act+Err_N_Perc.Act+Err_S_Perc.Act,
           data=pred_Typical_N.Act_low, family=binomial(link="logit"))
summary(fit)

# Remove Reads.Act
fit <- glm(Is_CrapLow_N.Act~dNdS.Exp+N.Exp+S.Exp+Entropy.Exp+Diverge.Exp+
             (Window_AmbigCodonFractn.Act:Reads.Act)+Window_AmbigCodonFractn.Act+Window_ErrPerc.Act+Window_Subst.Act+Window_Entropy.Act+Window_Diverge.Act+
             Subst.Act+
             AmbigBase.Act+Entropy.Act+Diverge.Act+
             ErrPerc.Act+Err_N_Perc.Act+Err_S_Perc.Act,
           data=pred_Typical_N.Act_low, family=binomial(link="logit"))
summary(fit)


# Remove Window_DIverge.Act
fit <- glm(Is_CrapLow_N.Act~dNdS.Exp+N.Exp+S.Exp+Entropy.Exp+Diverge.Exp+
             (Window_AmbigCodonFractn.Act:Reads.Act)+Window_AmbigCodonFractn.Act+Window_ErrPerc.Act+Window_Subst.Act+Window_Entropy.Act+
             Subst.Act+
             AmbigBase.Act+Entropy.Act+Diverge.Act+
             ErrPerc.Act+Err_N_Perc.Act+Err_S_Perc.Act,
           data=pred_Typical_N.Act_low, family=binomial(link="logit"))
summary(fit)

#' Remove ErrPerc.Act
fit <- glm(Is_CrapLow_N.Act~dNdS.Exp+N.Exp+S.Exp+Entropy.Exp+Diverge.Exp+
             (Window_AmbigCodonFractn.Act:Reads.Act)+Window_AmbigCodonFractn.Act+Window_ErrPerc.Act+Window_Subst.Act+Window_Entropy.Act+
             Subst.Act+
             AmbigBase.Act+Entropy.Act+Diverge.Act+
             Err_N_Perc.Act+Err_S_Perc.Act,
           data=pred_Typical_N.Act_low, family=binomial(link="logit"))
summary(fit)


#' Remove S.Exp
fit <- glm(Is_CrapLow_N.Act~dNdS.Exp+N.Exp+Entropy.Exp+Diverge.Exp+
             (Window_AmbigCodonFractn.Act:Reads.Act)+Window_AmbigCodonFractn.Act+Window_ErrPerc.Act+Window_Subst.Act+Window_Entropy.Act+
             Subst.Act+
             AmbigBase.Act+Entropy.Act+Diverge.Act+
             Err_N_Perc.Act+Err_S_Perc.Act,
           data=pred_Typical_N.Act_low, family=binomial(link="logit"))
summary(fit)


#' Remove Err_S_Perc.Act
fit <- glm(Is_CrapLow_N.Act~dNdS.Exp+N.Exp+Entropy.Exp+Diverge.Exp+
             (Window_AmbigCodonFractn.Act:Reads.Act)+Window_AmbigCodonFractn.Act+Window_ErrPerc.Act+Window_Subst.Act+Window_Entropy.Act+
             Subst.Act+
             AmbigBase.Act+Entropy.Act+Diverge.Act+
             Err_N_Perc.Act,
           data=pred_Typical_N.Act_low, family=binomial(link="logit"))
summary(fit)

#' Remove Err_N_Perc.Act
fit <- glm(Is_CrapLow_N.Act~dNdS.Exp+N.Exp+Entropy.Exp+Diverge.Exp+
             (Window_AmbigCodonFractn.Act:Reads.Act)+Window_AmbigCodonFractn.Act+Window_ErrPerc.Act+Window_Subst.Act+Window_Entropy.Act+
             Subst.Act+
             AmbigBase.Act+Entropy.Act+Diverge.Act,
           data=pred_Typical_N.Act_low, family=binomial(link="logit"))
summary(fit)
#' Best model so far.  
#' WTF????  Why does N.Act get crappier with lower Diverge.Exp?  but better with higher Entropy.Exp?
#' Plotting those variables shows that both have positive slopes

EDA(residuals(fit))  # too much data.  Use gamlssfit instead.

gamlssfit <- gamlss(Is_CrapLow_N.Act~dNdS.Exp+N.Exp+Entropy.Exp+Diverge.Exp+
                      (Window_AmbigCodonFractn.Act:Reads.Act)+Window_AmbigCodonFractn.Act+Window_ErrPerc.Act+Window_Subst.Act+Window_Entropy.Act+
                      Subst.Act+
                      AmbigBase.Act+Entropy.Act+Diverge.Act,
                    data=pred_Typical_N.Act_low, family=BI)
summary(gamlssfit)
plot(gamlssfit)



# Let's plot the probability
pred_Typical_N.Act_low$Is_CrapLow_N.Act_Prob <- fit$fitted.values
summary(pred_Typical_N.Act_low)
       
sapply(unlist(strsplit(names(fit$coefficients)[-1], ":")), # all coefficients except intercept.  Split up interactions
       function(colname) {
               fig <- ggplot(pred_Typical_N.Act_low, aes_string(x=colname, y="Is_CrapLow_N.Act_Prob")) + 
                 geom_smooth(method="glm", family="binomial") + 
                 geom_point(shape=1, alpha=0.5, aes(color=Is_CrapLow_N.Act)) + 
                 xlab(colname) + 
                 ylab("P(Crappy Low N.Act)") + 
                 ggtitle("Logistic regression for Crap Low N.Act for Typical Reads")
               print(fig)
             })




#' Compare this to stepwise AIC
stepfit <- glm(Is_CrapLow_N.Act~dNdS.Exp+N.Exp+S.Exp+Entropy.Exp+Diverge.Exp+
                 (Window_AmbigCodonFractn.Act*Reads.Act)+Window_ErrPerc.Act+Window_Subst.Act+Window_Entropy.Act+Window_Diverge.Act+
                 Subst.Act+
                 AmbigBase.Act+Entropy.Act+Diverge.Act+
                 ErrPerc.Act+Err_N_Perc.Act+Err_S_Perc.Act,
               data=pred_Typical_N.Act_low, family=binomial(link="logit"))
beststepfit <- stepAIC(stepfit, scope=list(lower=Is_CrapLow_N.Act~N.Exp), direction="both", trace=TRUE)
summary(beststepfit)
#' Resulting covariates are almost the same.  Except we discarded the main effects of Reads.Act and Window_AmbigCodonFractn.Act during manual



#' Compare this to stepwise AIC on everything
#' 
allstepfit <- glm(Is_CrapLow_N.Act~.,
               data=subset(pred_Typical_N.Act_low, 
                           select=-c(Source, Dist_dn_minus_dS, Dist_dNdS, Low_Dist_dNdS, Low_Dist_dn_minus_dS, LOD_dNdS, 
                                     LowLOD, Max.ErrFreeDownMask.N.Act, Min.ErrFreeDownMask.N.Act,
                                     Is_CrapLow_N.Act_Prob, CodonSite, Window_Start)), 
               family=binomial(link="logit"))
bestallstepfit <- stepAIC(allstepfit, scope=list(lower=Is_CrapLow_N.Act~N.Exp), direction="both", trace=TRUE)
summary(bestallstepfit)
plot(gamlss(bestallstepfit$formula, data=pred_Typical_N.Act_low, family=BI))
#' All this extra crap, even stuff that is barely insignificant at 0.0501...
#' Thinks that EN.*, ES.* matter.
#' Nice residuals though.

#' **For  N.Exp < 30, either the N.Act for Typical reads are good or they are way too high.  Why?**
#' 

pred_Typical_N.Act_hi <- merge(x=collate_exp_dnds_Ncheck, y=bound_errFreeDownMask, by="N.Exp", all=TRUE, sort=FALSE)
summary(pred_Typical_N.Act_hi)

pred_Typical_N.Act_hi$Is_CrapHi_N.Act <- as.factor(pred_Typical_N.Act_hi$Max.ErrFreeDownMask.N.Act < pred_Typical_N.Act_hi$N.Act)
pred_Typical_N.Act_hi <- subset(pred_Typical_N.Act_hi, subset=(Source == "Typical" & N.Act >= Min.ErrFreeDownMask.N.Act))
summary(pred_Typical_N.Act_hi)
dim(pred_Typical_N.Act_hi)

# Violin plots for the crappy high N.Act  
# Gap between true N.Exp and N.Act where it's too high widens when:
#' Entropy.Exp.  Low Expected Entropy = too crapy too high N.Act  
#' Diverge.Exp.  Low expected divergence = too crappy high N.Act  
#' dNdS.Exp.  Low expected dNdS = too crappy high N.Act  
#' Subst.Exp.  Low Subst = too crappy high N.Act  
#' N.Exp.  lower N.Exp = crappier higher N.Act 
#' Window_AmbigcodonFractn.Act = if we get less ambiguity in our window, we get too crappy high N.Act.  
#'      ??? Are we calculating window ambiguity properly here?  
#'      --> I think so, because less ambiguity can mean more errors.  Or more problematic sampling.
#' Window_ErrPerc.Act = if we get more error, then too crappy high N.Act
#' Window_Subst.Act = higher Actual substitutions in window (bigger window tree) = crappier higher N.Act
#' Window_Diverge.Act = lower divergence in window, crappier higher N.Act
#' Subst.Act = more subst per site, crappier higher N.Act
#' S.Act = more S actual per site, crappier higher N.Act
#' Diverge.Act = lower divergence per site, crappier higher N.Act
#' AmbigCodonFractn.Act = higher ambiguous codons, crappier higher N.Act
#' Err_N_Perc.Act = higher errors that lead to N, crappier Higer N.Act
#' ErrPerc.Act  - higher
#' AmbigBase.Act - higher
#' Entropy.Act - lower  ???? Why would lower randonness lead to higher N.Act than expected???
#' Diverge.Act - lower   ???Why would lower randonness lead to higher N.Act than expected???
sapply(COVAR_NAMES,
       function(colname) {
         fig <- ggplot(pred_Typical_N.Act_hi, aes_string(x="Is_CrapHi_N.Act", y=colname, color="Is_CrapHi_N.Act")) +
           geom_point(shape=1, alpha=0.5) + 
           geom_violin() +
           stat_summary(fun.ymin=function(y) {quantile(y, probs=c(0.25))},
                        fun.y=function(y){quantile(y, probs=c(0.5))},
                        fun.ymax=function(y) {quantile(y, probs=c(0.75))}, 
                        geom='crossbar') +
           xlab("Crappy and Hi N.Act?") + 
           ylab(colname) + 
           ggtitle(paste0("Violin plot Hi Crappy N.Act vs ", colname))
         print(fig)
       })

#' Stepwise AIC for explaining crappy hi N.Act
stephifit <- glm(Is_CrapHi_N.Act~(Reads.Act*Window_AmbigCodonFractn.Act)+Diverge.Act+Entropy.Act+AmbigBase.Act+
                   ErrPerc.Act+Err_N_Perc.Act+Err_S_Perc.Act+
                   AmbigCodonFractn.Act+Subst.Act+
                   Window_Diverge.Act+Window_Entropy.Act+Window_Subst.Act+
                   Window_ErrPerc.Act+Window_AmbigCodonFractn.Act+
                   N.Exp+S.Exp+Subst.Exp+EN.Exp+ES.Exp+dNdS.Exp+dN_minus_dS.Exp+
                   Diverge.Exp+Entropy.Exp, 
                 data=pred_Typical_N.Act_hi, family=binomial(link="logit"))
beststephifit <- stepAIC(stephifit, scope=list(lower=Is_CrapHi_N.Act~N.Exp), direction="both")
summary(beststephifit)

#' Manually remove insignificant var
#' Remove Window_ErrPerc.Act
beststephifit <- glm(Is_CrapHi_N.Act ~ Window_AmbigCodonFractn.Act + Diverge.Act + 
                       Entropy.Act + Err_N_Perc.Act + Err_S_Perc.Act + Subst.Act + 
                       Window_Subst.Act + N.Exp + S.Exp + EN.Exp + 
                       dNdS.Exp + dN_minus_dS.Exp + Diverge.Exp + Entropy.Exp,
                    data=pred_Typical_N.Act_hi, family=binomial )
summary(beststephifit)

#' Remove dN_minus_dS.Exp
beststephifit <- glm(Is_CrapHi_N.Act ~ Window_AmbigCodonFractn.Act + Diverge.Act + 
                       Entropy.Act + Err_N_Perc.Act + Err_S_Perc.Act + Subst.Act + 
                       Window_Subst.Act + N.Exp + S.Exp + EN.Exp + 
                       dNdS.Exp + Diverge.Exp + Entropy.Exp,
                     data=pred_Typical_N.Act_hi, family=binomial )
summary(beststephifit)


#' Remove Err_S_Perc.Act 
beststephifit <- glm(Is_CrapHi_N.Act ~ Window_AmbigCodonFractn.Act + Diverge.Act + 
                       Entropy.Act + Err_N_Perc.Act +  Subst.Act + 
                       Window_Subst.Act + N.Exp + S.Exp + EN.Exp + 
                       dNdS.Exp + Diverge.Exp + Entropy.Exp,
                     data=pred_Typical_N.Act_hi, family=binomial )
summary(beststephifit)

#' Remove Window_Subst.Act
beststephifit <- glm(Is_CrapHi_N.Act ~ Window_AmbigCodonFractn.Act + Diverge.Act + 
                       Entropy.Act + Err_N_Perc.Act +  Subst.Act + 
                       Window_Subst.Act + N.Exp + S.Exp + EN.Exp + 
                       dNdS.Exp + Diverge.Exp + Entropy.Exp,
                     data=pred_Typical_N.Act_hi, family=binomial )
summary(beststephifit)


#' Remove dN_minus_dS.Act
beststephifit <- glm(Is_CrapHi_N.Act ~ Window_AmbigCodonFractn.Act + Diverge.Act + 
                       Entropy.Act + ErrPerc.Act + Err_N_Perc.Act + Subst.Act + 
                       EN.Act + ES.Act +  
                       N.Exp + S.Exp + EN.Exp + ES.Exp + dNdS.Exp + 
                       Diverge.Exp + Entropy.Exp,
                     data=pred_Typical_N.Act_hi, family=binomial )
summary(beststephifit)
#' Best Model Fit where all variables significant

# Plot the probability fit
pred_Typical_N.Act_hi$Is_CrapHi_N.Act_Prob <- beststephifit$fitted.values
sapply(names(beststepfit$coefficients)[-1], function(colname) {
           fig <- ggplot(pred_Typical_N.Act_hi, aes_string(x=colname, y="Is_CrapHi_N.Act_Prob")) + 
             geom_smooth(method="glm", family="binomial") + 
             geom_point(shape=1, alpha=0.5, aes(color=Is_CrapHi_N.Act)) + 
             xlab(colname) + 
             ylab("P(Crappy Hi N.Act)") + 
             ggtitle("Logistic regression for Crap Hi N.Act for Typical Reads")
           print(fig)
         })


#' Violin plot inaccurate vs accurate dN/dS vs variables.
#' 
#' Gap between dNdS.Exp and dNdS.Act widens with:
#' Entropy.Exp - lower expected entropy
#' Diverge.Exp - lower expected divergence
#' dN_minus_dS.Exp - higher expected dN/dS
#' dNdS.Exp - higher expected dN/dS
#' Subst.Exp - lower expected substitutions
#' S.Exp - lower expected S
#' N.Exp - lower expected N
#' Window_AmbiguousCodonFractn.Act - higher ambiguous codons in window
#' Window_ErrPerc.Act - higher errors in window
#' Window_CodonDepth.Act - fewer codons in window
#' Window_Subst.Act - higher window substitutions
#' Window_Entropy.Act - higher window entropy
#' Window_Diverge.Act - higher window divergence
#' Subst.Act - fewer codon site subsitutions
#' S.Act - fewer codon site S
#' N.Act - more  codon site N
#' AmbigCodonFractn.Act - more ambiguous codons
#' AmgigBase.Act - more ambiguous bases
#' Entropy.Act - higher Entropy
#' Diverge.Act - higher diverge
#' CodonDepth.Act - lower codon depth
#' Reads.Act - fewer reads
#' 
sapply(COVAR_NAMES,
function(colname) {
  fig <- ggplot(collate_exp_dnds, aes_string(x="Low_Dist_dNdS", y=colname, color="Low_Dist_dNdS")) +    
    geom_violin() +
    stat_summary(fun.ymin=function(y) {quantile(y, probs=c(0.25))},
                 fun.y=function(y){quantile(y, probs=c(0.5))},
                 fun.ymax=function(y) {quantile(y, probs=c(0.75))}, 
                 geom='crossbar') +
    xlab("Accurate dN/dS?") + 
    ylab(colname) + 
    ggtitle(paste0("Violin plot Accurate vs Inaccurate dN/dS Variables ", colname))
  print(fig)
})

# # Create several training test sets with resampling
# trainRows <- createDataPartition(y=collate_exp_dnds$Dist_dNdS, times=5, list=FALSE, p=0.75)
# summary(trainRows)
# dim(trainRows)
# 
# # TODO:  Just do the first training test set for now.  Do them all in a loop
# trainSet <- collate_exp_dnds[trainRows[, 1],]
# dim(trainSet)
# summary(trainSet)
# 
# testSet <- collate_exp_dnds[-trainRows[, 1],]
# dim(testSet)
# 
# if (nrow(trainSet) + nrow(testSet) != nrow(collate_exp_dnds)) {
#   stop("Invalid training and test set")
# }
# 
# # Partial Least Squares Discriminant Analysis
# # Find the nubmer of Partial Least Squares to keep
# # I.e.  Find Least Number of orthogonal features that explain the most variance in response  (similar to PCA)
# plsFit <- train(form=Dist_dNdS~., data=trainSet, 
#                 method="pls",
#                 preProcess=c("center", "scale")  # centre and scale all features before fitting
#                 )
# 
# plsPredict <- predict(plsFit, newdata=testSet)
# summary(plsPredict)
# length(plsPredict)
# head(plsPredict)
# 
# 
# # How do Subst.Exp, Diverge.Act, dNdS.Exp affect distance between inferred-expected?
# lmfit <- glm(LOD_dNdS~Subst.Exp * Diverge.Act * dNdS.Exp, data=collate_exp_dnds)
# summary(lmfit)
# 
# # How do Diverge.Exp, Diverge.Act, dNdS.Exp affect distance between inferred-expected?
# lmfit <- glm(LOD_dNdS~Diverge.Exp * Diverge.Act * dNdS.Exp, data=collate_exp_dnds)
# summary(lmfit)
# 
# # How do Subst.Exp, Subst.Act, dNdS.Exp affect distance between inferred-expected?
# lmfit <- glm(LOD_dNdS~Subst.Exp * Subst.Act * dNdS.Exp, data=collate_exp_dnds)
# summary(lmfit)
# 
# # How do Subst.Exp, dNdS.Exp affect distance between inferred-expected?
# lmfit <- glm(LOD_dNdS~Subst.Exp * dNdS.Exp, data=collate_exp_dnds)
# summary(lmfit)
# 
# 
# # N.Exp, S.Exp
# lmfit <- glm(Dist_dNdS~N.Exp+S.Exp+
#                CodonSite+CodonDepth.Act+N.Act+S.Act+
#                Diverge.Exp+Entropy.Exp+Diverge.Act+Entropy.Act+
#                Window_CodonDepth.Act+Window_Diverge.Act+Window_Entropy.Act, 
#              data=collate_exp_dnds)
# bestfit <- stepAIC(lmfit)
# 
# lmfit <- glm(Dist_dNdS~N.Exp, data=collate_exp_dnds)
# summary(lmfit)
# 
# lmfit <- glm(Dist_dNdS~S.Exp, data=collate_exp_dnds)
# summary(lmfit)
# 
# 
# lmfit <- glm(Dist_dNdS~N.Exp+S.Exp, data=collate_exp_dnds)
# summary(lmfit)
# 
# # N.Exp, S.Exp
# lmfit <- glm(Dist_dNdS~N.Exp+S.Exp+Subst.Exp, data=collate_exp_dnds)
# summary(lmfit)
# 
# # CodonDepth, S.Act (most)
# lmfit <- glm(Dist_dNdS~CodonSite+CodonDepth.Act+N.Act+S.Act+Subst.Act, data=collate_exp_dnds)
# summary(lmfit)
# 
# collate_exp_dnds$Diverge.ExpPerc <- collate_exp_dnds$Diverge.Exp*100
# collate_exp_dnds$Entropy.ExpPerc <- collate_exp_dnds$Entropy.Exp*100
# lmfit <- glm(Dist_dNdS~Diverge.ExpPerc+Entropy.ExpPerc, data=collate_exp_dnds)
# summary(lmfit)
# 
# # Scaling the values by a straight factor does not change the p-values
# lmfit <- glm(Dist_dNdS~Diverge.Exp+Entropy.Exp, data=collate_exp_dnds)
# summary(lmfit)
# 
# lmfit <- glm(Dist_dNdS~Diverge.Act+Entropy.Act, data=collate_exp_dnds)
# summary(lmfit)
# 
# collate_exp_dnds$Entropy.ActLog <- log10(collate_exp_dnds$Entropy.Act)
# collate_exp_dnds$Entropy.ActLog[collate_exp_dnds$Entropy.Act==0] <- NA
# lmfit <- glm(Dist_dNdS~Entropy.ActLog, data=collate_exp_dnds)
# summary(lmfit)
# 
# 
# lmfit <- glm(Dist_dNdS~Entropy.Act, data=collate_exp_dnds)
# summary(lmfit)
# 
# 
# lmfit <- glm(Dist_dNdS~Window_CodonDepth.Act+Window_Diverge.Act+Window_Entropy.Act, data=collate_exp_dnds)
# summary(lmfit)
# 
# lmfit <- glm(Dist_dNdS~Window_CodonDepth.Act, data=collate_exp_dnds)
# summary(lmfit)
# 
# # When you include everything, only S.Exp, CodonSite, CodonDepth.Act, S.Act are sig
# # N.act no longer sig
# lmfit <- glm(Dist_dNdS~N.Exp+S.Exp+Subst.Exp+
#                CodonSite+CodonDepth.Act+N.Act+S.Act+Subst.Act+
#                Diverge.Exp+Entropy.Exp+Diverge.Act+Entropy.Act+
#                Window_CodonDepth.Act+Window_Diverge.Act+Window_Entropy.Act, 
#              data=collate_exp_dnds)
# summary(lmfit)
# 
# # Only include any variable that was significant in a previous fit.  they are all significant again.
# lmfit <- glm(Dist_dNdS~S.Exp+CodonSite+CodonDepth.Act+N.Act+S.Act,
#              data=collate_exp_dnds)
# summary(lmfit)
# 
# # Don't Keep N.Act
# anova(glm(Dist_dNdS~S.Exp+CodonSite+CodonDepth.Act+S.Act,data=collate_exp_dnds),
#       glm(Dist_dNdS~S.Exp+CodonSite+CodonDepth.Act+N.Act+S.Act,data=collate_exp_dnds),
#       test="Chisq")
# 
# 
# # Keep S.exp
# anova(glm(Dist_dNdS~CodonSite+CodonDepth.Act+S.Act,data=collate_exp_dnds),
#       glm(Dist_dNdS~S.Exp+CodonSite+CodonDepth.Act+S.Act,data=collate_exp_dnds),
#       test="Chisq")
# 
# # Interaction is better
# anova(glm(Dist_dNdS~S.Exp*CodonSite*CodonDepth.Act*S.Act,data=collate_exp_dnds),
#       glm(Dist_dNdS~S.Exp+CodonSite+CodonDepth.Act+S.Act,data=collate_exp_dnds),
#       test="Chisq")
# 
# anova(glm(Dist_dNdS~S.Exp+CodonSite+CodonDepth.Act+S.Act,data=collate_exp_dnds),
#       glm(Dist_dNdS~S.Exp*CodonSite*CodonDepth.Act*S.Act,data=collate_exp_dnds),
#       test="Chisq")
# 
# 
# # But the pvalues for interaction say only codonsite:codonDepth is significant
# lmfit <- glm(Dist_dNdS~S.Exp*CodonSite*CodonDepth.Act*S.Act,data=collate_exp_dnds)
# summary(lmfit)
# 
# anova(glm(Dist_dNdS~S.Exp+(CodonSite*CodonDepth.Act)+S.Act,data=collate_exp_dnds),
#       glm(Dist_dNdS~S.Exp*CodonSite*CodonDepth.Act*S.Act,data=collate_exp_dnds),
#       test="Chisq")
# 
# anova(glm(Dist_dNdS~CodonSite*CodonDepth.Act,data=collate_exp_dnds),
#       glm(Dist_dNdS~S.Exp*CodonSite*CodonDepth.Act*S.Act,data=collate_exp_dnds),
#       test="Chisq")
# 
# 
# # Plot Covariates
# ggplot(collate_exp_dnds, aes(x=CodonSite, y=log10(Dist_dNdS+0.000000000000000000000000000000000001))) + 
#   geom_point(shape=1, alpha=0.1) + 
#   geom_smooth(method="lm") + 
#   xlab("\nCodon Site") + 
#   ylab("Actual Window-Codon Site dN/dS - Expected Site dN/dS") + 
#   ggtitle("Plot Covariates Against Difference Between Actual Window-Site dN/dS and Expected Site dN/dS")
# 
# 
# ggplot(collate_exp_dnds, aes(x=CodonDepth.Act, y=Dist_dNdS)) + 
#   geom_point(shape=1, alpha=0.1) + 
#   geom_smooth(method="lm") + 
#   xlab("\nActual Codon Depth ") + 
#   ylab("Actual Window-Codon Site dN/dS - Expected Site dN/dS") + 
#   ggtitle("Plot Covariates Against Difference Between Actual Window-Site dN/dS and Expected Site dN/dS")
# 
# ggplot(collate_exp_dnds, aes(x=CodonDepth.Act, y=log10(Dist_dNdS+0.0001))) + 
#   geom_point(shape=1, alpha=0.1) + 
#   geom_smooth(method="lm") + 
#   xlab("\nActual Codon Depth ") + 
#   ylab("Actual Window-Codon Site dN/dS - Expected Site dN/dS") + 
#   ggtitle("Plot Covariates Against Difference Between Actual Window-Site dN/dS and Expected Site dN/dS")
# 
# 
# ggplot(collate_exp_dnds, aes(x=S.Exp, y=log10(Dist_dNdS+0.0001))) + 
#   geom_point(shape=1, alpha=0.1) + 
#   geom_smooth(method="lm") + 
#   xlab("\nExpected Synonymous Substitutions ") + 
#   ylab("Actual Window-Codon Site dN/dS - Expected Site dN/dS") + 
#   ggtitle("Plot Covariates Against Difference Between Actual Window-Site dN/dS and Expected Site dN/dS")
# 
# ggplot(collate_exp_dnds, aes(x=S.Exp, y=Dist_dNdS)) + 
#   geom_point(shape=1, alpha=0.1) + 
#   geom_smooth(method="lm") + 
#   xlab("\nExpected Synonymous Substitutions ") + 
#   ylab("Actual Window-Codon Site dN/dS - Expected Site dN/dS") + 
#   ggtitle("Plot Covariates Against Difference Between Actual Window-Site dN/dS and Expected Site dN/dS")
