# For multiple simulations, find the variables that affect accuracy

library(knitr)
opts_chunk$set(warning=FALSE, message=FALSE, width=1200, echo=TRUE)


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
library(speedglm)  # faster glms in parallel
library(caret)  # for finding correlations
library(gtools)  # for inv.logit
library(arm)  # for binned.plot for logistic regresison residuals
source ('./speedStepAIC.R')


DNDS_FILENAME <- "/home/thuy/gitrepo/Umberjack_Benchmark/simulations/out/collate_all.csv"

dnds <- read.table(DNDS_FILENAME, header=TRUE, sep=",", na.strings=c("", "None"))
dim(dnds)
summary(dnds)
head(dnds)
dnds$ErrBaseRate.Act <- dnds$ErrBase.Act/dnds$Reads.Act  # per-window-codonsite-read nucleotide errors
dnds$AmbigPadBaseRate.Act <- dnds$AmbigPadBase/dnds$Reads.Act  # per-window-codonsite-read nucleotide N's or gaps
dnds$UnambigCodonRate.Act <- dnds$UnambigCodons.Act/ dnds$Reads.Act  # per-window-codonsite fraction of unambiguous codons
dnds$Subst.Act <- dnds$N.Act + dnds$S.Act
dnds$Subst.Exp <- dnds$N.Exp + dnds$S.Exp
dnds$Cov.Act <- dnds$Reads.Act/dnds$PopSize.Act
dnds <- subset(dnds, select=-c(PopSize.Act, ErrBase.Act, AmbigPadBase.Act, UnambigCodons.Act))
dim(dnds)
summary(dnds)



# Average across all codon sites in a window
per_window_ave <- ddply(.data=dnds, .variables=c("File", "Window_Start"), 
                        .fun=function(x) {                            
                          data.frame(Window_Conserve.Act=mean(x$ConserveTrueBase.Act, na.rm=TRUE),
                                     Window_Entropy.Act=mean(x$EntropyTrueBase.Act, na.rm=TRUE),
                                     Window_UnambigCodonRate.Act=mean(x$UnambigCodonRate.Act, na.rm=TRUE),
                                     Window_ErrBaseRate.Act=mean(x$ErrBaseRate.Act, na.rm=TRUE))
                        })
dnds <- merge(x=dnds, y=per_window_ave, by=c("File", "Window_Start"), all=TRUE, sort=TRUE)
dim(dnds)
summary(dnds)



# Now remove window-codonsites where there is no dnds and no dn-ds information because of insufficient window sequences
dnds <- dnds[!is.na(dnds$dN_minus_dS.Exp) & !is.na(dnds$N.Act), ]
dim(dnds)
summary(dnds)


# Ignore dN/dS == 0 for numerical stability
dnds$LOD_dNdS <- log(dnds$dNdS.Act) - log(dnds$dNdS.Exp)
dnds$LOD_dNdS[dnds$dNdS.Exp == 0 | dnds$dNdS.Act == 0] <- NA
# dN_minus_dS actually refers to dN - dS / tree length.
#  Thus dN_minus_dS == 1 means that the difference between nonsyn substitution rate and syn substitution rate 
#  is 100% of the total substitutions in the tree.
dnds$Dist_dn_minus_dS <- dnds$dN_minus_dS.Act - dnds$dN_minus_dS.Exp
dnds$CrapLOD <- as.factor(abs(dnds$LOD_dNdS) > 1)
dnds$CrapDist <- as.factor(abs(dnds$Dist_dn_minus_dS) > 0.5)
summary(dnds)
dim(dnds)
head(dnds)

NUM_RESP_NAMES <- c("LOD_dNdS", "Dist_dn_minus_dS")
CAT_RESP_NAMES <- c("CrapLOD", "CrapDist")
COVAR_NAMES <- colnames(dnds[sapply(dnds,is.numeric)])[!colnames(dnds[sapply(dnds,is.numeric)]) %in% NUM_RESP_NAMES]
LM_COVAR_NAMES <- COVAR_NAMES[!(COVAR_NAMES %in% c("dNdS.Act", "dN_minus_dS.Act", "dN_minus_dS.Exp", 
                                                   # In separate analysis, conservation and entropy are highly correlated.
                                                   # When we use speedglm, it bugs out after it removes highly correlated variables.
                                                   # So we do it for them.
                                                   "ConserveTrueBase.Act", "ConserveTrueBase.Exp", "Window_Conserve.Act",
                                                   "UnambigCodonRate.Act", "Window_UnambigCodonRate.Act",
                                                   # These are highly correlated with N, S
                                                   "Subst.Act", "Subst.Exp",
                                                   "Window_Start", "Window_End", "CodonSite", "Reads.Act"
                                                   )
                                )]

#' Error Rate:
#' =============================================
#' 
#' **Nucleotide Error Rate After Umberjack Quality Masking = `r mean(dnds$ErrBaseRate.Act)/3 `**
#' 





#' Correlation Between Features
#' =============================================
#' 

# Find correlation between features.  Don't bother scaling.  The correlation heatmap is the same whether we scale & centre or not.
corMat <- cor(dnds[, c(NUM_RESP_NAMES, COVAR_NAMES)], use="complete.obs", method="spearman")
#+ fig.width=15, fig.height=15
heatmap.2(corMat,  col=bluered, density.info="none", trace="none", srtCol=45, main="Feature Correlation", margins=c(12,12), 
          colsep=(1:ncol(corMat)), rowsep=(1:nrow(corMat)), sepwidth=c(0.05, 0.05), sepcolor="black")





#' Density Plots of Each Numerical Column
#' =============================================
#' 
#' 
plot_density <- function(colname) {
  fig <- ggplot(dnds, aes_string(x=colname)) + 
    geom_density(na.rm=TRUE, color="black") + 
    geom_density(aes_string(x=colname), na.rm=TRUE) + 
    ggtitle(paste0("Density plot ", colname))
  print(fig)
}

sapresults <- sapply(COVAR_NAMES, plot_density)


#' Plot Category Response Vs Variable
#' =============================================
#' 
#' CrapLOD means a Log Odds Difference Between estimated dN/dS and expected dN/dS > 1.  
#' 
#' CrapDist means the distance between estimated scaled dN-dS and expected scaled dN-dS > 0.5.
#' 
#' 


plot_catresp_vs_var <- function(resp_colname) {
  sapresults <- sapply(COVAR_NAMES, 
         function(var_colname) {
           fig <- ggplot(dnds[!is.na(dnds[, resp_colname]),], aes_string(x=resp_colname, y=var_colname, color=resp_colname)) +             
             geom_boxplot() + 
             ggtitle(paste0("Box Plot ", resp_colname, " vs ", var_colname))
           print(fig)
         })
}
plot_catresp_vs_var("CrapLOD")
plot_catresp_vs_var("CrapDist")



#' What variables determine really bad accuracy?
#' =============================================
#' 


# # Only works for gamlss
# fit_and_plot_bin<- function(resp_colname, fit, df) {
#   bestfit <- stepAIC(fit, direction="both", trace=TRUE)
#   summ <- summary(bestfit)
#   plot(bestfit)
#   
#   coeffs <- data.frame(summ)
#   colnames(coeffs)[grep("Pr", colnames(coeffs))] <- "pval"
#   coeffs$adj.pval <- p.adjust(coeffs$pval, method="BH")
#   coeffs$name <- rownames(coeffs)
#   head(coeffs)
#   summary(coeffs)
#   
#   predictors <- coeffs$name[coeffs$name != "(Intercept)"]
#   df$fitted.values <- bestfit$mu.fv
#   figs <- sapply(predictors, 
#                  function(varname) {
#                    fig <- ggplot(df, aes_string(x=varname, y="fitted.values", color=resp_colname)) + 
#                      geom_point(alpha=0.5, shape=1) + 
#                      geom_smooth(method="glm", family="binomial", color="black") + 
#                      xlab(varname) + 
#                      ylab(paste0("Fitted Values of ", resp_colname, "\n")) + 
#                      ggtitle(paste0("Fitted Values of ", resp_colname, " Vs ", varname, " qval=", 
#                                     signif(coeffs[coeffs$name==varname, "adj.pval"], 2)))
#                    print(fig)
#                  })
#   return (bestfit)
#   
#  
# }


#' speedGLM  fit and plot
#' 
fit_and_plot_bin_fast<- function(resp_colname, fit, df) {
  
  bestfit <- stepAIC(fit, direction="both", trace=TRUE)
  print(summary(bestfit))
  
  # speedglm doesn't expose residuals or fitted values. Do it ourselves
  df_fit <- data.frame(Intercept=1, df[, colnames(df)[colnames(df) %in% names(coef(bestfit))]])
  df_fit$fitted.values <- inv.logit(as.vector(as.matrix(df_fit) %*% coef(bestfit)))  # [0, 1]
  df_fit[, resp_colname] <-  df[, resp_colname]  
  df_fit$NumResp <- ifelse(df_fit[, resp_colname] == TRUE, 1, 0)  #[0, 1]
  df_fit$residuals <- df_fit$NumResp-df_fit$fitted.values   #[-1, 1]
  summary(df_fit)
  
  binnedplot(x=df_fit$fitted.values, y=df_fit$residuals, nclass=NULL, 
             xlab="Predicted Values", ylab="Average residual", 
             main="Binned residual plot", 
             cex.pts=0.8, col.pts=1, col.int="gray")
  
  coeffs <- data.frame(summary(bestfit)$coefficients)
  colnames(coeffs)[grep("Pr", colnames(coeffs))] <- "pval"
  coeffs$pval <- as.numeric(as.character(coeffs$pval))
  coeffs$adj.pval <- p.adjust(coeffs$pval, method="BH")
  coeffs$name <- rownames(coeffs)
  head(coeffs)
  summary(coeffs)
  str(coeffs)
  
  # Print confusion matrix
  print(confusionMatrix(data=ifelse(df_fit$fitted.values>=0.5, TRUE, FALSE), reference=df_fit[, resp_colname]))
  predictors <- coeffs$name[coeffs$name != "(Intercept)"]
  
  figs <- sapply(predictors, 
                 function(varname) {
                   fig <- ggplot(df_fit, aes_string(x=varname, y="fitted.values", color=resp_colname)) + 
                     geom_point(alpha=0.5, shape=1) + 
                     geom_smooth(method="glm", family="binomial", color="black") + 
                     xlab(varname) + 
                     ylab(paste0("Predicted Values of ", resp_colname, "\n")) + 
                     ggtitle(paste0("Predicted Values of ", resp_colname, " Vs ", varname, " qval=", 
                                    signif(coeffs[coeffs$name==varname, "adj.pval"], 2)))
                   print(fig)
                 })
  
  return (bestfit)
}

cleandnds <- na.omit(dnds)
summary(cleandnds)
dim(cleandnds)



#' **Stepwise AIC to find what causes extremely inaccurate estimates of dN/dS**
crapLODFormula <- as.formula(paste0("CrapLOD~", paste0(LM_COVAR_NAMES, collapse=" + ")))
print(crapLODFormula)
allstepfitLOD <- speedglm(crapLODFormula, data=cleandnds, family=binomial(link='logit'))
bestallstepfitLOD <- fit_and_plot_bin_fast("CrapLOD", allstepfitLOD, df=cleandnds)


#' **Stepwise AIC to find what causes extremely inaccurate estimates of scaled dn-ds**
crapDistFormula <- as.formula(paste0("CrapDist~", paste0(LM_COVAR_NAMES, collapse=" + ")))
print(crapDistFormula)
allstepfitDist <- speedglm(crapDistFormula, data=cleandnds, family=binomial(link='logit'))
bestallstepfitDist <- fit_and_plot_bin_fast("CrapDist", allstepfitDist, df=cleandnds)



#' Plot Numerical Response Vs Numerical Variables
#' ==================================================
plot_resp_vs_var <- function(resp_colname) {
  figs <- sapply(COVAR_NAMES, 
         function(var_colname) {
           fig <- ggplot(dnds, aes_string(x=var_colname, y=resp_colname)) +             
             geom_point(shape=1, alpha=0.5, na.rm=TRUE) + 
             geom_smooth(method="lm") + 
             ggtitle(paste0(resp_colname, " vs ", var_colname))
           print(fig)
         })
}
plot_resp_vs_var("Dist_dn_minus_dS")
plot_resp_vs_var("LOD_dNdS")


#' Concordance
#' -----------------------------
#' 

#' **Find Concordance with all original data**
#' 
concord <- epi.ccc(dnds[!is.na(dnds$dNdS.Act) & !is.na(dnds$dNdS.Exp), ]$dNdS.Act, 
                   dnds[!is.na(dnds$dNdS.Act) & !is.na(dnds$dNdS.Exp), ]$dNdS.Exp)
print(concord$rho.c$est)
#' **Concordance for dN/dS when all considered = `r concord$rho.c$est`**
#' 

concord <- epi.ccc(dnds[!is.na(dnds$dN_minus_dS.Act) & !is.na(dnds$dN_minus_dS.Exp), ]$dN_minus_dS.Act, 
                   dnds[!is.na(dnds$dN_minus_dS.Act) & !is.na(dnds$dN_minus_dS.Exp), ]$dN_minus_dS.Exp)
print(concord$rho.c$est)
#' **Concordance for dN-dS when all considered = `r concord$rho.c$est`**
#' 

#' **Now Remove Window-Sites with Nonsynonymous or Synonymous Substitutions higher than zero but lower than 1.**  
#' 
#'  When sites have substitutions higher than zero but lower than 1, that means that there were no substitutions
#'  of that kind other than the ones transition to/from ambiguous codons.
#' 
dnds <- dnds[!(dnds$S.Act < 1 & dnds$S.Act > 0) & !(dnds$N.Act < 1 & dnds$N.Act > 0), ]
summary(dnds)
dim(dnds)

concord <- epi.ccc(dnds[!is.na(dnds$dNdS.Act) & !is.na(dnds$dNdS.Exp), ]$dNdS.Act, 
                   dnds[!is.na(dnds$dNdS.Act) & !is.na(dnds$dNdS.Exp), ]$dNdS.Exp)
print(concord$rho.c$est)
#' **Concordance for dN/dS when only good window-sites considered = `r concord$rho.c$est`**
#' 

concord <- epi.ccc(dnds[!is.na(dnds$dN_minus_dS.Act) & !is.na(dnds$dN_minus_dS.Exp), ]$dN_minus_dS.Act, 
                   dnds[!is.na(dnds$dN_minus_dS.Act) & !is.na(dnds$dN_minus_dS.Exp), ]$dN_minus_dS.Exp)
print(concord$rho.c$est)
#' **Concordance for dN-dS when only good window-sites considered = `r concord$rho.c$est`**
#' 

#' Plot actual versus expected
fig <- ggplot(dnds, aes(x=dNdS.Exp, y=dNdS.Act)) + 
  geom_point(alpha=0.5, shape=1) + 
  geom_smooth(method="lm", color="Red") + 
  xlab("\nExpected dN/dS") + 
  ylab("Inferred dN/dS\n") + 
  ggtitle("Scatterplot Expected vs Inferred dN/dS, Excl Low Syn")
print(fig)

fig <- ggplot(dnds, aes(x=dN_minus_dS.Exp, y=dN_minus_dS.Act)) + 
  geom_point(alpha=0.5, shape=1) + 
  geom_smooth(method="lm", color="Red") + 
  xlab("\nExpected dN/dS") + 
  ylab("Inferred dN/dS\n") + 
  ggtitle("Scatterplot Expected vs Inferred dN-dS, Exclude low Syn")
print(fig)


#' 
#' GLM for Predictors for Bad performance, Ignore Sites with N or S higher than zero but lower than 1.  
#' ========================================================
#' 

# Expects continuous response for speedglm
fit_and_plot_glm_fast <- function(resp_colname, fit, df) {
  
  bestfit <- stepAIC(fit, direction="both", trace=TRUE)
  print(summary(bestfit))
  
  # speedglm doesn't expose residuals or fitted values. Do it ourselves
  df_fit <- data.frame(Intercept=1, df[, colnames(df)[colnames(df) %in% names(coef(bestfit))]])
  df_fit$fitted.values <- as.vector(as.matrix(df_fit) %*% coef(bestfit))
  df_fit[, resp_colname] <-  df[, resp_colname]    
  df_fit$residuals <-  df_fit[, resp_colname] -df_fit$fitted.values
  summary(df_fit)
  
  fig <- ggplot(df_fit, aes(x=fitted.values, y=residuals)) + 
    geom_point(alpha=0.5, shape=1) + 
    geom_smooth(method="loess", color="Red") +
    xlab("\nPredicted Values") + 
    ylab("Residuals\n") + 
    ggtitle("Plot Predicted vs Residuals")
  print(fig)
  
  # residual normality
  qqnorm(df_fit$residuals)
  qqline(df_fit$residuals)
  
  coeffs <- data.frame(summary(bestfit)$coefficients)
  colnames(coeffs)[grep("Pr", colnames(coeffs))] <- "pval"
  coeffs$pval <- as.numeric(as.character(coeffs$pval))
  coeffs$adj.pval <- p.adjust(coeffs$pval, method="BH")
  coeffs$name <- rownames(coeffs)
  head(coeffs)
  summary(coeffs)
  str(coeffs)
  
  predictors <- coeffs$name[coeffs$name != "(Intercept)"]
  
  figs <- sapply(predictors, 
                 function(varname) {
                   fig <- ggplot(df_fit, aes_string(x=varname, y="residuals")) + 
                     geom_point(alpha=0.5, shape=1) + 
                     geom_smooth(method="lm", color="Red") + 
                     xlab(varname) + 
                     ylab("Residuals\n") + 
                     ggtitle(paste0("Residuals of  ", resp_colname, " Vs ", varname, " qval=", 
                                    signif(coeffs[coeffs$name==varname, "adj.pval"], 2)))
                   print(fig)
                 })
  
  return (bestfit)
}

# #' Stepwise AIC to find causes any difference between estimated and expected dn/ds
# fit_and_plot_glm <- function(resp_colname, fit) {
#   
#   
#   bestfit <- stepAIC(fit, direction="both", trace=TRUE)
#   print(summary(bestfit))
#   plot(bestfit)
#   
#   coeffs <- data.fcoeffs$adj.pval <- p.adjust(coeffs$Pr...t.., method="BH")
#   coeffs$name <- rownames(coeffs)
#   head(coeffs)
#   summary(coeffs)
#   
#   predictors <- coeffs$name[coeffs$name != "(Intercept)"]
#   
#   figs <- sapply(predictors, 
#          function(varname) {
#            fig <- ggplot(na.omit(dnds), aes_string(x=varname, y=resp_colname)) +             
#              geom_point(shape=1, alpha=0.5, na.rm=TRUE) + 
#              geom_smooth(method="lm", color="Red") + 
#              ggtitle(paste0(resp_colname, " vs ", varname,
#                             " qval=", signif(coeffs[coeffs$name==varname, "adj.pval"], 2)))
#            print(fig)
#          })
#   return (bestfit)
# }

cleandnds <- na.omit(dnds)
summary(cleandnds)
dim(cleandnds)

LODFormula <- as.formula(paste0("LOD_dNdS~", paste0(LM_COVAR_NAMES, collapse=" + ")))
print(LODFormula)
allfitLOD <- speedglm(LODFormula, data=cleandnds)
bestfit <- fit_and_plot_glm_fast("LOD_dNdS", allfitLOD, df=cleandnds)



DistFormula <- as.formula(paste0("Dist_dn_minus_dS~", paste0(LM_COVAR_NAMES, collapse=" + ")))
print(DistFormula)
allfitDist <- speedglm(DistFormula, data=cleandnds)
bestfit <- fit_and_plot_glm_fast("LOD_dNdS", allfitDist, df=cleandnds)

