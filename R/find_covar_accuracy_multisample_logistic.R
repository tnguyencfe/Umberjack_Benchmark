# For multiple simulations, find the variables that affect accuracy

library(knitr)
opts_chunk$set(warning=FALSE, message=FALSE, width=1200, echo=TRUE)

#+ message=FALSE
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
source('./load_all_sim_dnds.R')

dnds <- get_all_sim_dnds()
dim(dnds)
summary(dnds)

NUM_RESP_NAMES <- c("LOD_dNdS", "Dist_dn_minus_dS", "AbsLOD_dNdS", "AbsDist_dn_minus_dS", "SqDist_dn_minus_dS")
CAT_RESP_NAMES <- c("CrapLOD", "CrapDist")
NUM_NAMES <- colnames(dnds[sapply(dnds,is.numeric)])
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
                             "EN.Exp", "ES.Exp", "EN.Act", "ES.Act",
                             "Window_Start", "Window_End", "CodonSite", "Reads.Act", "PopSize.Act", "Is_Break"
                           )
                         ]
CAT_COVAR_NAMES <- c("IsLowSubst.Act")
LM_COVAR_NAMES <- c(CAT_COVAR_NAMES, COVAR_NAMES)

#' Error Rate:
#' =============================================
#' 
#' **Nucleotide Error Rate After Umberjack Quality Masking = `r mean(dnds$ErrPerCodon.Act)/3 `**
#' 





#' Correlation Between Features
#' =============================================
#' 

# Find correlation between features.  Don't bother scaling.  The correlation heatmap is the same whether we scale & centre or not.
corMat <- cor(dnds[, c(NUM_RESP_NAMES, COVAR_NAMES)], use="complete.obs", method="spearman")
#+ fig.width=15, fig.height=15
heatmap.2(corMat,  col=bluered, density.info="none", trace="none", srtCol=45, main="Feature Correlation", margins=c(12,12), 
          colsep=(1:ncol(corMat)), rowsep=(1:nrow(corMat)), sepwidth=c(0.05, 0.05), sepcolor="black")





#' Density Plots of Each Numerical Variable
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

#+ fig.width=5, fig.height=5
sapresults <- sapply(COVAR_NAMES, plot_density)


#' Plot Category Response Vs Variable
#' =============================================
#' 
#' CrapLOD means a Log Odds Difference Between estimated dN/dS and expected dN/dS > 1.  
#' 
#' CrapDist means the distance between estimated |scaled dN-dS and expected scaled dN-dS| > 1.
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

plot_catresp_vs_catvar <- function(resp_colname) {
  sapresults <- sapply(CAT_COVAR_NAMES, 
                       function(var_colname) {
                         fig <- ggplot(dnds[!is.na(dnds[, resp_colname]),], aes_string(x=resp_colname, fill=var_colname)) +             
                           geom_bar(color="Black") + 
                           ggtitle(paste0("Stacked BarPlot ", resp_colname, " vs ", var_colname))
                         print(fig)
                       })
}
plot_catresp_vs_catvar("CrapLOD")
plot_catresp_vs_catvar("CrapDist")

#' What variables determine really bad accuracy?
#' =============================================
#' 

#' speedGLM  fit and plot
#' 
fit_and_plot_bin_fast<- function(resp_colname, fit, df) {
  
  bestfit <- stepAIC(fit, direction="both", trace=FALSE)
  print(summary(bestfit))
  
  # speedglm doesn't expose residuals or fitted values. Do it ourselves
  df_fit <- data.frame(Intercept=1, df[, attributes(bestfit$terms)$term.labels])
  
  # TODO:  hack - this is temp hack to get factors to be numeric
  if (length(grep("IsLowSubst.Act", colnames(df_fit)) > 0)){
    df_fit$IsLowSubst.Act <- ifelse (df_fit$IsLowSubst.Act == TRUE, 1, 0)  
  }
  
  
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
  print(coeffs)
  
  
  # Print confusion matrix
  print("Confusion Matrix")
  print(confusionMatrix(data=ifelse(df_fit$fitted.values>=0.5, TRUE, FALSE), reference=df_fit[, resp_colname]))
  
  figs <- sapply(attributes(bestfit$terms)$term.labels, 
                 function(varname) {
                   fig <- ggplot(df_fit, aes_string(x=varname, y="fitted.values", color=resp_colname)) + 
                     geom_point(alpha=0.5, shape=1) + 
                     geom_smooth(method="glm", family="binomial", color="black") + 
                     xlab(varname) + 
                     ylab(paste0("Predicted Values of ", resp_colname, "\n")) + 
                     ggtitle(paste0("Predicted Values of ", resp_colname, " Vs ", varname))
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
print(summary(allstepfitLOD))
bestallstepfitLOD <- fit_and_plot_bin_fast(resp_colname="CrapLOD", fit=allstepfitLOD, df=cleandnds)


#' **Stepwise AIC to find what causes extremely inaccurate estimates of scaled dn-ds**
crapDistFormula <- as.formula(paste0("CrapDist~", paste0(LM_COVAR_NAMES, collapse=" + ")))
print(crapDistFormula)
allstepfitDist <- speedglm(crapDistFormula, data=cleandnds, family=binomial(link='logit'))
print(summary(allstepfitDist))
bestallstepfitDist <- fit_and_plot_bin_fast("CrapDist", allstepfitDist, df=cleandnds)



#' Plot Numerical Response Vs Numerical Variables
#' ==================================================
plot_resp_vs_var <- function(resp_colname) {
  figs <- sapply(COVAR_NAMES, 
                 function(var_colname) {
                   fig <- ggplot(dnds[!is.na(dnds[, resp_colname]) & abs(dnds[, resp_colname]) < 20,], aes_string(x=var_colname, y=resp_colname)) +             
                     geom_point(shape=1, alpha=0.5, na.rm=TRUE) + 
                     geom_smooth(method="lm") + 
                     ggtitle(paste0(resp_colname, " vs ", var_colname))
                   print(fig)
                 })
}
plot_resp_vs_var("Dist_dn_minus_dS")
plot_resp_vs_var("LOD_dNdS")


#' Concordance
#' ===========================================
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
gooddnds <- dnds[!(dnds$S.Act < 1 & dnds$S.Act > 0) & !(dnds$N.Act < 1 & dnds$N.Act > 0), ]
summary(gooddnds)
dim(gooddnds)

concord <- epi.ccc(gooddnds[!is.na(gooddnds$dNdS.Act) & !is.na(gooddnds$dNdS.Exp), ]$dNdS.Act, 
                   gooddnds[!is.na(gooddnds$dNdS.Act) & !is.na(gooddnds$dNdS.Exp), ]$dNdS.Exp)
print(concord$rho.c$est)
#' **Concordance for dN/dS when only good window-sites considered = `r concord$rho.c$est`**
#' 

concord <- epi.ccc(gooddnds[!is.na(gooddnds$dN_minus_dS.Act) & !is.na(gooddnds$dN_minus_dS.Exp), ]$dN_minus_dS.Act, 
                   gooddnds[!is.na(gooddnds$dN_minus_dS.Act) & !is.na(gooddnds$dN_minus_dS.Exp), ]$dN_minus_dS.Exp)
print(concord$rho.c$est)
#' **Concordance for dN-dS when only good window-sites considered = `r concord$rho.c$est`**
#' 

#' Plot actual versus expected
fig <- ggplot(gooddnds, aes(x=dNdS.Exp, y=dNdS.Act)) + 
  geom_point(alpha=0.5, shape=1) + 
  geom_smooth(method="lm", color="Red") + 
  xlab("\nExpected dN/dS") + 
  ylab("Inferred dN/dS\n") + 
  ggtitle("Scatterplot Expected vs Inferred dN/dS, Excl Low Syn")
print(fig)

fig <- ggplot(gooddnds, aes(x=dN_minus_dS.Exp, y=dN_minus_dS.Act)) + 
  geom_point(alpha=0.5, shape=1) + 
  geom_smooth(method="lm", color="Red") + 
  xlab("\nExpected dN/dS") + 
  ylab("Inferred dN/dS\n") + 
  ggtitle("Scatterplot Expected vs Inferred dN-dS, Exclude low Syn")
print(fig)

gooddnds$CovBin <- round(gooddnds$Cov.Act)
gooddnds$CovBin <- as.factor(gooddnds$CovBin)

fig <- ggplot(gooddnds, aes(x=EntropyCodon.Exp, y=LOD_dNdS, color=Cov.Act)) + 
  geom_point(alpha=0.5, shape=1) + 
  geom_smooth(method="lm", color="Red") + 
  xlab("\nSite Entropy") + 
  ylab("Log(Inferred dN/dS) - Log(Expected dN/dS)\n") + 
  ggtitle("Log Odds Difference dN/dS By Full Population Entropy and Coverage")
print(fig)

fig <- ggplot(gooddnds, aes(x=Subst.Exp, y=LOD_dNdS, color=Cov.Act)) + 
  geom_point(alpha=0.5, shape=1) + 
  geom_smooth(method="lm", color="Red") + 
  xlab("\nSite Entropy") + 
  ylab("Log(Inferred dN/dS) - Log(Expected dN/dS)\n") + 
  ggtitle("Log Odds Difference dN/dS By Full Population Site Substitutions and Coverage")
print(fig)


fig <- ggplot(gooddnds, aes(x=EntropyCodon.Exp, y=AbsLOD_dNdS, color=CovBin)) + 
  #geom_point(alpha=0.5, shape=1) + 
  geom_smooth(method="loess") + 
  xlab("\nSite Substitutions") + 
  ylab("|Log(Inferred dN/dS) - Log(Expected dN/dS)|\n") + 
  ggtitle("Log Odds Difference dN/dS By Full Population Entropy and Coverage")
print(fig)

fig <- ggplot(gooddnds, aes(x=Subst.Exp, y=AbsLOD_dNdS, color=CovBin)) + 
  #geom_point(alpha=0.5, shape=1) + 
  geom_smooth(method="loess") + 
  xlab("\nSite Entropy") + 
  ylab("|Log(Inferred dN/dS) - Log(Expected dN/dS)|\n") + 
  ggtitle("Log Odds Difference dN/dS By Full Population Site Substitutions and Coverage")
print(fig)




#' 
#' GLM for Predictors for Bad performance, Ignore Sites with N or S higher than zero but lower than 1.  
#' =========================================
#' 

# Expects continuous response for speedglm
fit_and_plot_glm_fast <- function(resp_colname, fit, df) {
  
  bestfit <- stepAIC(fit, direction="both", trace=FALSE)
  print(summary(bestfit))
  
  # speedglm doesn't expose residuals or fitted values. Do it ourselves
  # speedglm doesn't expose residuals or fitted values. Do it ourselves
  df_fit <- data.frame(Intercept=1, df[, attributes(bestfit$terms)$term.labels])
  
  # TODO:  hack - this is temp hack to get factors to be numeric
  if (length(grep("IsLowSubst.Act", colnames(df_fit)) > 0)){
    df_fit$IsLowSubst.Act <- ifelse (df_fit$IsLowSubst.Act == TRUE, 1, 0)  
  }
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
  print(coeffs)
  
  
  predictors <-  attributes(bestfit$terms)$term.labels
  
  figs <- sapply(predictors, 
                 function(varname) {
                   fig <- ggplot(df_fit, aes_string(x=varname, y="residuals")) + 
                     geom_point(alpha=0.5, shape=1) + 
                     geom_smooth(method="lm", color="Red") + 
                     xlab(varname) + 
                     ylab("Residuals\n") + 
                     ggtitle(paste0("Residuals of  ", resp_colname, " Vs ", varname))
                   print(fig)
                 })
  
  # Plot predicted values
  figs <- sapply(predictors, 
                 function(varname) {
                   fig <- ggplot(df_fit, aes_string(x=varname, y="fitted.values")) + 
                     geom_point(alpha=0.5, shape=1) + 
                     geom_smooth(method="lm", color="Red") + 
                     xlab(varname) + 
                     ylab("Predicted Values\n") + 
                     ggtitle(paste0("Predicted Values for  ", resp_colname, " From ", varname))
                   print(fig)
                 })
  
  return (bestfit)
}


cleandnds <- na.omit(dnds)
summary(cleandnds)
dim(cleandnds)

LODFormula <- as.formula(paste0("AbsLOD_dNdS~", paste0(LM_COVAR_NAMES, collapse=" + ")))
print(LODFormula)
allfitLOD <- speedglm(LODFormula, data=cleandnds)
print(summary(allfitLOD))
bestfit <- fit_and_plot_glm_fast("AbsLOD_dNdS", allfitLOD, df=cleandnds)



DistFormula <- as.formula(paste0("AbsDist_dn_minus_dS~", paste0(LM_COVAR_NAMES, collapse=" + ")))
print(DistFormula)
allfitDist <- speedglm(DistFormula, data=cleandnds)
print(summary(allfitDist))
bestfit <- fit_and_plot_glm_fast("AbsDist_dn_minus_dS", allfitDist, df=cleandnds)

