
# For multiple simulations, find the variables that affect accuracy

library(knitr)
opts_chunk$set(warning=FALSE, message=FALSE, width=1200, echo=TRUE)
options(width=150)

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

# Per Window-Site data
dnds <- get_all_sim_dnds()
dim(dnds)
summary(dnds)


# Per Window data
window <-  get_window_sim_dnds(dnds=dnds)
dim(window)
summary(window)





#' Error Rate:
#' =============================================
#' 
#' **Nucleotide Error Rate After Umberjack Quality Masking = `r mean(dnds$ErrPerCodon.Act)/3 `**
#' 





# #' Correlation Between Features
# #' =============================================
# #' 
# 
# # Find correlation between features.  Don't bother scaling.  The correlation heatmap is the same whether we scale & centre or not.
# corMat <- cor(dnds[, c(NUM_RESP_NAMES, COVAR_NAMES)], use="complete.obs", method="spearman")
# #+ fig.width=15, fig.height=15
# heatmap.2(corMat,  col=bluered, density.info="none", trace="none", srtCol=45, main="Feature Correlation", margins=c(12,12), 
#           colsep=(1:ncol(corMat)), rowsep=(1:nrow(corMat)), sepwidth=c(0.05, 0.05), sepcolor="black")
# 

#' Plot Inaccuracy Vs Numerical Variables That Apply at Window Level
#' ==================================================
#' 

# This removes outliers so that we can visualize the majority of points 
outlier_range <- function(x) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = TRUE, names=FALSE)
  fudge <- 1.5 * IQR(x, na.rm = TRUE)
  return (c(lower=qnt[1] - fudge, upper=qnt[2] + fudge))
}

plot_resp_vs_var <- function(data, resp_colname, var_colnames, color_colname=NULL) {
  resp_range <- outlier_range(data[, resp_colname])
  filter_data <- data[!is.na(data[, resp_colname]) & 
                        data[, resp_colname] >= resp_range["lower"] &
                        data[, resp_colname] <= resp_range["upper"], ]
  figs <- sapply(var_colnames, 
                 function(var_colname) {
                   if (!is.null(color_colname)) {
                     fig <- ggplot(filter_data, aes_string(x=var_colname, y=resp_colname, color=color_colname)) + 
                       guides(color=FALSE)
                   } else {
                     fig <- ggplot(filter_data, aes_string(x=var_colname, y=resp_colname))             
                   }
                   fig <- fig +             
                     xlab(nice(var_colname)) + 
                     ylab(nice(resp_colname)) + 
                     geom_point(shape=1, alpha=0.5, na.rm=TRUE) + 
                     geom_smooth(method="lm") + 
                     ggtitle("Inaccuracy Vs Covariate")
                   print(fig)
                 })
}
plot_resp_vs_var(data=window, resp_colname="WinSqDist_dn_minus_dS", var_colnames=WINDOW_COVAR_NAMES, color_colname="File")
#plot_resp_vs_var(data=window, resp_colname="WinAbsLOD_dNdS", var_colnames=WINDOW_COVAR_NAMES)

#' Plot robinson foulds with other confounding variables
plot_resp_vs_var(data=window, resp_colname="TreeDistPerRead.Act", var_colnames=WINDOW_COVAR_NAMES, color_colname="File")


#' Plot Inaccuracy Vs Numerical Variables That Apply at Window-Site Level
#' ==================================================

plot_resp_vs_var(data=dnds, resp_colname="SqDist_dn_minus_dS", var_colnames=WINDOW_SITE_COVAR_NAMES)
#plot_resp_vs_var(data=dnds, resp_colname="AbsLOD_dNdS", var_colnames=WINDOW_SITE_COVAR_NAMES)

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

#' **Now Remove Window-Sites with Phylogeny Substitutions Only Arising from Ambiguous Codons**  
#' 
gooddnds <- dnds[dnds$IsLowSubst.Act == FALSE, ]
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
fig <- ggplot(gooddnds[!is.na(gooddnds$dNdS.Exp) & !is.na(gooddnds$dNdS.Act), ], aes(x=dNdS.Exp, y=dNdS.Act)) + 
  geom_point(alpha=0.5, shape=1) + 
  geom_smooth(method="lm", color="Red") + 
  xlab("\nExpected dN/dS") + 
  ylab("Inferred dN/dS\n") + 
  ggtitle("Scatterplot Expected vs Inferred dN/dS, Excl Sites with Only Ambig Subs")
print(fig)

fig <- ggplot(gooddnds[!is.na(gooddnds$dN_minus_dS.Exp) & !is.na(gooddnds$dN_minus_dS.Act), ], aes(x=dN_minus_dS.Exp, y=dN_minus_dS.Act)) + 
  geom_point(alpha=0.5, shape=1) + 
  geom_smooth(method="lm", color="Red") + 
  xlab("\nExpected dN-dS") + 
  ylab("Inferred dN-dS\n") + 
  ggtitle("Scatterplot Expected vs Inferred dN-dS,Excl Sites with Only Ambig Subs")
print(fig)




#' 
#' GLM for Predictors for Umberjack Accuracy
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
    ggtitle("Predicted vs Residuals")
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
#   
#   figs <- sapply(predictors, 
#                  function(varname) {
#                    fig <- ggplot(df_fit, aes_string(x=varname, y="residuals")) + 
#                      geom_point(alpha=0.5, shape=1) + 
#                      geom_smooth(method="lm", color="Red") + 
#                      xlab(nice(varname)) + 
#                      ylab("Residuals\n") + 
#                      ggtitle(paste0("Residuals vs covariate"))
#                    print(fig)
#                  })
  
  return (bestfit)
}


cleandnds <- na.omit(dnds)
summary(cleandnds)
dim(cleandnds)

# LODFormula <- as.formula(paste0("AbsLOD_dNdS~", paste0(LM_COVAR_NAMES, collapse=" + ")))
# print(LODFormula)
# allfitLOD <- speedglm(LODFormula, data=cleandnds)
# print(summary(allfitLOD))
# bestfit <- fit_and_plot_glm_fast("AbsLOD_dNdS", allfitLOD, df=cleandnds)



DistFormula <- as.formula(paste0("AbsDist_dn_minus_dS~", paste0(LM_COVAR_NAMES, collapse=" + ")))
print(DistFormula)
allfitDist <- speedglm(DistFormula, data=cleandnds)
print(summary(allfitDist))
bestfit <- fit_and_plot_glm_fast("AbsDist_dn_minus_dS", allfitDist, df=cleandnds)


