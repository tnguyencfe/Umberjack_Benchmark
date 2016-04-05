
# For multiple simulations, find the variables that affect accuracy

library(knitr)
opts_chunk$set(warning=FALSE, message=FALSE, width=1200, echo=TRUE)
options(width=150)

#+ message=FALSE
# From data from all windows, aggregates by averaging over windows
library(ggplot2)
library(reshape2)
library(epiR) # concordance
library(plyr)
#library(stats)
#library(psych)
#library(GGally)
library(MASS)  # stepAIC
#library(gplots)
library(RColorBrewer)
library(energy)  # dcor()
#library(nortest)
#library(devtools)
#library(ggbiplot)
#library(directlabels)
#library(gamlss)
#library(logistf)
library(speedglm)  # faster glms in parallel
#library(caret)  # for finding correlations
#library(gtools)  # for inv.logit
#library(arm)  # for binned.plot for logistic regresison residuals
library(miscTools) # for rsquared
source ('./speedStepAIC.R')
source('./load_all_sim_dnds.R')
source("plot_helper.R")  # sliding window git repo

# Per Window-Site data
if (exists("dnds_filename")) {
  dnds <- get_all_sim_dnds(dnds_filename)  
} else {
  dnds <- get_all_sim_dnds()
}
# dim(dnds)
summary(dnds)

# Per Window data
window <-  get_window_sim_dnds(dnds=dnds)
dim(window)
summary(window)

cleandnds <- na.omit(dnds)
summary(cleandnds)
dim(cleandnds)

# hack for backwards compatibility when we didn't auto remove resolved codons from Umberjack
# Or there are no resolved substitutions, then don't bother using it in the GLM
if (!"ResolvedPerSub.Act" %in% colnames(dnds) |  sum(cleandnds$ResolvedPerSub.Act > 0, na.rm=TRUE) == 0) {
  LM_COVAR_NAMES <-  LM_COVAR_NAMES[!LM_COVAR_NAMES %in% c("ResolvedPerSub.Act" )]  
}



#' Univariate Relations
#' ========================

#' Individual predictors
#' 
for (predictor in LM_COVAR_NAMES) {
  
  oneform <- as.formula(paste0("GammaSqDist_dn_minus_dS ~ ", predictor))
  print(oneform)
  
  onefit <- glm(oneform, data=cleandnds, family=Gamma(), 
                start=rep(1, 2))
  print(summary(onefit))
  plot(onefit)
}

#' Individual predictors Gaussian.  These have large deviance than Gamma. Don't use.
#' 
for (predictor in LM_COVAR_NAMES) {
  
  oneform <- as.formula(paste0("SqDist_dn_minus_dS ~ ", predictor))
  print(oneform)
  
  onefit <- glm(oneform, data=cleandnds, family=gaussian, 
                start=rep(1, 2))
  print(summary(onefit))
  plot(onefit)
}

#' Find the spearman's correlation, distance correlation between individual predictor and the response
#' 
univar_rank <- data.frame(predictor=as.character(LM_COVAR_NAMES))
univar_rank$spear_cor <- sapply(as.character(univar_rank$predictor),
                                function(predictor) {
                                  print(predictor)
                                  cor(cleandnds$SqDist_dn_minus_dS, cleandnds[, predictor], method="spearman", use="complete.obs")
                                })
univar_rank$abs_spear_cor <- abs(univar_rank$spear_cor)
univar_rank <- univar_rank[order(-univar_rank$abs_spear_cor), ]

#+ results="asis"
kable(univar_rank, format="html", caption="spearman corr univarate ranking")

# Don't calculate distance correlation via dcor(), will crash due to excessive memory usage if >5000 datapoints
# for (predictor in LM_COVAR_NAMES) {
#   pred_dcor <- dcor(cleandnds$SqDist_dn_minus_dS, cleandnds[, predictor])
#   print(paste0(predictor, " dcor =", pred_dcor))
# }


#' 
#' GLM for Predictors for Umberjack Accuracy
#' =========================================
#' 

# Expects speedglm fit as input.
# Does stepwise AIC backwards selection.
fit_and_plot_glm_fast <- function(resp_colname, fit, df) {
  
  bestfit <- stepAIC(fit, direction="backward", trace=1)
  
  print(summary(bestfit))
  
  # speedglm doesn't expose residuals or fitted values. Do it ourselves
  # speedglm doesn't expose residuals or fitted values. Do it ourselves
  df_fit <- data.frame(Intercept=1, df[, attributes(bestfit$terms)$term.labels])
  
  # TODO:  hack - this is temp hack to get factors to be numeric
  if (length(grep("IsLowSubst.Act", colnames(df_fit)) > 0)){
    df_fit$IsLowSubst.Act <- ifelse (df_fit$IsLowSubst.Act == TRUE, 1, 0)  
  }
  
  if (!"residuals" %in% names(bestfit)) {
    df_fit$fitted.values <- as.vector(as.matrix(df_fit) %*% coef(bestfit))
    df_fit[, resp_colname] <-  df[, resp_colname]    
    df_fit$residuals <-  df_fit[, resp_colname] - df_fit$fitted.values  
  } else {
    df_fit$fitted.values <- bestfit$fitted.values
    df_fit[, resp_colname] <-  df[, resp_colname]    
    df_fit$residuals <-  bestfit$residuals  
  }
  
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
  
  # Print R-square values
  r2 <- rSquared(y=df_fit$fitted.values, resid=df_fit$residuals)
  print (paste0("rsqured = ", r2))
  
  coeffs <- data.frame(summary(bestfit)$coefficients)
  colnames(coeffs)[grep("Pr", colnames(coeffs))] <- "pval"
  coeffs$pval <- as.numeric(as.character(coeffs$pval))
  coeffs$adj.pval <- p.adjust(coeffs$pval, method="BH")
  coeffs$name <- rownames(coeffs)
  print(coeffs)
  
  
  predictors <-  attributes(bestfit$terms)$term.labels
  return (bestfit)
}

# Expects continuous response for Gamma GLM
# This was hacked up because Gamma GLM fits would crash unless we input start values.
# But we could never actually get the full backwards stepwise AIC to work because it couldn't fit certain nested models anyway.
fit_and_plot_glm_fast_gamma <- function(resp_colname, fit, df) {
  
  bestfit <- mystepAIC(fit, direction="backward", trace=TRUE, use.start=TRUE)
  
  print(summary(bestfit))
  
  # speedglm doesn't expose residuals or fitted values. Do it ourselves
  df_fit <- data.frame(Intercept=1, df[, attributes(bestfit$terms)$term.labels])
  
  # TODO:  hack - this is temp hack to get factors to be numeric
  if (length(grep("IsLowSubst.Act", colnames(df_fit)) > 0)){
    df_fit$IsLowSubst.Act <- ifelse (df_fit$IsLowSubst.Act == TRUE, 1, 0)  
  }
  if (!"residuals" %in% names(bestfit)) {
    df_fit$fitted.values <- as.vector(as.matrix(df_fit) %*% coef(bestfit))
    df_fit[, resp_colname] <-  df[, resp_colname]    
    df_fit$residuals <-  df_fit[, resp_colname] -df_fit$fitted.values  
  } else {
    df_fit$fitted.values <- bestfit$fitted.values
    df_fit[, resp_colname] <-  df[, resp_colname]    
    df_fit$residuals <-  bestfit$residuals  
  }
  
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
  
  # Print R-square values
  r2 <- rSquared(y=df_fit$fitted.values, resid=df_fit$residuals)  # from miscTools package
  print (paste0("rsqured = ", r2))
  
  
  coeffs <- data.frame(summary(bestfit)$coefficients)
  colnames(coeffs)[grep("Pr", colnames(coeffs))] <- "pval"
  coeffs$pval <- as.numeric(as.character(coeffs$pval))
  coeffs$adj.pval <- p.adjust(coeffs$pval, method="BH")
  coeffs$name <- rownames(coeffs)
  print(coeffs)
  
  
  predictors <-  attributes(bestfit$terms)$term.labels
  return (bestfit)
}





DistFormula <- as.formula(paste0("SqDist_dn_minus_dS~", paste0(LM_COVAR_NAMES, collapse=" + ")))
print(DistFormula)

#' Gaussian Fit, Backwards Feature Selection
#' 
allfitDist <- speedglm(DistFormula, data=cleandnds, family=gaussian(), fitted=TRUE)
bestfit <- fit_and_plot_glm_fast("SqDist_dn_minus_dS", allfitDist, df=cleandnds)



# #' Gamma fit, speedglm:  this will fail
# allfitDist <- speedglm(DistFormula, data=cleandnds, family=Gamma()
#                        #start=rep(1, length(LM_COVAR_NAMES) + 1)
#                        )
# 
# #' GLM gamma fit will succeed, but backwards stepwise AIC will crash because some GLM will be unable to fit model
# allfitDist <- glm(DistFormula, data=cleandnds,  family=Gamma(),
#                        start=rep(1, length(LM_COVAR_NAMES) + 1))
# 
# bestfit <- fit_and_plot_glm_fast_gamma("SqDist_dn_minus_dS", allfitDist, df=cleandnds)


#' In order to use Gamma family, the response must be in (0, inf).  That is, response can not be zero.  
#  We decided not to use tweedie or zero-inflated gamma response because there aren't enough datapoints where response = 0 to make a difference.
dim(cleandnds[cleandnds$SqDist_dn_minus_dS == 0, ])
dim(cleandnds)

#' Fraction of datapoints in which response = 0 = `r dim(cleandnds[cleandnds$SqDist_dn_minus_dS == 0, ])` /  `r dim(cleandnds) ` = 
#' `r dim(cleandnds[cleandnds$SqDist_dn_minus_dS == 0, ]) /   dim(cleandnds) `
#' 
cleandnds$GammaSqDist_dn_minus_dS <- cleandnds$SqDist_dn_minus_dS
cleandnds$GammaSqDist_dn_minus_dS[cleandnds$GammaSqDist_dn_minus_dS == 0] <- 1e-13
summary(cleandnds)

#' Gamma, normal GLM.  This will fail too.  
#' 
gammaDistFormula <- as.formula(paste0("GammaSqDist_dn_minus_dS ~", paste0(LM_COVAR_NAMES, collapse=" + ")))
print(gammaDistFormula)

allfitDist <- glm(gammaDistFormula, data=cleandnds, family=Gamma(),
                       start=rep(1, length(LM_COVAR_NAMES) + 1))

print(summary(allfitDist))
plot(allfitDist)

bestfit <- fit_and_plot_glm_fast("GammaSqDist_dn_minus_dS", allfitDist, df=cleandnds)



#library(biglm)
#bigglm(formula, data, family=gaussian(),...)
#bigfit <- bigglm(DistFormula, data=cleandnds, family=Gamma())
#bigfit <- bigglm(SqDist_dn_minus_dS ~ EntropyCodon.Act , data=cleandnds, family=Gamma())
# Gives:
# Error in coef.bigqr(object$qr) : 
#   NA/NaN/Inf in foreign function call (arg 3)



