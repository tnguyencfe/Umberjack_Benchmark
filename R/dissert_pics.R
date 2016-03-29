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
#library(GGally)
#library(MASS)
#library(caret)
#library(gplots)
library(RColorBrewer)
#library(nortest)
#library(devtools)
#library(ggbiplot)
#library(directlabels)
#library(gamlss)
#library(logistf)
#library(speedglm)  # faster glms in parallel
#library(caret)  # for finding correlations
#library(gtools)  # for inv.logit
#library(arm)  # for binned.plot for logistic regresison residuals
#source ('./speedStepAIC.R')
source('./load_all_sim_dnds.R')


# Per Window-Site data
dnds_filename <- "/home/thuy/gitrepo/Umberjack_Benchmark/simulations/out/collate_all.recombo.csv"
if (exists("dnds_filename")) {
  recombo_dnds <- get_all_sim_dnds(dnds_filename)  
} else {
  recombo_dnds <- get_all_sim_dnds()
}
dim(dnds)
summary(dnds)

# Per Window data
recombo_window <-  get_window_sim_dnds(dnds=recombo_dnds)
dim(recombo_window)
summary(recombo_window)


# This removes outliers so that we can visualize the majority of points 
outlier_range <- function(x) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = TRUE, names=FALSE)
  fudge <- 1.5 * IQR(x, na.rm = TRUE)
  return (c(lower=qnt[1] - fudge, upper=qnt[2] + fudge))
}


#' Investigate Tree Accuracy vs Umberjack dnds Accuracy
#' ==============
#' 
#' 
recombo_window$LogTreeLen.Act <- log10(recombo_window$TreeLen.Act)
recombo_window$LogReads.Act <- log10(1+recombo_window$Reads.Act)
recombo_window$LogPolytomy.Act <- log10(1+recombo_window$Polytomy.Act)
recombo_window$LogPolytomyPerTreeLen.Act <- log10(1+recombo_window$Polytomy.Act/recombo_window$TreeLen.Act)
recombo_window$PolytomyPerRead.Act <- recombo_window$Polytomy.Act/recombo_window$Reads.Act
recombo_window$LogTreeLenPerRead.Act <- log10(recombo_window$TreeLen.Act/recombo_window$Reads.Act)
recombo_window$LogWindow_Subst.Act <- log10(1+recombo_window$Window_Subst.Act)
#recombo_window$TreeLenPerRead.Act <- recombo_window$TreeLen.Act/recombo_window$Reads.Act

# Plot response vs single variable, color points using gradient
plot_gradient <- function(data, resp_colname, var_colname, color_colnames=NULL) {  
  figs <- sapply(color_colnames, 
                 function(color_colname) {                   
                   fig <- ggplot(data, aes_string(x=var_colname, y=resp_colname)) + 
                     xlab(nice(var_colname)) + 
                     ylab(nice(resp_colname)) +                     
                     theme_bw() + 
                     geom_point(aes_string(color=color_colname), alpha=0.7, na.rm=TRUE) +                      
                     geom_smooth(method="lm", se=FALSE, color="black", size=2) +                      
                     scale_colour_gradient2(low="blue", mid="grey", high="red", midpoint=mean(data[, color_colname], na.rm=TRUE)) +
                     ggtitle("Inaccuracy Vs Covariate")
                   print(fig)
                 })
}

#+ fig.width=10
plot_gradient(data=recombo_window, resp_colname="WinSqDist_dn_minus_dS", 
              var_colname="TreeDistPerRead.Act",
              #color_colnames=c("LogTreeLen.Act", "LogPolytomy.Act", "LogPolytomyPerTreeLen.Act", "PolytomyPerRead.Act", "LogTreeLenPerRead.Act",
                               #"LogReads.Act", "LogWindow_Subst.Act", "WinP_SameCodonFreq.Act")
              color_colnames=c("LogTreeLen.Act"))

THESIS_DIR <- '../../MutationPatterns/tex'

#+ fig.width=10, fig.height=7
fig <- ggplot(recombo_window, aes(x=TreeDistPerRead.Act, y=WinSqDist_dn_minus_dS)) + 
  xlab("Normalized Window Ave WRF") + 
  ylab(expression(paste("Window Average ", Delta))) + 
  theme_bw(base_size=12) + 
  geom_point(aes(color=recombo_window$TreeLen.Act)) +                      
  geom_smooth(method="lm", se=FALSE, color="black", size=2) +                      
  scale_colour_gradient2(name="Tree\nLength", low="blue", mid="darkgrey", high="red", 
                         midpoint=quantile(recombo_window$TreeLen.Act, probs=0.5, na.rm=TRUE)) + 
  theme(axis.title=element_text(size=24), axis.text=element_text(size=20),
        legend.title=element_text(size=20))
print(fig)

ggsave(filename=paste0(THESIS_DIR, "/umberjack/error_v_recombo.png"), plot=fig)



###############  ALL SIMULATED DATA FOR RANDOM FOREST

dnds_filename <- "/home/thuy/gitrepo/Umberjack_Benchmark/simulations/out/collate_all.treedist.csv"
if (exists("dnds_filename")) {
  dnds <- get_all_sim_dnds(dnds_filename)  
} else {
  dnds <- get_all_sim_dnds()
}
dim(dnds)
summary(dnds)

# Per Window data
window <-  get_window_sim_dnds(dnds=dnds)
dim(window)
summary(window)

#' Plot the univariate regression plots, with outliers removed
#+ fig.width=10, fig.height=7
for (predictor in LM_COVAR_NAMES) {
  
  fig <- ggplot(dnds, aes_string(x=predictor, y="SqDist_dn_minus_dS")) + 
    xlab(nice(predictor)) + 
    ylab(expression(Delta)) + 
    theme_bw(base_size=12) + 
    geom_point(alpha=0.5, alpha=0.5) +                               
    geom_smooth(method="lm", se=FALSE, color="blue", size=2) +                      
    scale_y_continuous(limits=outlier_range(dnds$SqDist_dn_minus_dS)) + 
    #scale_colour_gradient2(name="Tree\nLength", low="blue", mid="darkgrey", high="red", 
    #                       midpoint=quantile(recombo_window$TreeLen.Act, probs=0.5, na.rm=TRUE)) + 
    theme(axis.title=element_text(size=20), axis.text=element_text(size=20),
          legend.title=element_text(size=20))
  print(fig)
  
  #ggsave(filename=paste0(THESIS_DIR, "/umberjack/error_v_recombo.png"), plot=fig)

}


#' Plot the univariate regression plots, with outliers intact
#+ fig.width=10, fig.height=7
for (predictor in LM_COVAR_NAMES) {
  
  fig <- ggplot(dnds, aes_string(x=predictor, y="SqDist_dn_minus_dS")) + 
    xlab(nice(predictor)) + 
    ylab(expression(Delta)) + 
    theme_bw(base_size=12) + 
    geom_point(alpha=0.5, alpha=0.5) +                      
    geom_smooth(method="lm", se=FALSE, color="blue", size=2) +                      
    #scale_y_continuous(limits=outlier(dnds$SqDist_dn_minus_dS)) + 
    #scale_colour_gradient2(name="Tree\nLength", low="blue", mid="darkgrey", high="red", 
    #                       midpoint=quantile(recombo_window$TreeLen.Act, probs=0.5, na.rm=TRUE)) + 
    theme(axis.title=element_text(size=20), axis.text=element_text(size=20),
          legend.title=element_text(size=20))
  print(fig)
  
  #ggsave(filename=paste0(THESIS_DIR, "/umberjack/error_v_recombo.png"), plot=fig)
  
}

