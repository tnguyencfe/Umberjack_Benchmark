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
#library(psych)
#library(GGally)
#library(MASS)
#library(caret)
library(gplots)
library(RColorBrewer)
#library(nortest)
library(devtools)
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

dnds_filename <- "/home/thuy/gitrepo/Umberjack_Benchmark/simulations/out/collate_all.recombo.csv"
# Per Window-Site data
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


# This removes outliers so that we can visualize the majority of points 
outlier_range <- function(x) {
  qnt <- quantile(x, probs=c(.25, .75), na.rm = TRUE, names=FALSE)
  fudge <- 1.5 * IQR(x, na.rm = TRUE)
  return (c(lower=qnt[1] - fudge, upper=qnt[2] + fudge))
}

plot_resp_vs_var <- function(data, resp_colname, var_colnames, color_colname=NULL) {
  resp_range <- outlier_range(data[, resp_colname])
  figs <- sapply(var_colnames, 
                 function(var_colname) {
                   if (!is.null(color_colname)) { 
                     colourCount <- length(levels(data[, color_colname]))
                     getPalette <- colorRampPalette(brewer.pal(9, "Set1"))
                     
                     fig <- ggplot(data, aes_string(x=var_colname, y=resp_colname)) + 
                       xlab(nice(var_colname)) + 
                       ylab(nice(resp_colname)) +
                       scale_colour_manual(values = getPalette(colourCount)) +
                       theme_bw() + 
                       geom_point(aes_string(color=color_colname), alpha=0.7, na.rm=TRUE) + 
                       geom_smooth(aes_string(color=color_colname, group=color_colname), method="lm", se=FALSE) + 
                       geom_smooth(method="lm", se=FALSE, color="black", size=2) + 
#                        scale_y_continuous(limits=c(max(resp_range["lower"], min(data[, resp_colname], na.rm=TRUE)), 
#                                                    min(resp_range["upper"], max(data[, resp_colname], na.rm=TRUE)))) + 
                       ggtitle("Inaccuracy Vs Covariate")
                     
                     if (colourCount > 12) {  # If there are way too many factors, there's no point in coloring them.
                       fig <- fig + guides(color=FALSE)
                     }
                     print(fig)
                   } else {                          
                     fig <- ggplot(filter_data, aes_string(x=var_colname, y=resp_colname)) + 
                       theme_bw() + 
                       xlab(nice(var_colname)) + 
                       ylab(nice(resp_colname)) + 
                       geom_point(shape=1, alpha=0.5, na.rm=TRUE) +                        
                       geom_smooth(method="lm", se=FALSE) + 
                       ggtitle("Inaccuracy Vs Covariate")
                     print(fig)
                   }
                   
                 })
}



#' Plot robinson foulds with other confounding variables
#' 
#+ fig.width=10
plot_resp_vs_var(data=window, resp_colname="TreeDistPerRead.Act", var_colnames=WINDOW_COVAR_NAMES, color_colname="File")


#' Investigate Tree Accuracy vs Umberjack dnds Accuracy
#' ==============
#' 
#' 
window$LogTreeLen.Act <- log10(window$TreeLen.Act)
window$LogReads.Act <- log10(1+window$Reads.Act)
window$LogPolytomy.Act <- log10(1+window$Polytomy.Act)
window$LogPolytomyPerTreeLen.Act <- log10(1+window$Polytomy.Act/window$TreeLen.Act)
window$PolytomyPerRead.Act <- window$Polytomy.Act/window$Reads.Act
window$LogTreeLenPerRead.Act <- log10(window$TreeLen.Act/window$Reads.Act)
window$LogWindow_Subst.Act <- log10(1+window$Window_Subst.Act)

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
plot_gradient(data=window, resp_colname="WinSqDist_dn_minus_dS", 
              var_colname="TreeDist.Act",
              color_colnames=c("LogTreeLen.Act", "LogPolytomy.Act", "LogPolytomyPerTreeLen.Act", "PolytomyPerRead.Act", "LogTreeLenPerRead.Act",
                               "LogReads.Act", "LogWindow_Subst.Act", "WinP_SameCodonFreq.Act", "Window_Breaks"))

#+ fig.width=10
plot_gradient(data=window, resp_colname="WinSqDist_dn_minus_dS", 
              var_colname="TreeDistPerRead.Act",
              color_colnames=c("LogTreeLen.Act", "LogPolytomy.Act", "LogPolytomyPerTreeLen.Act", "PolytomyPerRead.Act", "LogTreeLenPerRead.Act",
                               "LogReads.Act", "LogWindow_Subst.Act", "WinP_SameCodonFreq.Act", "Window_Breaks"))

#' Plot the window average accuracy colored by the Simulated Dataset
#+ fig.width=15                                   
fig <- ggplot(window, aes(x=TreeDistPerRead.Act, y=WinSqDist_dn_minus_dS)) + 
 xlab(nice("TreeDistPerRead.Act")) + 
 ylab(nice("WinSqDist_dn_minus_dS")) +                     
 theme_bw() + 
 geom_point(aes(color=File), alpha=0.7, na.rm=TRUE) +                      
 geom_smooth(method="lm", se=FALSE, color="black", size=2) +                      
 #scale_colour_gradient2(low="blue", mid="grey", high="red", midpoint=mean(window$File, na.rm=TRUE)) +
 ggtitle("Inaccuracy By Dataset")
print(fig)

