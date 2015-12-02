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


#' Tree Topology Accuracy
#' ==================================
#' 
#' How does tree topology accuracy compare to location of breakpoints?
#' Break Ratio = sum across breakpoints [bases on left of breakpoint / bases on right of breakoint]

# window info is repeated for every codon site in collate_dnds.  Isolate just the window values
window <- aggregate(cbind(BreakRatio.Act, TreeDistPerRead.Act) ~ File + Window_Start, 
                    data=dnds, 
                    FUN=function(x) {
                      result <- unique(x)
                      if (length(result) > 1) {
                        stop("There should only be 1 unique value per window")
                      }
                      return (result)
                    })


fig <- ggplot(window, aes(x=BreakRatio.Act, y=TreeDistPerRead.Act)) + 
  geom_point(alpha=0.7, shape=1) + 
  geom_smooth(method="lm") + 
  xlab("Total Break Ratio in Window") + 
  ylab("Robinson Folds Distance from True Tree\n Normalized by Total Tips") + 
  ggtitle("Location of Breakpoint vs Tree Accuracy")
print (fig)


