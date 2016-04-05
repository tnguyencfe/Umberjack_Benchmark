
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





#' Error Rate:
#' =============================================
#' 
#' **Nucleotide Error Rate After Umberjack Quality Masking = `r mean(dnds$ErrPerCodon.Act)/3 `**
#' 

#' Plot Inaccuracy Vs Numerical Variables That Apply at Window Level
#' ==================================================
#' 


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
                       scale_y_continuous(limits=c(max(resp_range["lower"], min(data[, resp_colname], na.rm=TRUE)), 
                                                   min(resp_range["upper"], max(data[, resp_colname], na.rm=TRUE)))) + 
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
#+ fig.width=10
plot_resp_vs_var(data=window, resp_colname="WinSqDist_dn_minus_dS", var_colnames=WINDOW_COVAR_NAMES, color_colname="File")


#' Plot Inaccuracy Vs Numerical Variables That Apply at Window-Site Level
#' ==================================================
#+ fig.width=10
plot_resp_vs_var(data=dnds, resp_colname="SqDist_dn_minus_dS", var_colnames=WINDOW_SITE_COVAR_NAMES, color_colname="File")
#plot_resp_vs_var(data=dnds, resp_colname="AbsLOD_dNdS", var_colnames=WINDOW_SITE_COVAR_NAMES)




#' Concordance
#' ===========================================
#' 

# #' **Find Concordance with all original data**
# #' 
# concord <- epi.ccc(dnds[!is.na(dnds$dNdS.Act) & !is.na(dnds$dNdS.Exp), ]$dNdS.Act, 
#                    dnds[!is.na(dnds$dNdS.Act) & !is.na(dnds$dNdS.Exp), ]$dNdS.Exp)
# print(concord$rho.c)
# print(concord$rho.c$est)
# #' **Concordance for dN/dS when all considered = `r concord$rho.c$est`**
# #' 

concord <- epi.ccc(dnds[!is.na(dnds$dN_minus_dS.Act) & !is.na(dnds$dN_minus_dS.Exp), ]$dN_minus_dS.Act, 
                   dnds[!is.na(dnds$dN_minus_dS.Act) & !is.na(dnds$dN_minus_dS.Exp), ]$dN_minus_dS.Exp)
print(concord$rho.c)
print(concord$rho.c$est)
#' **Concordance for dN-dS when all considered = `r concord$rho.c$est`**
#' 

# #' **Now Remove Window-Sites with Phylogeny Substitutions Only Arising from Ambiguous Codons**  
# #' 
# gooddnds <- dnds[dnds$IsLowSubst.Act == FALSE, ]
# summary(gooddnds)
# dim(gooddnds)
# 
# concord <- epi.ccc(gooddnds[!is.na(gooddnds$dNdS.Act) & !is.na(gooddnds$dNdS.Exp), ]$dNdS.Act, 
#                    gooddnds[!is.na(gooddnds$dNdS.Act) & !is.na(gooddnds$dNdS.Exp), ]$dNdS.Exp)
# print(concord$rho.c)
# print(concord$rho.c$est)
# #' **Concordance for dN/dS when only good window-sites considered = `r concord$rho.c$est`**
# #' 
# 
# concord <- epi.ccc(gooddnds[!is.na(gooddnds$dN_minus_dS.Act) & !is.na(gooddnds$dN_minus_dS.Exp), ]$dN_minus_dS.Act, 
#                    gooddnds[!is.na(gooddnds$dN_minus_dS.Act) & !is.na(gooddnds$dN_minus_dS.Exp), ]$dN_minus_dS.Exp)
# print(concord$rho.c)
# print(concord$rho.c$est)
# #' **Concordance for dN-dS when only good window-sites considered = `r concord$rho.c$est`**
# #' 


#' **concordance by Dataset**

concord <- ddply(.data=dnds, .variables="File", .fun=function(x) {
  data.frame(concord=epi.ccc(x$dN_minus_dS.Act, x$dN_minus_dS.Exp)$rho.c$est)
})

#+ results="asis"
kable(concord, format="html", caption="Concordance by dataset")

# #'
# #' **concordance by file, ignore sites with only ambiguous substitutions**
# #' 
# concord <- ddply(.data=gooddnds, .variables="File", .fun=function(x) {
#   data.frame(concord=epi.ccc(x$dN_minus_dS.Act, x$dN_minus_dS.Exp)$rho.c$est)
# })
# 
# #+ results="asis"
# kable(concord, format="html", caption="Concordance by dataset, Excluding sites with Only  Ambiguous Subs")

#'
#' Plot expected values from each file
#+ fig.width=10
fig <- ggplot(dnds, aes(y=dN_minus_dS.Exp, x=File, color=File)) + 
  geom_boxplot() + 
  guides(color=FALSE) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  ylab("expected dn-ds") + 
  xlab("Dataset") + 
  ggtitle("Expected dnds from each dataset")
print (fig)

#' Plot diff from expected values from each file
#+ fig.width=12
fig <- ggplot(dnds, aes(x=(1+SqDist_dn_minus_dS), color=File)) + 
  geom_density() + 
  scale_x_log10() + 
  xlab("Sq diff from expected dn-ds") + 
  ggtitle("Sq diff from Expected dnds from each dataset")
print (fig)

fig <- ggplot(dnds, aes(y=(1+SqDist_dn_minus_dS), color=File, x=File)) + 
  geom_boxplot()  + 
  guides(color=FALSE) + 
  scale_y_log10() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  ylab("Sq diff from expected dn-ds") + 
  xlab("Dataset") + 
  ggtitle("sq diff from Expected dnds from each dataset")
print (fig)

#' Plot diff from expected values from each file
#+ fig.width=15, fig.height=12
fig <- ggplot(dnds, aes(x=CodonSite, y=(1+SqDist_dn_minus_dS), color=File)) +   
  geom_line() + 
  xlab("CodonSite") + 
  ylab("Sq Diff from Expected dn-ds")+ 
  scale_y_log10() + 
  ggtitle("Sq diff from Expected dnds from each dataset")

sm_fig <- ggplot(dnds, aes(x=CodonSite, y=(1+SqDist_dn_minus_dS), color=File)) +   
  geom_smooth(se=FALSE) + 
  xlab("CodonSite") + 
  ylab("Sq Diff from Expected dn-ds")+ 
  scale_y_log10() + 
  ggtitle("Sq diff from Expected dnds from each dataset, smoothed")

exp_fig <- ggplot(dnds, aes(x=CodonSite, y=dN_minus_dS.Exp, color=File)) +   
  geom_smooth(se=FALSE) + 
  xlab("CodonSite") + 
  ylab("Expected dn-ds")+ 
  ggtitle("Expected dn-ds from each dataset, smoothed")


list_gps <- AlignPlots(fig, sm_fig, exp_fig)
do.call(grid.arrange, args=c(list_gps, ncol=1))

#' Plot actual versus expected
# fig <- ggplot(gooddnds[!is.na(gooddnds$dNdS.Exp) & !is.na(gooddnds$dNdS.Act), ], aes(x=dNdS.Exp, y=dNdS.Act)) + 
#   geom_point(alpha=0.5, shape=1) + 
#   geom_smooth(method="lm", color="Red") + 
#   xlab("\nExpected dN/dS") + 
#   ylab("Inferred dN/dS\n") + 
#   ggtitle("Scatterplot Expected vs Inferred dN/dS, Excl Sites with Only Ambig Subs")
# print(fig)

#+ fig.width=12, fig.height=12
# fig <- ggplot(gooddnds[!is.na(gooddnds$dN_minus_dS.Exp) & !is.na(gooddnds$dN_minus_dS.Act), ], aes(x=dN_minus_dS.Exp, y=dN_minus_dS.Act, color=File)) + 
#   geom_point(alpha=0.5, shape=1) + 
#   #geom_smooth(method="lm", color="Red") + 
#   geom_smooth(method="lm") + 
#   xlab("\nExpected dN-dS") + 
#   ylab("Inferred dN-dS\n") + 
#   ggtitle("Scatterplot Expected vs Inferred dN-dS,Excl Sites with Only Ambig Subs")
# print(fig)


#+ fig.width=12, fig.height=12
fig <- ggplot(dnds[!is.na(dnds$dN_minus_dS.Exp) & !is.na(dnds$dN_minus_dS.Act), ], aes(x=dN_minus_dS.Exp, y=dN_minus_dS.Act, color=File)) + 
  geom_point(alpha=0.5, shape=1) + 
  #geom_smooth(method="lm", color="Red") + 
  geom_smooth(method="lm", se=FALSE) + 
  xlab("\nExpected dN-dS") + 
  ylab("Inferred dN-dS\n") + 
  ggtitle("Scatterplot Expected vs Inferred dN-dS")
print(fig)
