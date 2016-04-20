
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
dim(dnds)
head(dnds)
summary(dnds)

# Per Window data
window <-  get_window_sim_dnds(dnds=dnds)
dim(window)
summary(window)



ave_dnds <- ddply(.data=dnds,
              .variables=c("File", "CodonSite"),
              .fun=function(x) {
                if (length(unique(x$dN_minus_dS.Exp)) > 1) {
                  stop("there should only be one expected dn-ds per site")
                }
                data.frame(Ave.dN_minus_dS.Act = weighted.mean(x$dN_minus_dS.Act,
                                                        x$UnambigCodonRate.Act * x$Reads.Act, na.rm=TRUE),
                           dN_minus_dS.Exp = x$dN_minus_dS.Exp[1])
              })
summary(ave_dnds)
head(ave_dnds)
dim(ave_dnds)



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

nona_dnds <- subset(dnds, !is.na(dnds$dN_minus_dS.Exp) & !is.na(dnds$dN_minus_dS.Act))
rsq <- rSquared(y=nona_dnds$dN_minus_dS.Exp, resid = nona_dnds$dN_minus_dS.Act - nona_dnds$dN_minus_dS.Exp)
print(rsq)
#' ** Explained variance  = `r rsq `**
#

#' Averaged site dnds
ave_site_dnds <- ddply(.data=nona_dnds[nona_dnds$UnambigCodonRate.Act > 0, ],
                      .variables=c("File", "CodonSite"),
                      .fun=function(x) {
                        if (length(unique(x$dN_minus_dS.Exp)) > 1) {
                          stop("there should only be one expected dn-ds per site")
                        }
                        data.frame(Ave.dN_minus_dS.Act = weighted.mean(x$dN_minus_dS.Act,
                                                                       x$UnambigCodonRate.Act * x$Reads.Act, na.rm=TRUE),
                                   Windows = length(x$Window_Start),
                                   dN_minus_dS.Exp = x$dN_minus_dS.Exp[1])
                      })

summary(ave_site_dnds)
head(ave_site_dnds)
dim(ave_site_dnds)

rsq <- rSquared(y=ave_site_dnds$Ave.dN_minus_dS.Act, resid = (ave_site_dnds$Ave.dN_minus_dS.Act - ave_site_dnds$dN_minus_dS.Exp))
print(rsq)
concord <- epi.ccc(ave_site_dnds$Ave.dN_minus_dS.Act, ave_site_dnds$dN_minus_dS.Exp)
print(concord$rho.c)

aved_ideal <- ddply(ave_site_dnds, "File", function(x) {
  concord <- epi.ccc(x$Ave.dN_minus_dS.Act, x$dN_minus_dS.Exp)
  data.frame(concord_low=concord$rho$lower,
             concord_est=concord$rho$est,
             concord_hi=concord$rho$upper)
})
write.table(aved_ideal, file="site_ave_dnds_concord_by_dataset.csv", sep=",", quote=FALSE)

#' **Now Remove Window-Sites with Phylogeny Substitutions Only Arising from Ambiguous Codons, etc**  
#' 
gooddnds <- dnds[dnds$Subst.Act >= 5 & 
                   dnds$UnambigCodonRate.Act >= 0.8 &
                   dnds$Reads.Act >= 50 &
                   dnds$ES.Act >= 0.5 & 
                   #dnds$Window_Subst.Act > 10 &
                   !is.na(dnds$dN_minus_dS.Act) & !is.na(dnds$dN_minus_dS.Exp),]
summary(gooddnds)
dim(gooddnds)
head(gooddnds)

#' Missing window-site predictions = `r length(gooddnds$dN_minus_dS.Act) / sum( !is.na(dnds$dN_minus_dS.Act)) `
#' 
length(gooddnds$dN_minus_dS.Act)
sum( !is.na(dnds$dN_minus_dS.Act)) 
length(gooddnds$dN_minus_dS.Act)/sum( !is.na(dnds$dN_minus_dS.Act)) 

#' Missing site predictions = `r nrow(unique(gooddnds[, c("File", "CodonSite")]))/nrow(unique(dnds[!is.na(dnds$dN_minus_dS.Act), c("File", "CodonSite")]))`
nrow(unique(gooddnds[, c("File", "CodonSite")]))
nrow(unique(dnds[!is.na(dnds$dN_minus_dS.Act), c("File", "CodonSite")]))
nrow(unique(gooddnds[, c("File", "CodonSite")]))/nrow(unique(dnds[!is.na(dnds$dN_minus_dS.Act), c("File", "CodonSite")]))

rsq <- rSquared(y=gooddnds$dN_minus_dS.Exp, resid = gooddnds$dN_minus_dS.Act - gooddnds$dN_minus_dS.Exp)
print(rsq)

#' **Concordance for dn-ds when only good window-sites considered = `r concord$rho.c$est`**
#' 

concord <- epi.ccc(gooddnds$dN_minus_dS.Act, gooddnds$dN_minus_dS.Exp)
print(concord$rho.c)
#' **Concordance for dN-dS when only good window-sites considered = `r concord$rho.c$est`**
#' 

cor(gooddnds$dN_minus_dS.Act, gooddnds$dN_minus_dS.Exp, use="complete.obs", method="spearman")

#' Averaged site dnds
ave_gooddnds <- ddply(.data=gooddnds,
                      .variables=c("File", "CodonSite"),
                      .fun=function(x) {
                        if (length(unique(x$dN_minus_dS.Exp)) > 1) {
                          stop("there should only be one expected dn-ds per site")
                        }
                        data.frame(Ave.dN_minus_dS.Act = weighted.mean(x$dN_minus_dS.Act,
                                                                       x$UnambigCodonRate.Act * x$Reads.Act, na.rm=TRUE),
                                   dN_minus_dS.Exp = x$dN_minus_dS.Exp[1])
                      })
summary(ave_gooddnds)
head(ave_gooddnds)
dim(ave_gooddnds)

#' Filtered ave site dn-ds versus expected
#' 
fig <- ggplot(ave_gooddnds, aes(x=dN_minus_dS.Exp, y=Ave.dN_minus_dS.Act)) + 
  geom_abline(slope=1, color="red") + 
  geom_point(alpha=0.2) + 
  geom_smooth(method="lm") +   
  ylab("Inferred dN-dS") + 
  xlab("Expected dN-dS") + 
  #ggtitle("Filtered Ave Site dnds vs Exp") + 
  theme_bw() + 
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=16)) 
print(fig)

#' Rsquared
rsq <- rSquared(y=ave_gooddnds$dN_minus_dS.Exp, resid = ave_gooddnds$Ave.dN_minus_dS.Act - ave_gooddnds$dN_minus_dS.Exp)
print(rsq)

#' **Concordance for ave dn-ds when only good window-sites considered = `r concord$rho.c$est`**
#' 
concord <- epi.ccc(ave_gooddnds$Ave.dN_minus_dS.Act, ave_gooddnds$dN_minus_dS.Exp)
print(concord$rho.c)

cor(ave_gooddnds$Ave.dN_minus_dS.Act, ave_gooddnds$dN_minus_dS.Exp, use="complete.obs", method="spearman")



#' Choose best filter
#' ------------------- 
#' 
#' See affect of rsquared and lins concordance on different filters
#' Choose filter based on elbow of rsqred vs non-missing predictions.
#' 
        
subst <- c(1, 2, 5, 10)
unambig <-  c(0.7, 0.8, 0.9)
reads <- c(10, 50, 100)
es <- c(0.01, 0.5, 1)
filters <- expand.grid(subst=subst, unambig=unambig, reads=reads, es=es)        

do_filter_perf <- function(f) {  
  filterdnds <- dnds[dnds$Subst.Act >= f$subst & 
                       dnds$UnambigCodonRate.Act >= f$unambig &
                       dnds$Reads.Act >= f$reads &
                       dnds$ES.Act >= f$es & 
                       !is.na(dnds$dN_minus_dS.Act) & !is.na(dnds$dN_minus_dS.Exp),]
  
  filter_rsq <- rSquared(y=filterdnds$dN_minus_dS.Exp, resid = filterdnds$dN_minus_dS.Act - filterdnds$dN_minus_dS.Exp)
  filter_concord <- epi.ccc(filterdnds$dN_minus_dS.Act, filterdnds$dN_minus_dS.Exp)
  filter_has_predict <- length(filterdnds$dN_minus_dS.Act)/sum( !is.na(dnds$dN_minus_dS.Act)) 
    
  filterave_gooddnds <- ddply(.data=filterdnds,
                              .variables=c("File", "CodonSite"),
                              .fun=function(x) {
                                if (length(unique(x$dN_minus_dS.Exp)) > 1) {
                                  stop("there should only be one expected dn-ds per site")
                                }
                                data.frame(Ave.dN_minus_dS.Act = weighted.mean(x$dN_minus_dS.Act,
                                                                               x$UnambigCodonRate.Act * x$Reads.Act, na.rm=TRUE),
                                           dN_minus_dS.Exp = x$dN_minus_dS.Exp[1])
                              })
  
  filterave_rsq <- rSquared(y=filterave_gooddnds$dN_minus_dS.Exp, resid = filterave_gooddnds$Ave.dN_minus_dS.Act - filterave_gooddnds$dN_minus_dS.Exp)  
  filterave_concord <- epi.ccc(filterave_gooddnds$Ave.dN_minus_dS.Act, filterave_gooddnds$dN_minus_dS.Exp)
 
  orig_num_siteave_rows <- nrow(unique(dnds[!is.na(dnds$dN_minus_dS.Act), c("File", "CodonSite")]))
  filterave_has_predict <- length(filterave_gooddnds$Ave.dN_minus_dS.Act) / orig_num_siteave_rows
  
  return (data.frame(filter_rsq = filter_rsq,
                     filter_concord = filter_concord$rho.c$est,                     
                     filter_has_predict = filter_has_predict,
                     filterave_rsq = filterave_rsq,
                     filterave_concord = filterave_concord$rho.c$est,
                     filterave_has_predict = filterave_has_predict
    ))
}

filename <- "./filter_perf.csv"
if (file.exists(filename)) {
  warning(paste0("Not regenerating filter performance ", filename))
  filter_perf <- read.table(filename, sep=",", header=TRUE)
} else {
  filter_perf <-  adply(.data=filters,
                        .margins=1,
                        .fun=do_filter_perf)
  
  write.table(filter_perf, filename, sep=",", row.names=FALSE, quote=FALSE)
}
filter_perf <- filter_perf[filter_perf$subst > 1,] # this was retarded.  1 sub is obviously not enough.
head(filter_perf)
summary(filter_perf)
dim(filter_perf)


filter_perf$score_rsq <- filter_perf$filter_rsq + filter_perf$filter_has_predict
head(filter_perf[order(-filter_perf$score_rsq),])
best_score_rsq <- filter_perf[order(-filter_perf$score_rsq),][1,]

filter_perf$score_concord <- filter_perf$filter_concord + filter_perf$filter_has_predict
head(filter_perf[order(-filter_perf$score_concord),])
best_score_concord <- filter_perf[order(-filter_perf$score_concord),][1,]

#filter_perf$score_ave_rsq <- filter_perf$filterave_rsq + filter_perf$filterave_has_predict
filter_perf$score_ave_rsq <- filter_perf$filterave_rsq
head(filter_perf[order(-filter_perf$score_ave_rsq),])
best_score_ave_rsq <- filter_perf[order(-filter_perf$score_ave_rsq),][1,]

filter_perf$score_ave_concord <- filter_perf$filterave_concord * filter_perf$filterave_has_predict
head(filter_perf[order(-filter_perf$score_ave_concord),])
best_score_ave_concord <- filter_perf[order(-filter_perf$score_ave_concord),][1,]

#' Elbow plot
fig <- ggplot(filter_perf, aes(x=filter_has_predict, y=filter_rsq)) + 
  geom_point() + 
  geom_point(data=best_score_rsq, color="red", size=3) + 
  geom_line() + 
  theme_bw() + 
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=16)) 
print(fig)

fig <- ggplot(filter_perf, aes(x=filter_has_predict, y=filter_concord)) + 
  geom_point() + 
  geom_line() + 
  geom_point(data=best_score_concord, color="red", size=3) + 
  theme_bw() + 
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=16)) 
print(fig)

fig <- ggplot(filter_perf, aes(x=filterave_has_predict, y=filterave_rsq)) + 
  geom_point() + 
  geom_line() + 
  geom_point(data=best_score_ave_rsq, color="red", size=3) + 
  theme_bw() + 
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=16)) 
print(fig)

fig <- ggplot(filter_perf, aes(x=filterave_has_predict, y=filterave_concord)) + 
  geom_point() + 
  geom_line() + 
  geom_point(data=best_score_ave_concord, color="red", size=3) + 
  theme_bw() + 
  theme(axis.title = element_text(size=20),
        axis.text = element_text(size=16)) 
print(fig)

#' **concordance by Dataset**
#' -----------
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


realdnds <- read.table('/home/thuy/gitrepo/MutationPatterns/R/longshot/per_sample_dnds.csv', sep=",", header=TRUE)
realdnds$Subst.Act <- realdnds$N.Act+ realdnds$S.Act
summary(realdnds)
dim(realdnds)

dim(realdnds[realdnds$Subst.Act >= 2, ])
nrow(realdnds[realdnds$Subst.Act >= 2, ]) / nrow(realdnds)

dim(realdnds[realdnds$Subst.Act >= 5, ])
nrow(realdnds[realdnds$Subst.Act >= 5, ]) / nrow(realdnds)

realsubs <- read.table( '/home/thuy/gitrepo/MutationPatterns/R/longshot/timing/across_samples_sub.winsite.csv' , sep=",", header=TRUE)
realsubs$Subst.Act <- realsubs$N.Act+ realsubs$S.Act
summary(realsubs)
dim(realsubs)

dim(realsubs[realsubs$WinSiteSub >= 2, ])
nrow(realsubs[realsubs$WinSiteSub >= 2, ]) / nrow(realsubs)

dim(realsubs[realsubs$WinSiteSub >= 5, ])
nrow(realsubs[realsubs$WinSiteSub >= 5, ]) / nrow(realsubs)

#'ave dnds
#'
ave_dnds <- 
  
# rf predicted stuff
predicted <- read.table("lhs_regression_real_fix/umberjack_accuracy_predict.real.csv", sep=",", header=TRUE)
summary(predicted)
dim(predicted)


# rsquared of the true dnds vs expected dns
rSquared(y=predicted$dN_minus_dS.Act, resid=(predicted$dN_minus_dS.Act - predicted$dN_minus_dS.Exp))
# concordance  of the true dnds vs expected dns
epi.ccc(predicted$dN_minus_dS.Act, predicted$dN_minus_dS.Exp)$rho.c
# squared correlation.  ie true r2
cor(predicted$dN_minus_dS.Act, predicted$dN_minus_dS.Exp)^2

ave_predicted <- ddply(predicted, c("File", "CodonSite"), 
                       function(x) {
                         data.frame(Ave.dN_minus_dS.Act = weighted.mean(x$dN_minus_dS.Act, 
                                                                        x$UnambigCodonRate.Act * x$Reads.Act, na.rm=TRUE),
                                    Windows = length(x$Window_Start),
                                    dN_minus_dS.Exp = x$dN_minus_dS.Act[1])
                       })
head(ave_predicted)
dim(ave_predicted)

rSquared(y=ave_predicted$Ave.dN_minus_dS.Act, resid=(ave_predicted$Ave.dN_minus_dS.Act - ave_predicted$dN_minus_dS.Exp))
# concordance  of the true dnds vs expected dns
epi.ccc(ave_predicted$Ave.dN_minus_dS.Act, ave_predicted$dN_minus_dS.Exp)$rho.c
# squared correlation.  ie true r2
cor(ave_predicted$Ave.dN_minus_dS.Act, ave_predicted$dN_minus_dS.Exp)^2


# > IQR(gold$Scaled.dN.dS[gold$Sub > 0], na.rm=TRUE)
#[1] 3.28176122093

# rsquared of the true dnds vs expected dns when the random forest predicted an error less than 3.2
predicted_filt <- predicted[sqrt(predicted$oob_pred) <= 3.28176122093, ]
summary(predicted_filt)
dim(predicted_filt)

rSquared(y=predicted_filt$dN_minus_dS.Act, resid=(predicted_filt$dN_minus_dS.Act - predicted_filt$dN_minus_dS.Exp))
# concordance  of the true dnds vs expected dns when the random forest predicted an error less than 3.2
epi.ccc(predicted_filt$dN_minus_dS.Act, predicted_filt$dN_minus_dS.Exp)$rho.c
# true r2
cor(predicted_filt$dN_minus_dS.Act, predicted_filt$dN_minus_dS.Exp)^2

ave_predicted_filt <- ddply(predicted_filt, c("File", "CodonSite"), 
                       function(x) {
                         data.frame(Ave.dN_minus_dS.Act = weighted.mean(x$dN_minus_dS.Act, 
                                                                        x$UnambigCodonRate.Act * x$Reads.Act, na.rm=TRUE),
                                    dN_minus_dS.Exp = x$dN_minus_dS.Act[1])
                         })
dim(ave_predicted_filt)

rSquared(y=ave_predicted_filt$Ave.dN_minus_dS.Act, resid=(ave_predicted_filt$Ave.dN_minus_dS.Act - ave_predicted_filt$dN_minus_dS.Exp))
# concordance  of the true dnds vs expected dns
epi.ccc(ave_predicted_filt$Ave.dN_minus_dS.Act, ave_predicted_filt$dN_minus_dS.Exp)$rho.c
# squared correlation.  ie true r2
cor(ave_predicted_filt$Ave.dN_minus_dS.Act, ave_predicted_filt$dN_minus_dS.Exp)^2

fig <- ggplot(ave_predicted_filt,   aes(x=dN_minus_dS.Exp, Ave.dN_minus_dS.Act)) + 
  geom_abline(slope=1, color="red") + 
  geom_point(alpha=0.4) + 
  xlab("True dN-dS") + 
  ylab("Inferred dN-dS") + 
  theme_bw() + 
  theme(axis.text=element_text(size=16),
        axis.title=element_text(size=20)) + 
  geom_smooth(method="lm", se=FALSE)
print(fig)

filename <- "/home/thuy/gitrepo/MutationPatterns/tex/longshot/qc/cleaned_umberjack_vs_true.png"
if (!file.exists(filename)) {
  ggsave(filename=filename, plot=fig)
} else {
  warning(paste0("Not resaving", filename))
}
