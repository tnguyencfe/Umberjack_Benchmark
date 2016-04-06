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
library(xtable)
source('./load_all_sim_dnds.R')
source("plot_helper.R")

THESIS_DIR <- '../../MutationPatterns/tex'


latex_nice <-  function(name) {
  if (name == "BreakRatio.Act") {
    return ("Window Total Breakpoint Ratio")
  } else if (name == "Subst.Act") {
    return ("Window-Site Substitutions")
  } else if (name == "UnambigCodonRate.Act") {
    return ("Window-Site Unambiguous Codon Rate")
  } else if (name == "AADepth.Act") {
    return ("Window-Site Amino Acid Depth")
  } else if (name == "EntropyCodon.Act") {
    return ("Window-Site Codon Entropy")
  } else if (name == "UnknownPerCodon.Act") {
    return ("Window-Site Unknown Bases Per Read")
  } else if (name == "ErrPerCodon.Act") {
    return ("Window-Site Sequence Error Rate")
  } else if (name == "N.Act") {
    return ("Window-Site Nonsynonymous Substitutions")
  } else if (name == "S.Act") {
    return ("Window-Site Synonymous Substitutions")
  } else if (name == "TreeLen.Act") {
    return ("Window Tree Length")
  } else if (name == "TreeDepth.Act") {
    return ("Window Tree Depth (Subs/Site)")
  } else if (name == "TreeDistPerRead.Act") {
    return ("Normalized $\\overline{WRF}$")
  } else if (name == "TreeDist.Act") {
    return ("Window Weighted Robinson Foulds")
  } else if (name == "EntropyCodon.Exp") {
    return ("True Site Codon Entropy")
  } else if (name == "WinP_SameCodonFreq.Act") {
    return ("log10 P(Window-Site Codon Distro = True Site Codon Distro), Ave Across Window")
  } else if (name == "P_SameCodonFreq.Act") {
    return ("Codon Distribution P-value")
  } else if (name == "ResolvedPerSub.Act ") {
    return ("Fraction of Substitutions from Ambiguous Codons")
  } else if (name == "EN.Act") {
    return ("Normalized Window-Site Expected Nonsynonymous Substitutions")
  } else if (name == "ES.Act") {
    return ("Normalized Window-Site Expected Synonymous Substitutions")
  } else if (name == "N.Exp") {
    return ("True Site Nonsyn Subs")
  } else if (name == "S.Exp") {
    return ("True Site Syn Subs")
  } else if (name == "Cov.Act") {
    return ("Window Read Coverage Per Individual")
  } else if (name == "Window_Entropy.Act") {
    return ("Window-Site Codon Entropy, Ave Across Window")
  } else if (name == "Window_UnambigCodonRate.Act") {
    return ("Window-Site Unambig Codons Per Read, Ave Across Window")
  } else if (name == "Window_ErrPerCodon.Act") {
    return ("Window-Site Sequence Errors Per Read, Ave Across Window")
  } else if (name == "Window_Subst.Act") {
    return ("Site Substitutions Ave Across Window")
  } else if (name == "LOD_dNdS") {
    return ("log2(inferred window site dn/ds) - log2(expected site dn/ds)")
  } else if (name == "Dist_dn_minus_dS") {
    return ("(inferred window site dn-ds) - (expected site dn-ds)")
  } else if (name == "AbsLOD_dNdS") {
    return ("|log2(inferred window site dn/ds) - log2(expected site dn/ds)|")
  } else if (name == "AbsDist_dn_minus_dS") {
    return ("|(inferred window site dn-ds) - (expected site dn-ds)|")
  } else if (name == "SqDist_dn_minus_dS") {
    return ("[(inferred window site dn-ds) - (expected site dn-ds)]^2")
  } else if (name == "WinAbsLOD_dNdS") {
    return ("|log2(inferred window site dn/ds) - log2(expected site dn/ds)| \n Ave Across Window")
  } else if (name == "WinAbsDist_dn_minus_dS") {
    return ("|(inferred window site dn-ds) - (expected site dn-ds)|\nAve Across Window")
  } else if (name == "WinSqDist_dn_minus_dS") {
    return ("[(inferred window site dn-ds) - (expected site dn-ds)] ^2\n Ave Across Window")
  } else if (name == "TreeLenPerRead.Act") {
    return ("Normalized Window Tree Length")
  } else if (name == "PolytomyPerRead.Act") {
    return ("Normalized Polytomies")
  } else {
    return (name)
  }
}


# Per Window-Site data
dnds_filename <- "../simulations/out/collate_all.recombo.csv"
if (exists("dnds_filename")) {
  recombo_dnds <- get_all_sim_dnds(dnds_filename)  
} else {
  recombo_dnds <- get_all_sim_dnds()
}

# Per Window data
recombo_window <-  get_window_sim_dnds(dnds=recombo_dnds)
recombo_rate <- read.table('/home/thuy/gitrepo/Umberjack_Benchmark/sim_config/sim_args.recombo.tsv', header=TRUE, sep="\t")
recombo_rate$rate <- recombo_rate$NumBreakpoints / (recombo_rate$CodonSites * recombo_rate$Generations)
recombo_rate$BreakPerCodon <- recombo_rate$CodonSites / recombo_rate$NumBreakpoints

recombo_rate$rateFactor <- paste(format((signif(recombo_rate$rate, 2)), scientific=TRUE), 
                                 format((signif(recombo_rate$BreakPerCodon, 1)), scientific=FALSE),
                                 sep=", ")
recombo_rate$rateFactor[is.infinite(recombo_rate$BreakPerCodon)] <- "0.0, "

rate_str <- unique(recombo_rate[, c("rate", "BreakPerCodon", "rateFactor")])
rate_str <- rate_str[order(rate_str$rate), ]

recombo_rate$rateFactor <- factor(recombo_rate$rateFactor, levels=rate_str$rateFactor)

head(recombo_rate)
summary(recombo_rate)
dim(recombo_rate)

recombo_window <- merge(x=recombo_window, 
                        y=subset(recombo_rate, select=c(Name, rate, rateFactor)),
                        by.x=c("File"), by.y=c("Name"),
                        all.x=TRUE, all.y=FALSE)
head(recombo_window)
dim(recombo_window)
summary(recombo_window)





#' Output simluated data params as latex files
#' ===========
sim_args <- read.table("/home/thuy/gitrepo/Umberjack_Benchmark/sim_config/sim_args.tsv", header=TRUE, sep="\t")
#Name  PopSize	CodonSites	ART_Profile	Seed	SelectionRate	FragLenStd	FragLenAve	Cover	NumBreakpoints	Generations	TreeLen	WindowSize	MinWinWidth	MinWinDepth	MinQual	recomborate
sim_args$recomborate <- sim_args$NumBreakpoints / (sim_args$CodonSites * sim_args$Generations)

sim_args <- adply(.data=sim_args,
                  .margins=1,
                  .fun=function(row) {
  n_treelens <- as.numeric(unlist(strsplit(as.character(row$TreeLen), split=",")))
  mut_rates <- n_treelens * row$SelectionRate/row$Generations
  #s_mut_rates <- paste0(format(mut_rates, digits=2, scientific=TRUE), collapse=", ")
  
  return (data.frame(MutRate1=mut_rates[1], MutRate2=mut_rates[2], MutRate3=mut_rates[3]))
})

summary(sim_args)
head(sim_args)
dim(sim_args)


# PopSize
# Cover
# NumBreakpoints
# Generations
# TreeLen
# FragLenAve
# Don't use kable since we cant do per-column formatting
recombo_rate$MutRate <- recombo_rate$TreeLen / (recombo_rate$Generations * recombo_rate$PopSize)
cols <- c("Name", "rate")
nice_cols <- c("Name", "Recombination Rate")
tab <- xtable(recombo_rate[, cols], 
       caption=paste0("Parameters for generating simulated datasets to predict Umberjack accuracy under varying recombination.  ", 
                      "Each row represents a simulated population and its simulated paired-end MiSeq sequence library.  ",
                      "Recombination rate units in recombinations/bp/generation.  ",
                      "Mutation rate units in mutations/bp/generation.  ",
                      "Read coverage is per extent individual in the population.  ",
                      "Sequencing fragment size in bp.  ",
                      "Genome size = 930bp. Selection rate = 0.01/generation. ",                      
                      "Mutation Rate = 4e-5 mutations/bp/generation.  ", 
                      "Generations = 5000.  ", 
                      "Extent population size = 100.  ", 
                      "2x250bp MiSeq paired-end reads.  ", 
                      "Read coverage per individual = 2x.  ", 
                      "Sequencing fragment size mean = 375bp.  ",
                      "Sequencing fragment size standard deviation = 75bp.  ",
                      "Umberjack Configuration b = {Window size = 300bp, Min window width coverage = 0.875, Min window read depth = 10, Min phred quality score = 20}."),
       label="tab:sim_args_recombo",
       align=c(rep("c|", length(cols)), "c"),
       display=c("s", "s", "e"),
       )
names(tab) <- nice_cols
filename <- paste0(THESIS_DIR, "/umberjack/sim_args_recombo.tex")
if (file.exists(filename)) {
  warning(paste0("Not regenerating", filename))
} else {
  write(print(tab, include.rownames=FALSE),  file=filename)  
}



#' Print out latin hypercube values
#FieldName  Min	Max	Is_Int
sim_args_range <- read.table("/home/thuy/gitrepo/Umberjack_Benchmark/sim_config/sim_args_range.csv", sep=",", header=TRUE)
summary(sim_args_range)
head(sim_args_range)
dim(sim_args_range)


#' Print out recombo sim args
#' 
cols <- c("Generations", "MutRate1", "MutRate2", "MutRate3", "recomborate", "NumBreakpoints", "Cover", "FragLenAve")
nice_cols <- c("Generations", "Mutation Rate 1", "Mutation Rate 2", "Mutation Rate 3", 
               "Recombination Rate", "Breakpoints", "Coverage", "Mean Fragment Size")
tab <- xtable(sim_args[, cols], 
              caption=paste0("Parameters for generating simulated datasets to predict Umberjack accuracy.  ", 
                             "Each row represents a simulated population and its simulated paired-end MiSeq sequence library.  ",
                             "Recombination rate units in recombinations/bp/generation.  ",
                             "Mutation rate units in mutations/bp/generation.  ",
                             "Read coverage is per extent individual in the population.  ",
                             "Sequencing fragment size in bp.  ",
                             "Genome size = 900bp. Selection rate = 0.01/generation. ",                      
                             "Extent population size = 1000.  ", 
                             "2x250bp MiSeq paired-end reads.  ", 
                             "Sequencing fragment size standard deviation = 100bp.  ",
                             "Umberjack Configuration a = {Window size = 150bp, Min window width coverage = 0.7, Min window read depth = 10, Min phred quality score = 15}, ",
                             "Umberjack Configuration b = {Window size = 300bp, Min window width coverage = 0.875, Min window read depth = 10, Min phred quality score = 20}."),
              label="tab:sim_args_umberjack_accuracy",
              align=c(rep("c|", length(cols)), "c"),
              display=c("s", "d", "e", "e", "e", "e", "d", "d", "d"),
)
names(tab) <- nice_cols
filename <- paste0(THESIS_DIR, "/umberjack/sim_args.tex")
if (file.exists(filename)) {
  warning(paste0("Not regenerating", filename))
} else {
  write(print(tab, include.rownames=FALSE),  file=filename)  
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

# # Plot response vs single variable, color points using gradient
# plot_gradient <- function(data, resp_colname, var_colname, color_colnames=NULL) {  
#   figs <- sapply(color_colnames, 
#                  function(color_colname) {                   
#                    fig <- ggplot(data, aes_string(x=var_colname, y=resp_colname)) + 
#                      xlab(nice(var_colname)) + 
#                      ylab(nice(resp_colname)) +                     
#                      theme_bw() + 
#                      geom_point(aes_string(color=color_colname), alpha=0.7, na.rm=TRUE) +                      
#                      geom_smooth(method="lm", se=FALSE, color="black", size=2) +                      
#                      scale_colour_gradient2(low="blue", mid="grey", high="red", midpoint=mean(data[, color_colname], na.rm=TRUE)) +
#                      ggtitle("Inaccuracy Vs Covariate")
#                    print(fig)
#                  })
# }
# 
# #+ fig.width=10
# plot_gradient(data=recombo_window, resp_colname="WinSqDist_dn_minus_dS", 
#               var_colname="TreeDistPerRead.Act",
#               #color_colnames=c("LogTreeLen.Act", "LogPolytomy.Act", "LogPolytomyPerTreeLen.Act", "PolytomyPerRead.Act", "LogTreeLenPerRead.Act",
#                                #"LogReads.Act", "LogWindow_Subst.Act", "WinP_SameCodonFreq.Act")
#               color_colnames=c("LogTreeLen.Act"))



#' Show effect of recombination alone on umberjack inaccuracy, highlighting effects of low diversity
#' 
#+ fig.width=10, fig.height=7
fig <- ggplot(recombo_window, aes(x=TreeDistPerRead.Act, y=WinSqDist_dn_minus_dS)) + 
  xlab(expression(paste("Normalized ", bar("WRF")))) + 
  ylab(expression(paste("Window Average ", Delta))) + 
  theme_bw(base_size=12) + 
  geom_point(aes(color=recombo_window$TreeLen.Act)) +                      
  geom_smooth(method="lm", se=FALSE, color="black", size=2) +                      
  scale_colour_gradient2(name="Tree\nLength", low="blue", mid="darkgrey", high="red", 
                         midpoint=quantile(recombo_window$TreeLen.Act, probs=0.5, na.rm=TRUE))  + 
  theme(axis.title=element_text(size=20), axis.text=element_text(size=20),
      legend.title=element_text(size=20),
      legend.text=element_text(size=16),
      legend.position=c(0.8, 0.7))
print(fig)

picname <- paste0(THESIS_DIR, "/umberjack/delta_v_recombo.png")
if (file.exists(picname)) {
  warning(paste0("Not replacing ", picname))
} else {
  ggsave(filename=picname, plot=fig, width=10, height=7, units="in")  
} 


#' Show effect of recombination alone on umberjack inaccuracy, highlighting effects of recombination rate
#' 
head(recombo_window)
summary(recombo_window)
dim(recombo_window)

#+ fig.width=10, fig.height=7
fig <- ggplot(recombo_window, aes(x=TreeDistPerRead.Act, y=WinSqDist_dn_minus_dS)) + 
  xlab(expression(paste("Normalized ", bar("WRF")))) + 
  ylab(expression(paste("Window Average ", Delta))) + 
  theme_bw(base_size=12) + 
  geom_point(aes(color=rateFactor), alpha=0.7, size=2) +                      
  geom_smooth(method="lm", se=FALSE, color="black", size=2) +          
  scale_color_discrete(name=expression(paste(rho, ", Interval"))) +
  theme(axis.title=element_text(size=20), axis.text=element_text(size=20),
        legend.title=element_text(size=20),
        legend.text=element_text(size=16),
        legend.position=c(0.8, 0.7))
print(fig)

picname <- paste0(THESIS_DIR, "/umberjack/delta_v_recombo_showRate.png")
if (file.exists(picname)) {
  warning(paste0("Not replacing ", picname))
} else {
  ggsave(filename=picname, plot=fig, width=10, height=7, units="in")  
} 


#' Show effect of recombination breakpoints on robinson foulds
fig <- ggplot(recombo_window, aes(x=Window_Breaks, y=TreeDistPerRead.Act)) + 
  ylab(expression(paste("Normalized ", bar("WRF")))) + 
  xlab("Breakpoints in Window") +
  theme_bw(base_size=12) + 
  geom_point(alpha=0.7, size=2) +                      
  geom_smooth(method="lm", se=FALSE, color="blue", size=2) +          
  scale_colour_gradient2(name="Tree\nLength", low="blue", mid="darkgrey", high="red", 
                         midpoint=quantile(recombo_window$TreeLen.Act, probs=0.5, na.rm=TRUE))  + 
  theme(axis.title=element_text(size=20), axis.text=element_text(size=20),
        legend.title=element_text(size=20),
        legend.text=element_text(size=16),
        legend.position=c(0.2, 0.7))
print(fig)

picname <- paste0(THESIS_DIR, "/umberjack/breaks_v_robinsonfoulds.png")
if (file.exists(picname)) {
  warning(paste0("Not replacing ", picname))
} else {
  ggsave(filename=picname, plot=fig, width=10, height=7, units="in")  
} 

###############  ALL SIMULATED DATA FOR RANDOM FOREST

dnds_filename <- "../simulations/out/collate_all.treedist.csv"
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
#+ fig.width=2.5, fig.height=2.5
dndsrange <- outlier_range(dnds$SqDist_dn_minus_dS)
load("./lhs_regression_fix/allfitDist_gamma_glm.RData")


for (predictor in LM_COVAR_NAMES) {
  fig <- ggplot(dnds, aes_string(x=predictor, y="SqDist_dn_minus_dS")) + 
    xlab(nice(predictor)) + 
    ylab(expression(Delta)) + 
    theme_bw(base_size=12) + 
    geom_point(alpha=0.1, shape=1) +                               
    #geom_smooth(method="lm", se=FALSE, color="blue") +                      
    #geom_abline(intercept=allfitDist$coeff["(Intercept)"],
    #            slope=allfitDist$coeff[predictor], color="blue") + 
    scale_y_continuous(limits=dndsrange)
    #scale_colour_gradient2(name="Tree\nLength", low="blue", mid="darkgrey", high="red", 
    #                       midpoint=quantile(recombo_window$TreeLen.Act, probs=0.5, na.rm=TRUE)) + 
#     theme(axis.title=element_text(size=20), axis.text=element_text(size=20),
#           legend.title=element_text(size=20))
  print(fig)
  
  picname <- paste0(THESIS_DIR, "/umberjack/univar_blowup_delta_v_", predictor, ".png")
  if (file.exists(picname)) {
    warning(paste0("Not replacing ", picname))
  } else {
    ggsave(filename=picname, plot=fig, width=2.5, height=2.5, units="in")  
  } 

}


#' Plot the univariate regression plots, with outliers intact
#+ fig.width=3, fig.height=2.1
for (predictor in LM_COVAR_NAMES) {
  
  fig <- ggplot(dnds, aes_string(x=predictor, y="SqDist_dn_minus_dS")) + 
    xlab(nice(predictor)) + 
    ylab(expression(Delta)) + 
    theme_bw(base_size=12) + 
    geom_point(alpha=0.1)       
    #geom_smooth(method="lm", se=FALSE, color="blue")                  
    #scale_y_continuous(limits=outlier(dnds$SqDist_dn_minus_dS)) + 
    #scale_colour_gradient2(name="Tree\nLength", low="blue", mid="darkgrey", high="red", 
    #                       midpoint=quantile(recombo_window$TreeLen.Act, probs=0.5, na.rm=TRUE)) + 
#     theme(axis.title=element_text(size=20), axis.text=element_text(size=20),
#           legend.title=element_text(size=20))
  print(fig)
  
  picname <- paste0(THESIS_DIR, "/umberjack/univar_delta_v_", predictor, ".png")
  if (file.exists(picname)) {
    warning(paste0("Not replacing ", picname))
  } else {
    ggsave(filename=picname, plot=fig, width=2.5, height=2.5, units="in")  
  } 
  
}



#' Find the spearman's correlation, distance correlation between individual predictor and the response
#' 
univar_rank <- data.frame(predictor=as.character(LM_COVAR_NAMES))
univar_rank <- adply(.data=univar_rank,
                     .margins=1,
                     .fun=function(row) {
                       predictor <- as.character(row$predictor)
                       result <- cor.test(dnds$SqDist_dn_minus_dS, 
                                          dnds[, predictor], method="spearman", na.action=na.exclude, exact=FALSE)
                       return (data.frame(spear_cor=result$estimate, p.value=result$p.value))
                     })
univar_rank$abs_spear_cor <- abs(univar_rank$spear_cor)
univar_rank <- univar_rank[order(-univar_rank$abs_spear_cor), ]
univar_rank$latex_name <- sapply(as.character(univar_rank$predictor), latex_nice)

cols <- c("latex_name", "spear_cor", "p.value")
nice_cols <- c("Feature", "Corr", "P-value")
tab <- xtable(univar_rank[, cols]
       , caption="Spearman Correlation Between Feature and $\\Delta$"
       , label="tab:univar_corr_delta"
       , align=c("c", "l|", "c|", "c")
       , display=c("s", "e", "e", "e")
       , digits=2)
names(tab) <- nice_cols
filename <- paste0(THESIS_DIR, "/umberjack/univar_corr_delta.tex")
if (file.exists(filename)) {
  warning(paste0("Not regenerating", filename))
} else {
  write(print(tab, include.rownames=FALSE, sanitize.text.function=identity), file=filename)  
}


#' Print out variable importance table for accuracy random forest
imp <- read.table("./lhs_regression_fix/importance.csv", sep=",", header=TRUE)
head(imp)
summary(imp)
dim(imp)

#imp <- imp[imp$Feature %in% LM_COVAR_NAMES, ]
imp$p.value <- 2*pnorm(-abs(imp$IncMSE))
imp$latex_feature <- sapply(imp$Feature, latex_nice)

cols <- c("latex_feature", "IncMSE", "p.value")
nice_cols <- c("Feature", "Feature Importance", "Approx p-value")
tab <- xtable(imp[, cols]
              , caption=paste0("Feature Ranking for Final Random Forest Predicting Umberjack Accuracy.")
              , label="tab:randomforest_feature_imp"
              , align=c("c", "l|", "c|", "c")
              , display=c("s", "s", "fg", "e")
              , digits=2)
names(tab) <- nice_cols
filename <- paste0(THESIS_DIR, "/umberjack/randomforest_feature_imp.tex")
if (!file.exists(filename)) {
  write(print(tab, include.rownames=FALSE), file=filename)  
} else {
  warning(paste0("Not regenerating ", filename))
}


