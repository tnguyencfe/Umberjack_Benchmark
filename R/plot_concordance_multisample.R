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

NUM_RESP_NAMES <- c("LOD_dNdS", "Dist_dn_minus_dS", "AbsLOD_dNdS", "AbsDist_dn_minus_dS")
CAT_RESP_NAMES <- c("CrapLOD", "CrapDist")
COVAR_NAMES <- colnames(dnds[sapply(dnds,is.numeric)])[!colnames(dnds[sapply(dnds,is.numeric)]) %in% NUM_RESP_NAMES]
CAT_COVAR_NAMES <- c("IsLowSubst.Act")
LM_COVAR_NAMES <- c(CAT_COVAR_NAMES, 
                    COVAR_NAMES[!(COVAR_NAMES %in% c("dNdS.Act", "dN_minus_dS.Act", "dN_minus_dS.Exp", 
                                                     # In separate analysis, conservation and entropy are highly correlated.
                                                     # When we use speedglm, it bugs out after it removes highly correlated variables.
                                                     # So we do it for them.
                                                     "ConserveTrueBase.Act", "ConserveTrueBase.Exp", "Window_Conserve.Act",
                                                     "UnambigCodonRate.Act", "Window_UnambigCodonRate.Act",
                                                     # These are highly correlated with N, S
                                                     "Subst.Act", "Subst.Exp",
                                                     "EN.Exp", "ES.Exp", "EN.Act", "ES.Act",
                                                     "Window_Start", "Window_End", "CodonSite", "Reads.Act"
                    )
                    )])

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
#gooddnds <- dnds[!(dnds$S.Act < 1 & dnds$S.Act > 0), ]
gooddnds$CovBin <- cut(gooddnds$Cov.Act, breaks=0:5, right=TRUE)
gooddnds$WinSize <- as.factor(gooddnds$Window_End - gooddnds$Window_Start + 1)
gooddnds$BreadthThresh <- as.factor(gooddnds$BreadthThresh)
gooddnds$DepthThresh <- as.factor(gooddnds$DepthThresh)
gooddnds$ErrFree <- FALSE
gooddnds$ErrFree[grep("errFree", gooddnds$File)] <- TRUE
gooddnds$ErrFree <- as.factor(gooddnds$ErrFree)

summary(gooddnds)
dim(gooddnds)



# concord <- epi.ccc(gooddnds[gooddnds$File == "/home/thuy/gitrepo/Umberjack_Benchmark/simulations/out/small.cov2.indiv1000.codon400.bwa.rand/consensus/window300.breadth0.875.depth10.0/collate_dnds.csv",]$dNdS.Act, 
#                    gooddnds[gooddnds$File =="/home/thuy/gitrepo/Umberjack_Benchmark/simulations/out/small.cov2.indiv1000.codon400.bwa.rand/consensus/window300.breadth0.875.depth10.0/collate_dnds.csv",]$dNdS.Exp)
# print(concord$rho.c$est)


concord <- epi.ccc(gooddnds$dNdS.Act, 
                   gooddnds$dNdS.Exp)
print(concord$rho.c$est)
#' **Concordance for dN/dS when only good window-sites considered = `r concord$rho.c$est`**
#' 

concord <- epi.ccc(gooddnds$dN_minus_dS.Act, 
                   gooddnds$dN_minus_dS.Exp)
print(concord$rho.c$est)
#' **Concordance for dN-dS when only good window-sites considered = `r concord$rho.c$est`**
#' 

aves <- ddply(.data=gooddnds, .variables=c("File", "CodonSite", "ErrFree"), .fun=function(x) {
  if (length(unique(x$dNdS.Exp)) > 1 | length(unique(x$dN_minus_dS.Exp)) > 1) {
    print(unique(x$dNdS.Exp))
    print(unique(x$dN_minus_dS.Exp))
    stop("Too many dnds.exp or dn-ds.exp")
  }
  data.frame(Ave_dNdS.Act=weighted.mean(x$dNdS.Act, x$UnambigCodonRate*x$Reads.Act, na.rm=TRUE),  # weighted by unambiguous codon depth
             Ave_dN_minus_dS.Act=weighted.mean(x$dN_minus_dS.Act, x$UnambigCodonRate*x$Reads.Act, na.rm=TRUE),
             dNdS.Exp=x$dNdS.Exp[1],
             dN_minus_dS.Exp=x$dN_minus_dS.Exp[1])
})
summary(aves)
head(aves)
dim(aves)



concord <- epi.ccc(aves[aves$ErrFree == FALSE, ]$Ave_dNdS.Act, aves[aves$ErrFree == FALSE, ]$dNdS.Exp)
print(concord$rho.c$est)

concord <- epi.ccc(aves[aves$ErrFree == TRUE, ]$Ave_dN_minus_dS.Act, aves[aves$ErrFree == TRUE, ]$dN_minus_dS.Exp)
print(concord$rho.c$est)

concord <- epi.ccc(aves[aves$ErrFree == TRUE, ]$Ave_dNdS.Act, aves[aves$ErrFree == TRUE, ]$dNdS.Exp)
print(concord$rho.c$est)
#' **Concordance for Ave dN/dS when only good window-sites considered = `r concord$rho.c$est`**
#' 

concord <- epi.ccc(aves[aves$ErrFree == TRUE, ]$Ave_dN_minus_dS.Act, aves[aves$ErrFree == TRUE, ]$dN_minus_dS.Exp)
print(concord$rho.c$est)
#' **Concordance for Ave dN-dS when only good window-sites considered, err free = `r concord$rho.c$est`**
#' 
#' 

# concord <- epi.ccc(aves[aves$File=="/home/thuy/gitrepo/Umberjack_Benchmark/simulations/out/small.cov2.indiv1000.codon400.bwa.rand/consensus/window300.breadth0.875.depth10.0/collate_dnds.csv" & aves$ErrFree == FALSE, ]$Ave_dNdS.Act, 
#                    aves[aves$File=="/home/thuy/gitrepo/Umberjack_Benchmark/simulations/out/small.cov2.indiv1000.codon400.bwa.rand/consensus/window300.breadth0.875.depth10.0/collate_dnds.csv" & aves$ErrFree == FALSE, ]$dNdS.Exp)
# print(concord$rho.c$est)
# 
# concord <- epi.ccc(aves[aves$File=="/home/thuy/gitrepo/Umberjack_Benchmark/simulations/out/small.cov2.indiv1000.codon1600/consensus/window300.breadth0.875.depth10.0/collate_dnds.csv" & aves$ErrFree == FALSE, ]$Ave_dNdS.Act, 
#                    aves[aves$File=="/home/thuy/gitrepo/Umberjack_Benchmark/simulations/out/small.cov2.indiv1000.codon1600/consensus/window300.breadth0.875.depth10.0/collate_dnds.csv" & aves$ErrFree == FALSE, ]$dNdS.Exp)
# print(concord$rho.c$est)
# 
# orig_act_dnds <- read.table("/home/thuy/gitrepo/SlidingWindow/test/simulations/out_hyphyfix/small.cov2.indiv1k.codon400.bwa.rand/consensus/window300/actual_dnds_by_site.csv",
#                             header=TRUE, sep=",")
# head(orig_act_dnds)
# dim(orig_act_dnds)
# summary(orig_act_dnds)
# 
# orig_act_dnds <- merge(x=orig_act_dnds, y=aves[aves$File=="/home/thuy/gitrepo/Umberjack_Benchmark/simulations/out/small.cov2.indiv1000.codon400.bwa.rand/consensus/window300.breadth0.875.depth10.0/collate_dnds.csv" & aves$ErrFree == FALSE, ],
#                                                by.x="Site", by.y="CodonSite", all=TRUE)
# 
# head(orig_act_dnds)
# dim(orig_act_dnds)
# 
# sum(is.na(orig_act_dnds$Ave_dNdS.Act))
# sum(is.na(orig_act_dnds$dNdSWeightByReadsNoLowSyn))
# 
# sum(abs(orig_act_dnds$Ave_dNdS.Act - orig_act_dnds$dNdSWeightByReadsNoLowSyn) < 0.1, na.rm=TRUE)
# sum(!is.na(orig_act_dnds$Ave_dNdS.Act) & !is.na(orig_act_dnds$dNdSWeightByReadsNoLowSyn))
# head(orig_act_dnds[!is.na(orig_act_dnds$Ave_dNdS.Act) & !is.na(orig_act_dnds$dNdSWeightByReadsNoLowSyn) &
#                 abs(orig_act_dnds$Ave_dNdS.Act - orig_act_dnds$dNdSWeightByReadsNoLowSyn) > 0.1,])
# 
# 
# 
# concord <- epi.ccc(orig_act_dnds$dNdSWeightByReadsNoLowSyn, orig_act_dnds$dNdS.Exp)
# print(concord$rho.c$est)
# 
# concord <- epi.ccc(orig_act_dnds$Ave_dNdS.Act, orig_act_dnds$dNdS.Exp)
# print(concord$rho.c$est)

#' **Concordance for Ave dN/dS when only good window-sites considered = `r concord$rho.c$est`**
#' 

concord <- epi.ccc(aves[aves$ErrFree == FALSE, ]$Ave_dN_minus_dS.Act, aves[aves$ErrFree == FALSE, ]$dN_minus_dS.Exp)
print(concord$rho.c$est)
#' **Concordance for Ave dN-dS when only good window-sites considered, typical error `r concord$rho.c$est`**
#' 
#' 

#' **Concordance by window size, breadth, depth**
#'
concord_by_windowTrait <- ddply(.data=gooddnds, .variables=c("WinSize", "BreadthThresh", "DepthThresh", "CovBin", "ErrFree"), 
                        .fun=function(x) {                            
                          data.frame(concord_dN_minus_dS=epi.ccc(x$dN_minus_dS.Act, x$dN_minus_dS.Exp)$rho.c$est,
                                     concord_dNdS=epi.ccc(x$dNdS.Act, x$dNdS.Exp)$rho.c$est,
                                     total=nrow(x))
                        })
summary(concord_by_windowTrait)
print(concord_by_windowTrait)

# Plot concordance by window trait
fig <- ggplot(concord_by_windowTrait, aes(x=DepthThresh, y=concord_dN_minus_dS, color=CovBin)) + 
  geom_point() + 
  geom_line() + 
  xlab("\nDepth Threshold as Fraction of Population Size") + 
  ylab("Concorance for dN-dS\n") + 
  ggtitle("Concordance by Window Trait") + 
  facet_grid(WinSize~BreadthThresh)
print(fig)


fig <- ggplot(concord_by_windowTrait, aes(x=DepthThresh, y=concord_dNdS, color=CovBin)) + 
  geom_point() + 
  geom_line() + 
  xlab("\nDepth Threshold as Fraction of Population Size") + 
  ylab("Concorance for dN/dS\n") + 
  ggtitle("Concordance by Window Trait") + 
  facet_grid(WinSize+ErrFree~BreadthThresh)
print(fig)

# Plot total windows by trait
fig <- ggplot(concord_by_windowTrait, aes(x=DepthThresh, y=total, color=CovBin)) + 
  geom_point() +   
  xlab("\nDepth Threshold as Fraction of Population Size") + 
  ylab("Total Windows\n") + 
  ggtitle("Total Windows by Trait") + 
  facet_grid(WinSize+ErrFree~BreadthThresh)
print(fig)



#' Plot wndow una mbig codon by site
fig <- ggplot(unique(gooddnds[gooddnds$File=="/home/thuy/gitrepo/Umberjack_Benchmark/simulations/out/small.cov2.indiv1000.codon400.bwa.rand/consensus/window300.breadth0.875.depth10.0/collate_dnds.csv", 
                              c("Window_Start", "ErrFree", "Window_UnambigCodonRate.Act"),]), 
              aes(x=Window_Start, y=Window_UnambigCodonRate.Act, color=ErrFree)) + 
  geom_point(alpha=0.5, shape=1) + 
  geom_smooth(method="loess") + 
  xlab("\nsite") + 
  ylab("Window Unambig Codon Rate\n") + 
  ggtitle("Window Unambi Codon Across Sites")
print(fig)

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


fig <- ggplot(gooddnds, aes(x=EntropyTrueBase.Exp, y=LOD_dNdS, color=Cov.Act)) + 
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


fig <- ggplot(gooddnds, aes(x=EntropyTrueBase.Exp, y=AbsLOD_dNdS, color=CovBin)) + 
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


fig <- ggplot(gooddnds, aes(x= Window_Entropy.Act, y=AbsLOD_dNdS, color=CovBin)) + 
  #geom_point(alpha=0.5, shape=1) + 
  geom_smooth(method="loess") + 
  xlab("\nWindow Entropy") + 
  ylab("|Log(Inferred dN/dS) - Log(Expected dN/dS)|\n") + 
  ggtitle("Log Odds Difference dN/dS By Window Entropy and Coverage")
print(fig)


fig <- ggplot(gooddnds, aes(x=as.factor(WinSize), y=AbsLOD_dNdS, color=CovBin)) + 
  #geom_point(alpha=0.5, shape=1) + 
  #geom_smooth(method="loess") + 
  geom_boxplot() + 
  xlab("\nWindow Size") + 
  ylab("|Log(Inferred dN/dS) - Log(Expected dN/dS)|\n") + 
  ggtitle("Log Odds Difference dN/dS By Window Size and Coverage")
print(fig)


fig <- ggplot(gooddnds, 
              aes(x=Window_UnambigCodonRate.Act, y=AbsLOD_dNdS, color=CovBin)) + 
  geom_point(alpha=0.5, shape=1) + 
  geom_smooth(method="loess", na.rm=TRUE) + 
  #scale_x_log10() + 
  xlab("\nWindow Unambiguous Codon Rate") + 
  ylab("|Log(Inferred dN/dS) - Log(Expected dN/dS)|\n") + 
  ggtitle("Log Odds Difference dN/dS By Window Size and Coverage")
print(fig)
