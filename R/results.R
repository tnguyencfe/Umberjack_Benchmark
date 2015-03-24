library(limma)
library(ggplot2)
library(reshape2)
library(knitr)
library(epiR)
library(plyr)


# results <- read.table("./results.small.cov2x.indiv1k.codon400.bwa.rand.csv", header=TRUE, sep=",")
# summary(results)
# head(results)
# 
# results$MutationFactor <- factor(results$Mutation)
# results$QualityFactor <- paste0(results$Quality, " Read Quality")
# 
# ggplot(results, aes(x=Mutation, y=Condordance)) + 
#   geom_point() + 
#   ylab("Concordance\n") + 
#   xlab("\nSubstitutions/Site") + 
#   #scale_x_continuous(breaks=seq(c(0.004, 0.06, 0.02)) + 
#   # scale_y_continuous(breaks=seq(0, 1.25, 0.25)) + 
#    #scale_linetype_manual(values = c(rep("solid", 10), rep("dashed", 6))) +
#    #scale_color_manual(breaks=c("F1","F2","F3"),
#    #theme(axis.title=element_text(size=32), axis.text=element_text(size=16) ) + 
#    theme(legend.position="bottom") + 
#    scale_color_discrete(name="Genome Coverage") + 
#    theme(strip.text.x = element_text(size = 14), axis.title=element_text(size=14), axis.text=element_text(size=14)) + 
#    facet_wrap(~QualityFactor)
                       
ACTUAL_DNDS_FILENAME<-"/home/thuy/gitrepo/SlidingWindow/test/simulations/out_hyphyfix/small.cov2.indiv1k.codon400.bwa.rand/consensus/window300/actual_dnds_by_site.csv"
ACTUAL_ERRFREE_DNDS_FILENAME<-"/home/thuy/gitrepo/SlidingWindow/test/simulations/out_hyphyfix/small.cov2.indiv1k.codon400.bwa.rand/consensus/window300.errFree/actual_dnds_by_site.csv"
EXPECTED_DNDS_FILENAME<-"/home/thuy/gitrepo/SlidingWindow/test/simulations/data/small.cov2.indiv1k.codon400.bwa.rand/mixed/small.cov2.indiv1k.codon400.bwa.rand.mixed.dnds.tsv"
INDELIBLE_DNDS_FILENAME<-"/home/thuy/gitrepo/SlidingWindow/test/simulations/data/small.cov2.indiv1k.codon400.bwa.rand/mixed/small.cov2.indiv1k.codon400.bwa.rand.mixed.rates.csv"

# Read in Indelible's intended dN/dS
# Cols: Site,Interval,Scaling_factor,Rate_class,Omega
indelible_dnds <- read.table(INDELIBLE_DNDS_FILENAME, header=TRUE, sep=",")
summary(indelible_dnds)
head(indelible_dnds)


# Cols:  Ref  Site  aveDnDs  dNdSWeightBySubst  dN_minus_dS	Windows	Codons	NonSyn	Syn	Subst	dNdSWeightByReads	multisiteAvedNdS	multisitedNdSWeightBySubst	simpleDnDs  dNdSWeightByReadsNoLowSyn
 actual_dnds <- read.table(ACTUAL_DNDS_FILENAME, header=TRUE, na.strings="None", comment.char = "#", sep=",")
 dim(actual_dnds)
 head(actual_dnds)
 str(actual_dnds)
 summary(actual_dnds)

actual_errFree_dnds <- read.table(ACTUAL_ERRFREE_DNDS_FILENAME, header=TRUE, na.strings="None", comment.char = "#", sep=",")
dim(actual_errFree_dnds)
head(actual_errFree_dnds)
str(actual_errFree_dnds)
summary(actual_errFree_dnds)
 
 # Cols: Observed S Changes  Observed NS Changes	E[S Sites]	E[NS Sites]	Observed S. Prop.	P{S}	dS	dN	dN-dS	P{S leq. observed}	P{S geq. observed}	Scaled dN-dS	dn/ds
 # Parse the expected dnds filename for it start and end nucleotide positions (1-based)
 
 expected_dnds <- read.table(EXPECTED_DNDS_FILENAME, header=TRUE, sep="\t")  # Site  Interval	Scaling_factor	Rate_class	Omega
expected_dnds_start_nuc_pos <- 1
 expected_dnds_start_codon <- (expected_dnds_start_nuc_pos %/% 3) + 1
 expected_dnds$Site <- as.numeric(rownames(expected_dnds)) + expected_dnds_start_codon -1
 
 expected_dnds$Omega <- expected_dnds$dN/expected_dnds$dS
 expected_dnds$Omega[expected_dnds$dS == 0] <- NA
expected_dnds$Omega[expected_dnds$Observed.S.Changes < 1 | expected_dnds$Observed.NS.Changes < 1] <- NA
 dim(expected_dnds)
 head(expected_dnds)
 str(expected_dnds)
 summary(expected_dnds)
 
 
  
 # check consistency
 all(expected_dnds$Site == actual_dnds$Site)
 
 
 

fullDat <- data.frame(Site=expected_dnds$Site,
                      ActualTypical=actual_dnds$dNdSWeightByReadsNoLowSyn, 
                      ActualPerfect=actual_errFree_dnds$dNdSWeightByReadsNoLowSyn, 
                      Expected=expected_dnds$Omega, 
                      Scaling=indelible_dnds$Scaling_factor)

summary(fullDat)
head(fullDat)



fullDatMelt <- reshape2:::melt.data.frame(data=fullDat, na.rm = FALSE, id.vars=c("Site", "Expected", "Scaling"),
                                          measure.vars=c("ActualPerfect", "ActualTypical"),
                                          variable.name="Quality", value.name="Actual")

summary(fullDatMelt)
head(fullDatMelt)


fullDatMelt$Quality <- revalue(fullDatMelt$Quality, c("ActualPerfect"="Without Read Error", "ActualTypical"="With Read Error"))

summary(fullDatMelt)

fig <- ggplot(fullDatMelt  , aes(x=Expected, y=Actual)) + 
  #geom_point(shape=1, alpha=0.5, size=3) +
  geom_point(alpha=0.3, size=2) +
  geom_abline(slope=1, color="Red") +   
  #geom_smooth(method=lm, size=1, se=FALSE) +  
  #geom_smooth(size=1) +  
  #geom_smooth(method="lm", size=1, fill="cadetblue") +  
  #geom_smooth(size=1, fill="cadetblue") +  
  ylab("Inferred Site dN/dS\n") + 
  xlab("\nExpected Site dN/dS") + 
  #coord_fixed(ratio=1) + 
  scale_y_continuous(breaks=0:5, limits=c(0, 5)) + 
  scale_x_continuous(limits=c(0, 3)) + 
  theme_bw(base_size = 12) + 
  theme(panel.grid.major = element_line(size = .3, color = "lightgray")) + 
  theme(strip.text.x = element_text(size=rel(1.3)), axis.title=element_text(size=rel(1)), axis.text=element_text(size=rel(1)),
        panel.margin=unit(c(1,1,1,1),"mm"), plot.margin=unit(c(0,0,0,0),"mm")) + 
  #ggtitle("Inferred Ave dN/dS vs Expected dn/ds") + 
  facet_wrap(~Quality, nrow=1)
print(fig)


# tiff(filename="/home/thuy/gitrepo/SlidingWindowAppNote/scatterplot_typical_v_perfect.tiff",
#      width=4, height=2.6, units="in", compression="none", type="cairo", pointsize=12,
#      res=1200)
svg(filename="/home/thuy/gitrepo/SlidingWindowAppNote/scatterplot_typical_v_perfect.svg",
    width=4, height=2.6, pointsize=1, antialias = "gray")
plot(fig)
dev.off()

cairo_ps(filename="/home/thuy/gitrepo/SlidingWindowAppNote/scatterplot_typical_v_perfect.eps", width=4, height=2.6)
plot(fig)
dev.off()

fullDatNoNA <- fullDat[!is.na(fullDat$ActualTypical) & !is.na(fullDat$ActualPerfect),]
dnds_ccc <- epi.ccc(fullDatNoNA$ActualTypical, fullDatNoNA$Expected)
print(dnds_ccc$rho.c)

dnds_ccc <- epi.ccc(fullDatNoNA$ActualPerfect, fullDatNoNA$Expected)
print(dnds_ccc$rho.c)

dnds_ccc <- epi.ccc(fullDat$ActualTypical, fullDat$Expected)
print(dnds_ccc$rho.c)

dnds_ccc <- epi.ccc(fullDat$ActualPerfect, fullDat$Expected)
print(dnds_ccc$rho.c)

# results <- read.table("./results.csv", header=TRUE, sep=",")
# summary(results)
# head(results)
# 
# lmfit1 <- glm(formula=Condordance~Coverage+Quality+Mutation, data=results)
# summary(lmfit1)
# pvalue <- data.frame(summary(lmfit1)$coefficient)$Pr...t..
# adjpvalue <- p.adjust(pvalue,  method="BH")
# 
# lmfit2 <- glm(formula=Condordance~(Coverage*Mutation) + Quality, data=results)
# summary(lmfit2)
# 
# 
# 
# 
# 
#                      
# results$CoverageFactor <- factor(results$Coverage)
# results$QualityFactor <- paste0(results$Quality, " Read Quality")
# results$CoverQual <- paste0(results$Coverage, results$Quality)
# summary(results)
# head(results)
# ggplot(results, aes(x=Mutation, y=Condordance)) + 
#   geom_line(aes(color=CoverQual, linetype=CoverQual)) +
#   geom_point() + 
#   ylab("Concordance\n") + 
#   xlab("\nSite-Mutation Rate") + 
#   scale_x_continuous(breaks=c(0.004, 0.008, 0.03, 0.06)) + 
#   theme(axis.title=element_text(size=32), axis.text=element_text(size=24) )
# 
#   ggplot(results, aes(x=Mutation, y=Condordance)) + 
#     geom_line(aes(color=CoverageFactor), size=1) +
#     geom_point() + 
#     ylab("Concordance\n") + 
#     xlab("\nSite-Mutation Rate") + 
#     #scale_x_continuous(breaks=c(0.004, 0.008, 0.03, 0.06)) + 
#     scale_x_continuous(breaks=seq(c(0.004, 0.06, 0.02)) + 
#     scale_y_continuous(breaks=seq(0, 1.25, 0.25)) + 
#     #scale_linetype_manual(values = c(rep("solid", 10), rep("dashed", 6))) +
#     #scale_color_manual(breaks=c("F1","F2","F3"),
#     #theme(axis.title=element_text(size=32), axis.text=element_text(size=16) ) + 
#     theme(legend.position="bottom") + 
#     scale_color_discrete(name="Fold Coverage") + 
#     theme(strip.text.x = element_text(size = 14), axis.title=element_text(size=14), axis.text=element_text(size=14)) + 
#     facet_wrap(~QualityFactor)
#   
# 
# 
# 
