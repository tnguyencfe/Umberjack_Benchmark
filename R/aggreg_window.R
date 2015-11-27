#' Looks at causes of window dnds inaccuracy for a single dataset 
#'
#'#+ setup, include=FALSE
library(knitr)
opts_chunk$set(progress=FALSE, verbose=FALSE, warning=FALSE, message=FALSE, width=1200)
options(width=100)

# From data from all windows, aggregates by averaging over windows
library(ggplot2)
library(reshape2)
library(epiR)
library(plyr)
library(stats)
library(RColorBrewer)

NUC_PER_CODON <- 3


# Read locations of input files from local sliding_window_tree_unit_test.config file
CONFIG_FILENAME <- "./aggreg_window.config"
config<-read.table(CONFIG_FILENAME, sep="=", col.names=c("key","value"), as.is=c(1,2))
COLLATE_DNDS_FILENAME <- config[config$key=="COLLATE_DNDS_FILENAME",]$val
EXPECTED_DNDS_FILENAME <-  config[config$key=="EXPECTED_DNDS_FILENAME",]$val
EXPECTED_DNDS_START_NUC_POS <-  as.numeric(config[config$key=="EXPECTED_DNDS_START_NUC_POS",]$val)
EXPECTED_DNDS_END_NUC_POS <-  as.numeric(config[config$key=="EXPECTED_DNDS_END_NUC_POS",]$val)

FULL_POPN_CONSERVE_CSV <- config[config$key=="FULL_POPN_CONSERVE_CSV",]$val
ORIG_CONSERVE_CSV <- config[config$key=="ORIG_CONSERVE_CSV",]$val
ALN_CONSERVE_CSV <- config[config$key=="ALN_CONSERVE_CSV",]$val

#+ results='asis'
kable(config, format="html", caption="config")

#+
# Cols:  Window_Start,Window_End,Reads,CodonSite,CodonDepth,Conserve,Entropy,N,S,EN,ES,dN,dS,dN_minus_dS
collate_dnds <- read.table(COLLATE_DNDS_FILENAME, header=TRUE, na.strings="None", comment.char = "#", sep=",")
collate_dnds$Subst <- collate_dnds$N + collate_dnds$S
colnames(collate_dnds)[grep("^Codons$", colnames(collate_dnds), perl=TRUE)] <- "CodonDepth"
collate_dnds$dNdS <- collate_dnds$dN/collate_dnds$dS
collate_dnds$dNdS[collate_dnds$dS==0] <- NA
collate_dnds$PadPerCodon <- collate_dnds$Pad/(collate_dnds$Reads)
collate_dnds$GapPerCodon <- collate_dnds$Gap/(collate_dnds$Reads)
collate_dnds$AmbigPerCodon <- collate_dnds$Ambig/(collate_dnds$Reads)
collate_dnds$UnknownPerCodon <- (collate_dnds$Pad + collate_dnds$Ambig + collate_dnds$Gap)/collate_dnds$Reads
collate_dnds$ErrPerCodon <- collate_dnds$Err/(collate_dnds$Reads)                                 
collate_dnds$ErrNPerCodon <- collate_dnds$Err_N/collate_dnds$Reads
collate_dnds$ErrSPerCodon <- collate_dnds$Err_S/collate_dnds$Reads 
collate_dnds$TreeDistFractn <- collate_dnds$TreeDist/collate_dnds$Reads
collate_dnds$Is_Break <-   as.factor(collate_dnds$Is_Break)
  
# Average across all codon sites in a window
per_window_ave <- ddply(.data=collate_dnds, .variables=c("Window_Start"), 
                        .fun=function(x) {                            
                          data.frame(Window_Conserve=mean(x$ConserveCodon, na.rm=TRUE),
                                     Window_Entropy=mean(x$EntropyCodon, na.rm=TRUE),
                                     Window_Subst=mean(x$Subst, na.rm=TRUE),
                                     Window_CodonDepth=mean(x$CodonDepth, na.rm=TRUE),
                                     Window_ErrPerCodon=mean(x$ErrPerCodon, na.rm=TRUE),
                                     Window_UnknownPerCodon=mean(x$UnknownPerCodon, na.rm=TRUE),
                                     Window_ErrNPerCodon=mean(x$ErrNPerCodon, na.rm=TRUE),
                                     Window_ErrSPerCodon=mean(x$ErrSPerCodon, na.rm=TRUE)
                          )
                        })
collate_dnds <- merge(x=collate_dnds, y=per_window_ave, by="Window_Start", all=TRUE, sort=TRUE)

dim(collate_dnds)
head(collate_dnds)

#' **Summary Per-Window-Codon Stats**
#' 
summary(collate_dnds)


#' Isolate breakpoints for each window
#' 
num_breaks_per_site <- aggregate(Is_Break ~ CodonSite, data=collate_dnds, FUN=function(x){length(unique(x))})
if (!all(num_breaks_per_site$Is_Break== 1)) {
  stop("There are codon sites that are both breakpoints and not breakpoints at same time.  Impossible!")
}
codon_site_breaks <- aggregate(Is_Break ~ CodonSite, data=collate_dnds, FUN=function(x){unique(x)})
breakpoints <- codon_site_breaks$CodonSite[codon_site_breaks$Is_Break == 1]

#'
#' Read in Expected dN/dS
#' ==============================
#'  

# Cols: Observed S Changes  Observed NS Changes  E[S Sites]  E[NS Sites]	Observed S. Prop.	P{S}	dS	dN	dN-dS	P{S leq. observed}	P{S geq. observed}	Scaled dN-dS	dn/ds
# Parse the expected dnds filename for it start and end nucleotide positions (1-based)
expected_dnds <- read.table(EXPECTED_DNDS_FILENAME, header=TRUE, sep="\t")  # Site  Interval	Scaling_factor	Rate_class	Omega
expected_dnds$CodonSite <- as.numeric(expected_dnds$Site)  + 1  # Site is 0-based codon site, CodonSite is 1-based codon site
expected_dnds$Omega <- expected_dnds$dN/expected_dnds$dS
expected_dnds$Omega[expected_dnds$dS == 0] <- NA
colnames(expected_dnds)[grep("Observed.S.Changes", colnames(expected_dnds))] <- "S"
colnames(expected_dnds)[grep("Observed.NS.Changes", colnames(expected_dnds))] <- "N"
expected_dnds$Subst <- expected_dnds$S + expected_dnds$N


# check consistency.  Because the windows slide right, it's possible that the end of the genome might not be covered by any windows
# depending on window size and slide.
if (min(collate_dnds$CodonSite, na.rm=TRUE) != min(expected_dnds$CodonSite, na.rm=TRUE)) {
  stop("Expected dN/dS Sites don't line up with Actual dN/dS Sites")
}

#' **Summary Expected dN/dS Stats**
#' 
summary(expected_dnds)



#' **How is our coverage?**
#' -------------------------------------

# Plots the window stats for each window at each codon site
plot_window_nums <- function(collate_dnds, colname, descname) {
  fig <- ggplot(collate_dnds, aes(x=CodonSite)) + 
    geom_point( aes_string(y=colname), alpha=0.5, shape=1) +  
    geom_smooth(aes_string(y=colname)) + 
    ylab(paste0(descname, "\n")) + 
    xlab("\nCodon Site") + 
    #theme(axis.title=element_text(size=32), axis.text=element_text(size=24) ) + 
    ggtitle(paste0(descname, " From Each Window"))
  
  if (length(breakpoints) > 0) {
    fig <- fig + geom_vline(xintercept=breakpoints, color="red")
  }
    
  print(fig)
}


plot_window_nums(collate_dnds, "Reads", "Max Read Depth")
plot_window_nums(collate_dnds, "CodonDepth", "Unambiguous Codon Count")
plot_window_nums(collate_dnds, "ConserveCodon", "Conservation with Consensus Codon")
plot_window_nums(collate_dnds, "EntropyCodon", "Codon Entropy")
plot_window_nums(collate_dnds, "Subst", "Substitutions")
plot_window_nums(collate_dnds, "N", "Nonsynonymous Substitutions")
plot_window_nums(collate_dnds, "S", "Synonymous Substitutions")
plot_window_nums(collate_dnds, "ErrPerCodon", "Errors per Codon")
plot_window_nums(collate_dnds, "UnknownPerCodon", "Unknown Bases Per Codon")
plot_window_nums(collate_dnds, "PadPerCodon", "Left/Right Pad Per Codon")
plot_window_nums(collate_dnds, "AmbigPerCodon", "N Per Codon")
plot_window_nums(collate_dnds, "GapPerCodon", "Deletion Per Codon")
plot_window_nums(collate_dnds, "TreeDist", "Diff From True Tree")
plot_window_nums(collate_dnds, "TreeDistFractn", "Diff From True Tree Normalized by Tips")



#' How much variation is there from window to window for the same codon site?
#' -------------------------------------------------------------------------------
#' 
#' The coefficient of variation (CV) is defined as the ratio of the standard deviation to the mean
#' 

# Plots the coefficient of variance across windows at each codon site
# Assumes that collate_dnds has data
plot_coeff_var <- function(collate_dnds, colname, descname) {
  coeff_var <- ddply(.data=collate_dnds, 
                     .variables="CodonSite",                          
                     .fun=function(x) {
                       ave <- mean(x[, colname], na.rm=TRUE)
                       stddev <- sd(x[, colname], na.rm=TRUE)
                       cv <- stddev/ave
                       return (data.frame(#ave=ave, stddev=stddev, 
                         cv=cv))
                     })
  
  coeff_var <- coeff_var[order(coeff_var$CodonSite), ]
  
  # Check that we have the correct dimensions
  if(nrow(coeff_var) != length(unique(collate_dnds$CodonSite))) {
    stop("Invalid dimensions.  Number of codon sites for coeff_var != number of codon sites for collate_dnds")
  }
  
#   print(summary(coeff_var))
  
  fig <- ggplot(coeff_var, aes(x=CodonSite, y=cv)) + 
    geom_line(shape=1, alpha=0.5) +  
    geom_smooth() + 
    ylab(paste0("Coeff of Var of ", descname, "\n")) + 
    xlab("\nCodon Site") + 
    #theme(axis.title=element_text(size=32), axis.text=element_text(size=24) ) + 
    ggtitle(paste0("Coefficient of Variance Across Windows \n for ", descname))

  if (length(breakpoints) > 0) {
    fig <- fig + geom_vline(xintercept=breakpoints, color="red")
  }
  print(fig)
}

plot_coeff_var(collate_dnds, "ConserveCodon", "Codon Conservation")
plot_coeff_var(collate_dnds, "Reads", "Max Reads")
plot_coeff_var(collate_dnds, "EntropyCodon", "Codon Entropy")
plot_coeff_var(collate_dnds, "N", "Nonsynonymous Substitutions")
plot_coeff_var(collate_dnds, "S", "Synonymous Substitutions")
plot_coeff_var(collate_dnds, "Subst", "Substitutions")
plot_coeff_var(collate_dnds, "dNdS", "dN/dS")
plot_coeff_var(collate_dnds, "dN", "dN")
plot_coeff_var(collate_dnds, "dS", "dS")
plot_coeff_var(collate_dnds, "dN_minus_dS", "Scaled dN-dS")
plot_coeff_var(collate_dnds, "TreeDistFractn", "Diff From True Tree Normalized by Tips")


#' What is the Relationship Between Sequence Diversity And Full Population Phylogeny Substitution Count?
#' ------------------------------------------
#' 

# Takes a *.conserve.csv file with nucleotide stats for unsliced sequences
# Expects cols: "NucSite", "Conserve", "Entropy", "NucDepth", "CodonDepth"  
# Returns a dataframe with elements:
# - values average per-codon values for Conserve, ConserveCodon, Entropy, CodonDepth
fill_nowindow_codon_conserve_dat <- function(conserve_csv) {
  # Per-nucleotide conservation, ConserveCodonnce, Entropy
  conserve <- read.table(conserve_csv, header=TRUE, sep=",")
  return (conserve)
}



# Scatterplots Phylogeny Substitution Count Vs Sequence Diversity
# INPUT:
# - codon_dat:  dataframe with cols: CodonSite,CodonDepth,Entropy,N,S,Subst,ConserveCodon
# - diversity_title:  x-axes title prefix
# - subst_title:  y-axes title prefix
# - diversity_measure_suffix:  diversity measure suffix
# - subst_measure_suffix:  subsitution measure suffix
# - codon_diversity:  dataframe with cols: CodonSite,CodonDepth,Entropy,ConserveCodon.  Only used if codon_dat is NULL
# - codon_diversity:  dataframe with cols: CodonSite,N,S, Subst.  Only used if codon_dat is NULL
plot_subst_vs_diversity <- function(codon_dat, diversity_title, subst_title, diversity_measure_suffix, subst_measure_suffix, 
                                    codon_diversity=NULL, codon_subst=NULL) {
  if (is.null(codon_dat) & (is.null(codon_diversity)  | is.null(codon_subst))) {
    stop("Either define codon_dat or both codon_diversity and codon_subst")
  } else if (is.null(codon_dat))  {
    codon_dat <- merge(x=subset(codon_subst, select=c(CodonSite, S, N, Subst)),
                       y=subset(codon_diversity, select=c(CodonSite, CodonDepth, ConserveCodon, EntropyCodon)),
                       by="CodonSite", all=TRUE)
  }
  phylo_subst_melt <- reshape2:::melt(data=subset(codon_dat, select=c(CodonSite, S, N, Subst)),
                                      id.vars="CodonSite",
                                      measure.vars=c("S", "N", "Subst"),
                                      variable.name="PhyloMeasure", 
                                      value.name="PhyloVal")
  phylo_subst_melt$PhyloMeasure <- paste0(phylo_subst_melt$PhyloMeasure, subst_measure_suffix)
  
  leaf_diversity_melt <- reshape2:::melt(data=subset(codon_dat, select=c(CodonSite, CodonDepth, ConserveCodon, EntropyCodon)),
                                         id.vars="CodonSite",
                                         measure.vars=c("CodonDepth", "ConserveCodon", "EntropyCodon"),
                                         variable.name="LeafMeasure", 
                                         value.name="LeafVal")
  leaf_diversity_melt$LeafMeasure <- paste0(leaf_diversity_melt$LeafMeasure, diversity_measure_suffix)
  
  combo <- merge(x=phylo_subst_melt, y=leaf_diversity_melt, by="CodonSite", all=TRUE)
  
  fig <- ggplot(combo, aes(x=LeafVal, y=PhyloVal)) + 
    facet_grid(PhyloMeasure~LeafMeasure, scales="free_x") + 
    geom_point(shape=1, alpha=0.1) + 
    geom_smooth(method="lm") + 
    xlab(paste0("\n", diversity_title, " Sequence Diversity")) + 
    ylab(paste0(subst_title, " SLAC Count\n")) + 
    ggtitle(paste0("Scatterplot ", subst_title, " SLAC Count Vs ", diversity_title, " Sequence Diversity"))
  print(fig)
}

#'
#' Window Sequence Diversity Vs Window Phylogeny Substitution Count
#' =====================================================================
#'  

#+ fig.width=12, fig.height=12
plot_subst_vs_diversity(collate_dnds, diversity_title="Window", subst_title="Window", 
                        diversity_measure_suffix=".Win", subst_measure_suffix=".Win")

#' Full Population Sequence Diversity Vs Full Population Phylogeny Substitution Count
#' =====================================================================
#'  
full_popn_conserve_dat <- fill_nowindow_codon_conserve_dat(FULL_POPN_CONSERVE_CSV)
#+ results='asis'
kable(summary(full_popn_conserve_dat), format="html", caption="Full Population Codon Stats")

#+ fig.width=12, fig.height=12
plot_subst_vs_diversity(codon_dat=NULL, diversity_title="Full Population", subst_title="Full Population", 
                        diversity_measure_suffix=".Full", subst_measure_suffix=".Full", 
                        codon_diversity=full_popn_conserve_dat, codon_subst=expected_dnds)


#'
#' Original Read Diversity Vs Full Population Phylogeny Substitution Count
#' =====================================================================
#'  

orig_conserve_dat <- fill_nowindow_codon_conserve_dat(ORIG_CONSERVE_CSV)

#+ results='asis'
kable(summary(orig_conserve_dat), format="html", caption="Original Reads Codon Stats")

#+ fig.width=12, fig.height=12
plot_subst_vs_diversity(codon_dat=NULL, diversity_title="Original Read", subst_title="Full Population", 
                        diversity_measure_suffix=".OrigRead", subst_measure_suffix=".Full", 
                        codon_diversity=orig_conserve_dat, codon_subst=expected_dnds)

#'
#' Aligned Read Diversity Vs Full Population Phylogeny Substitution Count
#' =====================================================================
#'  
aln_conserve_dat <- fill_nowindow_codon_conserve_dat(ALN_CONSERVE_CSV)

#+ results='asis'
kable(summary(aln_conserve_dat), format="html", caption="Aligned Reads Codon Stats")

#+ fig.width=12, fig.height=12
plot_subst_vs_diversity(codon_dat=NULL, diversity_title="Aligned Read", subst_title="Full Population", 
                        diversity_measure_suffix="Aligned", subst_measure_suffix=".Full", 
                        codon_diversity=aln_conserve_dat, codon_subst=expected_dnds)

#'
#' What is the Effect of Window-ing?
#' ------------------------------------------
#' 
# Scatterplots Window Sequence Diversity Vs Non Window Sequence Diversity
# Scaterplots Window SLAC Counts Vs Non Window Sequence Diversity
# I.e.  Compares effect of Window-ing on Diversity.  Compares effect of Window-ing on Substitutions.
# INPUT:
# - nowindow_leaf_diversity:  dataframe with cols:  CodonSite, CodonDepth, ConserveCodon, Entropy.  
#     Values are calculated per-nucleotide and averaged per-codon.  No windows used in calculations.
# - nowindow_phylo_subst:  dataframe with cols: N, S, Subst
#     Values are calculated per codon against full population.  No windows used in calculations.  No reads used in calculations.
# - win_codon_dat:  dataframe with cols: CodonSite,CodonDepth,Entropy,N,S,Subst,ConserveCodon 
# - nowindow_title:  nowindow axes title
# - window_title:  window_title axes title
# - nowindow_measure_suffix:  nowindow measure suffix
# - window_measure_suffix:  window measure suffix
plot_win_subst_div_vs_nowin_div <- function(nowindow_leaf_diversity, nowindow_phylo_subst, win_codon_dat,
                                            nowindow_title, window_title, nowindow_measure_suffix, window_measure_suffix) {
  
  # Just compare full population entropy vs window SLAC counts
  # No need to compare both full population entropy and and conservation against window SLAC counts once
  # we verify that full population conservation and entropy are legit
  win_subst_nowindow_leaf_combo <- merge(x=win_codon_dat, y=nowindow_leaf_diversity, by="CodonSite", all=TRUE, 
                                         suffixes=c(window_measure_suffix, nowindow_measure_suffix))
  
  
  fig <- ggplot(win_subst_nowindow_leaf_combo, aes_string(x=paste0("EntropyCodon", nowindow_measure_suffix), 
                                                            y=paste0("EntropyCodon", window_measure_suffix))) + 
    geom_point(shape=1, alpha=0.1) + 
    geom_smooth(method="lm") + 
    xlab("Full Popn Codon Entropy") + 
    ylab ("Window Codon Entropy") + 
    ggtitle(paste0("Scatterplot Window Entropy vs Full Popn Entropy"))
  print(fig)
  
  fig <- ggplot(win_subst_nowindow_leaf_combo, aes_string(x=paste0("EntropyCodon", nowindow_measure_suffix), y="N")) +     
    geom_point(shape=1, alpha=0.1) + 
    geom_smooth(method="lm") + 
    xlab("Full Popn Codon Entropy") + 
    ylab ("Window N") + 
    ggtitle(paste0("Scatterplot Window N vs Full Popn Entropy"))
  print(fig)
  
  fig <- ggplot(win_subst_nowindow_leaf_combo, aes_string(x=paste0("EntropyCodon", nowindow_measure_suffix), y="S")) +     
    geom_point(shape=1, alpha=0.1) + 
    geom_smooth(method="lm") + 
    xlab("Full Popn Codon Entropy") + 
    ylab ("Window S") + 
    ggtitle(paste0("Scatterplot Window S vs Full Popn Entropy"))
  print(fig)
  
  fig <- ggplot(win_subst_nowindow_leaf_combo, aes_string(x=paste0("EntropyCodon", nowindow_measure_suffix), y="Subst")) +     
    geom_point(shape=1, alpha=0.1) + 
    geom_smooth(method="lm") + 
    xlab("Full Popn Codon Entropy") + 
    ylab ("Window Subs") + 
    ggtitle(paste0("Scatterplot Window Subs vs Full Popn Entropy"))
  print(fig)
  
}

# Scatterplots Window SLAC Counts vs Full Population SLAC Counts
plot_win_subst_vs_full_subst <- function(nowindow_phylo_subst, win_codon_dat,
                                         nowindow_title, window_title, nowindow_measure_suffix, window_measure_suffix) {
  
  # Density Plot of Window-Codon Sites at Each Full Population N, S, Subst Count
  # Scatterplot Window SLAC Vs  Full Population SLAC
  win_nowindow_subst_combo <- merge(x=subset(nowindow_phylo_subst, select=c(CodonSite, S, N, Subst)),
                                    y=subset(win_codon_dat, select=c(CodonSite, S, N, Subst)), 
                                    by="CodonSite", all=TRUE, suffixes=c(nowindow_measure_suffix, window_measure_suffix))
  
  sapply(c("N", "S", "Subst"), function(col) {
    
#     # Cumulative Frequency plot Window-Codon Sites at Each Full POpulation N, S, Subst Count
#     fig <- ggplot(win_nowindow_subst_combo, aes_string(x=paste0(col, nowindow_measure_suffix))) +
#       stat_ecdf() + 
#       xlab(paste0("\n", nowindow_title, " " , col)) + 
#       ylab("Cumulative Frequency Window-Codon Sites\n") + 
#       ggtitle(paste0("Cumulative Frequency Plot Window-Codon Sites by  ", nowindow_title, " ", col))
#     print(fig)
    
    
    # Scatterplot Window N, S, Subst Vs Full Population N, S, Subst
    fig <- ggplot(win_nowindow_subst_combo, aes_string(x=paste0(col, nowindow_measure_suffix), 
                                                       y=paste0(col, window_measure_suffix))) + 
      geom_point(shape=1, alpha=0.1) + 
      geom_smooth(method="lm") + 
      geom_abline(slope=1, color="red") + 
      xlab(paste0("\n", nowindow_title, " " , col)) + 
      ylab(paste0(window_title, " ", col, "\n")) + 
      ggtitle(paste0("Scatterplot ", window_title, " ", col, " Vs ", nowindow_title, " ", col))
    print(fig)
  })
  
}

#'
#' Effect of Windowing on Phylogeny Compared to Full Population
#' =================================================
#' 

#+ fig.width=12, fig.height=12
plot_win_subst_vs_full_subst(nowindow_phylo_subst=expected_dnds, win_codon_dat=collate_dnds,
                             nowindow_title="Full Population", window_title="Window", 
                             nowindow_measure_suffix=".Full", window_measure_suffix=".Win")

#'
#' Effect on Windowing in Sequencing Diversity and Phylogeny Substitutions Compared to Full Population
#' =================================================
#' 

#+ fig.width=12, fig.height=12
plot_win_subst_div_vs_nowin_div(nowindow_leaf_diversity=full_popn_conserve_dat, 
                                nowindow_phylo_subst=expected_dnds, win_codon_dat=collate_dnds,
                                nowindow_title="Full Population", window_title="Window", 
                                nowindow_measure_suffix=".Full", window_measure_suffix=".Win")

#'
#' Effect on Windowing in Sequencing Diversity and Phylogeny Substitutions Compared to Original ART reads
#' =================================================
#' 

#+ fig.width=12, fig.height=12
plot_win_subst_div_vs_nowin_div(nowindow_leaf_diversity=orig_conserve_dat, 
                                nowindow_phylo_subst=expected_dnds, win_codon_dat=collate_dnds,
                                nowindow_title="Original Reads", window_title="Window", 
                                nowindow_measure_suffix=".Orig", window_measure_suffix=".Win")


# #' Effect on Windowing in Sequencing Diversity and Phylogeny Substitutions Compared to Aligned ART reads
# #' =================================================
# #' 
# 
# #+ results='asis'
# kable(summary(aln_conserve_dat), format="html", caption="Aligned Reads Codon Stats")
# 
# #+ fig.width=12, fig.height=12
# plot_win_subst_div_vs_nowin_div(nowindow_leaf_diversity=aln_conserve_dat, 
#                                 nowindow_phylo_subst=expected_dnds, win_codon_dat=collate_dnds,
#                                 nowindow_title="Aligned Reads", window_title="Window", 
#                                 nowindow_measure_suffix=".Aln", window_measure_suffix=".Win")
# 
# 



#'
#' What is the Effect of Window-ing on dN-dS?
#' ------------------------------------------
#' 

#' **Scatterplot actual vs expected dn ds together**
fullDat <- merge(x=subset(collate_dnds, select=c(CodonSite, Window_Start, dN, dS, dN_minus_dS, 
                                                 TreeDist, TreeDistFractn, TreeDepth, TreeLen, ErrPerCodon, UnknownPerCodon,
                                                 N, S)), 
                 y=subset(expected_dnds, select=c(CodonSite, dN, dS, Omega, Scaled.dN.dS)),
                 by=c("CodonSite"),
                 suffixes=c(".Act", ".Exp"))
colnames(fullDat)[grep("Omega", colnames(fullDat))] <- "dNdS.Exp"
colnames(fullDat)[grep("dN_minus_dS", colnames(fullDat))] <- "dN_minus_dS.Act"
colnames(fullDat)[grep("Scaled.dN.dS", colnames(fullDat))] <- "dN_minus_dS.Exp"

fullDat$dNdS.Act <- fullDat$dN.Act/fullDat$dS.Act
fullDat$dNdS.Act[fullDat$S == 0] <- NA

fullDat$IsAmbigSub <- as.factor((fullDat$N > 0 & fullDat$N < 1) | (fullDat$S > 0 & fullDat$S < 1))

summary(fullDat)

#'Scatterplot Expected vs Inferred Averaged Across Windows
#'  
#+ fig.width=10, fig.height=10
fig <- ggplot(fullDat , aes(x=dN_minus_dS.Exp, y=dN_minus_dS.Act, color=IsAmbigSub)) + 
  geom_smooth(method=lm, size=4, color="#A30052", fill="#FF99CC") +
  geom_point(alpha=0.5) +
  geom_abline(slope=1) + 
  ylab("Inferred\n") + 
  xlab("\nExpected") + 
  ggtitle("dN_minus_dS.Act Vs dN_minus_dS.Exp")
print(fig)

fig <- ggplot(fullDat , aes(x=dNdS.Exp, y=dNdS.Act, color=IsAmbigSub)) + 
  geom_smooth(method=lm, size=4, color="#A30052", fill="#FF99CC") +
  geom_point(alpha=0.5) +
  geom_abline(slope=1) + 
  scale_y_log10() + 
  ylab("Inferred\n") + 
  xlab("\nExpected") + 
  ggtitle("dNdS.Act Vs dNdS.Exp")
print(fig)



#' What is the dn-dS along Genome?  Check off-by-one errors
#' 
#+ fig.width=20, fig.height=30
fig <- ggplot(fullDat) + 
  geom_point(aes(x=CodonSite, y=abs(dN_minus_dS.Act-dN_minus_dS.Exp)), color="black") + 
  geom_line(aes(x=CodonSite, y=abs(dN_minus_dS.Act-dN_minus_dS.Exp)), color="black") + 
  geom_vline(data=as.data.frame(breakpoints), xintercept=breakpoints, color="blue") + 
  
  facet_wrap(~Window_Start, ncol=1)
print(fig)

# #' What is the dn/dS along Genome?  Check of off-by-one
# #' 
# #+ fig.width=20, fig.height=20
# fig <- ggplot(fullDat) + 
#   geom_point(aes(x=CodonSite, y=dNdS.Act), color="black") + 
#   geom_point(aes(x=CodonSite, y=dNdS.Exp), color="red") + 
#   geom_line(aes(x=CodonSite, y=dNdS.Act), color="black") + 
#   geom_line(aes(x=CodonSite, y=dNdS.Exp), color="red") + 
#   geom_vline(data=as.data.frame(CODON_BREAKPOINTS), xintercept=CODON_BREAKPOINTS, color="blue") + 
#   facet_wrap(~Window_Start, ncol=1)
# print(fig)


#'**Concordance**
#'

# Returns a table of Lin's concordance correlation values
corr_vals <- data.frame(dNdS = epi.ccc(fullDat$dNdS.Act, fullDat$dNdS.Exp)$rho.c$est,
                        dNdS_NoAmbig = epi.ccc(fullDat$dNdS.Act[fullDat$IsAmbigSub == FALSE], 
                                               fullDat$dNdS.Exp[fullDat$IsAmbigSub == FALSE])$rho.c$est,
                        dN_minus_dS = epi.ccc(fullDat$dN_minus_dS.Act, fullDat$dN_minus_dS.Exp)$rho.c$est,
                        dN_minus_dS_NoAmbig = epi.ccc(fullDat$dN_minus_dS.Act[fullDat$IsAmbigSub == FALSE], 
                                                      fullDat$dN_minus_dS.Exp[fullDat$IsAmbigSub == FALSE])$rho.c$est
)

kable(corr_vals, format="html", caption="Concordance Correlation")


# Plots Lin's concordance correlation values by Non-Window Sequence Conservation
# Plot count of dn/ds values by Sequence Conservation
plot_table_corr_by_div <- function(nowindow_leaf_diversity) {
  
  nowindow_leaf_diversity$ConserveCodonBin <- cut(nowindow_leaf_diversity$ConserveCodon, 
                                            breaks=seq(0, max(nowindow_leaf_diversity$ConserveCodon, na.rm=TRUE)+0.1, 0.1),
                                            include.lowest = TRUE, right=FALSE)
  full_div_combo <- merge(x=fullDat, 
                          y=subset(nowindow_leaf_diversity, select=c(CodonSite, ConserveCodon, ConserveCodonBin)), 
                          by="CodonSite", all=TRUE, suffixes=c(".Win", ".Full"))
  
  corr_vals <-ddply(.data=full_div_combo,
                    .variables="ConserveCodonBin",
                    .fun=function(x) {
                      data.frame(                                               
                        dNdS = epi.ccc(x$dNdS.Act, x$dNdS.Exp)$rho.c$est,
                        dNdS_NoAmbig = epi.ccc(x$dNdS.Act[x$IsAmbigSub == FALSE], 
                                               x$dNdS.Exp[x$IsAmbigSub == FALSE])$rho.c$est,
                        dN_minus_dS = epi.ccc(x$dN_minus_dS.Act, x$dN_minus_dS.Exp)$rho.c$est,
                        dN_minus_dS_NoAmbig = epi.ccc(x$dN_minus_dS.Act[x$IsAmbigSub == FALSE], 
                                                      x$dN_minus_dS.Exp[x$IsAmbigSub == FALSE])$rho.c$est
                      )
                    })
  
  corr_vals_melt <- reshape2:::melt(data=corr_vals, id.vars="ConserveCodonBin", variable.name="dnds_calc", value.name="Concordance")
  #corr_vals_melt$Concordance[is.na(corr_vals_melt$Concordance)] <-  0
  fig <- ggplot(corr_vals_melt, aes(x=ConserveCodonBin, y=Concordance)) + 
    geom_point() + 
    geom_line(aes(x=ConserveCodonBin, y=Concordance, group=dnds_calc)) + 
    theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
    facet_wrap(~dnds_calc)
  print(fig)
  
}

#' What is the Concordance By Full Population Diversity
#' =======================================================================================
#' 
#+ fig.width=15, fig.height=12
plot_table_corr_by_div(full_popn_conserve_dat)



# #'
# #' What is the Concordance By original Read Diversity
# #' =======================================================================================
# #' 
# #+ fig.width=15, fig.height=12
# plot_table_corr_by_div(orig_conserve_dat)

#'
#' What is the Concordance By Tree Accuracy (weighted robinson foulds)
#' =======================================================================================
#' 
#+ fig.width=15, fig.height=12

fullDat$TreeDistFractnBin <- cut(fullDat$TreeDistFractn, 
                                          breaks=seq(0.013, 0.023, 0.001),
                                          include.lowest = TRUE, right=FALSE)
corr_vals <-ddply(.data=fullDat,
                  .variables="TreeDistFractnBin",
                  .fun=function(x) {
                    data.frame(                                               
                      dNdS = epi.ccc(x$dNdS.Act, x$dNdS.Exp)$rho.c$est,
                      dNdS_NoAmbig = epi.ccc(x$dNdS.Act[x$IsAmbigSub == FALSE], 
                                             x$dNdS.Exp[x$IsAmbigSub == FALSE])$rho.c$est,
                      dN_minus_dS = epi.ccc(x$dN_minus_dS.Act, x$dN_minus_dS.Exp)$rho.c$est,
                      dN_minus_dS_NoAmbig = epi.ccc(x$dN_minus_dS.Act[x$IsAmbigSub == FALSE], 
                                                    x$dN_minus_dS.Exp[x$IsAmbigSub == FALSE])$rho.c$est
                    )
                  })

corr_vals_melt <- reshape2:::melt(data=corr_vals, id.vars="TreeDistFractnBin", variable.name="dnds_calc", value.name="Concordance")

fig <- ggplot(corr_vals_melt, aes(x=TreeDistFractnBin, y=Concordance)) + 
  geom_point() + 
  geom_line(aes(x=TreeDistFractnBin, y=Concordance, group=dnds_calc)) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  facet_wrap(~dnds_calc)
print(fig)

