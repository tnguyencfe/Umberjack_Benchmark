library(ggplot)

frag <- read.table('/home/thuy/gitrepo/Umberjack_Benchmark/simulations/data/benchmark_art_profile/seq_error_rate.fraglen.csv', sep=",", header=TRUE)
head(frag)
summary(frag)


#_indel0.00045_qs0_cover9_fragmean136_fragstd94_seed6298981545792848188.sam
samsplit <- strsplit(as.character(frag$Sam), "_")
frag$indel <- unlist(lapply(samsplit, 
                function(x) {
                  indelstr <- x[grep("indel", x)]
                  indelrate <- as.numeric(sub("indel", "", indelstr))
                  }))
frag$qualshift <- unlist(lapply(samsplit, 
                            function(x) {
                              qualshiftstr <- x[grep("qs", x)]
                              rate <- as.numeric(sub("qs", "", qualshiftstr))
                            }))
frag$cover <- unlist(lapply(samsplit, 
                                function(x) {
                                  coverstr <- x[grep("cover", x)]
                                  rate <- as.numeric(sub("cover", "", coverstr))
                                }))
frag$fragmean <- unlist(lapply(samsplit, 
                                function(x) {
                                  fragmeanstr <- x[grep("fragmean", x)]
                                  rate <- as.numeric(sub("fragmean", "", fragmeanstr))
                                }))
frag$fragstd <- unlist(lapply(samsplit, 
                               function(x) {
                                 fragstdstr <- x[grep("fragstd", x)]
                                 rate <- as.numeric(sub("fragstd", "", fragstdstr))
                               }))
shortfrag <- frag[frag$fragmean < 200, ]
head(shortfrag)
summary(shortfrag)
shortfrag$cover <- as.factor(shortfrag$cover)

fig <- ggplot(shortfrag, aes(x=Fraglen, weight=Count)) + 
  #geom_density(color="black", fill="blue") + 
  geom_histogram(color="black", fill="blue", binwidth=50) + 
  facet_wrap(~cover)
print(fig)
