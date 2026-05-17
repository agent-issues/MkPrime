.libPaths(c("/nobackup/pjjg18/mkp-sim3-multirep-v3/lib", .libPaths()))
library(MkPrime)
s <- as.data.frame(ReadMkLog("/nobackup/pjjg18/mkp-sim3-phi2/results/aware.log"))
cat("pi0:  min=", round(min(s$pi0),3), " med=", round(median(s$pi0),3),
    " max=", round(max(s$pi0),3), " mean=", round(mean(s$pi0),3), "\n")
cat("phi:  min=", round(min(s$phi),3), " med=", round(median(s$phi),3),
    " max=", round(max(s$phi),3), "\n")
cat("TL:   min=", round(min(s$tree_length),2), " med=", round(median(s$tree_length),2),
    " mean=", round(mean(s$tree_length),2), "\n")
cat("nSamples:", nrow(s), "  truth pi0=0.75  truth phi=2\n")
