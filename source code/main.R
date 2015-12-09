# 
require('clue')
# KCCA lib.
require('kernlab')
# Plot.
require('ggplot2')
# Set the path of environment.
setwd('/Users/Salmon/Dropbox/NTUST/1041 - 隨班附讀/Data mining 戴碧如老師/Projects/source code/')
# Import other source files.
source('src/sliptAndMissDataset.R')
source('src/CoKL.R')
source('src/NMI.R')

# Read the 'seeds' data from UCI datasets.
seeds           <- read.table('datasets/seeds_dataset.txt', sep='\t', header=FALSE)
#seeds           <- read.table('datasets/seeds_less_dataset.txt', sep='\t', header=FALSE)
colnames(seeds) <- c('area', 'perimeter', 'compactness', 'length_of_kernel', 'width_of_kernel', 'asymmetry_coefficient', 'length_of_kernel_groove', 'class_of_wheat')
# Step 0.
missingRate       <- seq(0.1, 0.9, by=0.025)
# missingRate       <- 0.5
seedSlipts        <- lapply(missingRate, function (rate) { return (sliptAndMissDataset(seeds[-ncol(seeds)], rate)) })
names(seedSlipts) <- paste("missingRate_", missingRate, sep="")
# Seeds get kernel via CoKL algorithm.
CoKL_Kernels  <- lapply(seedSlipts, function (slipt) { CoKL(slipt) })
# CoKL + KCCA and NMI.
NMI_results  <- lapply(CoKL_Kernels, function (kernel) {
  mr <- kernel$missingRate
  or <- kernel$order
  # KCCA_res   <- kcca(x=kernel$X$all, y=kernel$Y$all, kernel="vanilladot", kpar=list(), ncomps=7)
  KCCA_res   <- kcca(x=kernel$X$all, y=kernel$Y$all, kernel="rbfdot", kpar=list(sigma=0.1), gamma=0.1, ncomps=7)
  kmeans_res <- kmeans(cbind(KCCA_res@xcoef, KCCA_res@ycoef), 3)
  nmi_res    <- NMI(kmeans_res$cluster[or], seeds$class_of_wheat)
  return (nmi_res)
})
# 
missingRate <- c()
nmiValue <- c()
for (mr in names(NMI_results)) {
  missingRate <- c(missingRate, as.numeric(strsplit(mr, split="missingRate_")[[1]][2]))
  nmiValue    <- c(nmiValue, as.numeric(NMI_results[mr]))
}
# Final, plot.
qplot(x=missingRate, y=nmiValue, xlim=c(0.1, 0.9), ylim=c(0, 1), xlab="Missing Rate", ylab="NMI", colour="red") + geom_line(colour="red") + theme(legend.position="none")
