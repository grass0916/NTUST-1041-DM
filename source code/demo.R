# 
require('clue')
# KCCA lib.
require('kernlab')
# Set the path of environment.
setwd('/Users/Salmon/Dropbox/NTUST/1041 - 隨班附讀/Data mining 戴碧如老師/Projects/source code/')
# Import other source files.
source('src/sliptAndMissDataset.R')
source('src/CoKL.R')

testM <- matrix(c(
  c(15.26, 14.84, 0.871, 5.763, 3.312, 2.221, 5.22),
  c(17.63, 15.98, 0.8673, 6.191, 3.561, 4.076, 6.06),
  c(16.84, 15.67, 0.8623, 5.998, 3.484, 4.675, 5.877),
  c(13.07, 13.92, 0.848, 5.472, 2.994, 5.304, 5.395),
  c(14.88, 14.57, 0.8811, 5.554, NA, NA, NA),
  c(13.32, 13.94, 0.8613, 5.541, NA, NA, NA),
  c(NA, NA, NA, NA, 3.337, 2.699, 4.825)),
  , byrow=TRUE, ncol=7)

testSet <- list()
testSet$missingRate_0.5 <- list()
testSet$missingRate_0.5$X  <- data.frame(testM[,1:4])
testSet$missingRate_0.5$Y  <- data.frame(testM[,5:7])
testSet$missingRate_0.5$c  <- 4
testSet$missingRate_0.5$m1 <- 2
testSet$missingRate_0.5$m2 <- 1

CoKL_Kernels <- lapply(testSet, function (slipt) { CoKL(slipt) })
KCCA_results <- kcca(CoKL_Kernels$missingRate_0.5$X$all, CoKL_Kernels$missingRate_0.5$Y$all)
kmeans(cbind(KCCA_results@xcoef, KCCA_results@ycoef), 3)

print(cbind(KCCA_results@xcoef, KCCA_results@ycoef))