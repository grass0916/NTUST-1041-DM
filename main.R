# Set the path of environment.
setwd('/Projects/source code/datasets/')
# CoKL algorithm.
CoKL <- function (slipt) {
  c  <- slipt$c
  m1 <- slipt$m1
  m2 <- slipt$m2
  # --- Inner fuctions --- #
  # Step 1: Fix the missing values by means of each column.
  fixMissingValues <- function (s) {
    # s: set.
    s$X[is.na(s$X[,1]),] <- as.list(apply(s$X, 2, function (column) { mean(na.omit(column)) }))
    s$Y[is.na(s$Y[,1]),] <- as.list(apply(s$Y, 2, function (column) { mean(na.omit(column)) }))
    return (s)
  }
  # Step 2: Calculate initial kernel matrix Kx = BB' and Ky = AA'.
  kernelMatrix <- function (s) {
    # s: set.
    kernelX <- as.matrix(s$Y) %*% t(as.matrix(s$Y))
    kernelY <- as.matrix(s$X) %*% t(as.matrix(s$X))
    print(kernelY)
    print(chol(kernelY))
    kernel  <- list(X=matrixSeperator(kernelX), Y=matrixSeperator(kernelY))
    return (kernel)
  }
  # Step 3 & 4: Calculate Laplacian Lx and Ly.
  graphLaplacian <- function (k) {
    # k: kernel.
    diagonalX  <- diag(apply(k$X$all, 1, function (row) { sum(row) }), nrow(k$X$all), ncol(k$X$all))
    diagonalY  <- diag(apply(k$Y$all, 1, function (row) { sum(row) }), nrow(k$Y$all), ncol(k$Y$all))
    laplacianX <- diagonalX - k$X$all
    laplacianY <- diagonalY - k$Y$all
    laplacian  <- list(X=matrixSeperator(laplacianX), Y=matrixSeperator(laplacianY))
    return (laplacian)
  }
  matrixSeperator <- function (m) {
    # s: set ; m: original matrix.
    s <- list()
    s$all  <- matrix(m, ncol=c+m1+m2, byrow=TRUE)
    s$cc   <- matrix(m[1:c, 1:c], ncol=c, byrow=TRUE)
    s$cm1  <- if (m1 == 0) matrix() else matrix(m[1:c, (c+1):(c+m1)], ncol=m1, byrow=TRUE)
    s$cm2  <- if (m2 == 0) matrix() else matrix(m[1:c, (c+1+m1):(c+m1+m2)], ncol=m2, byrow=TRUE)
    s$m1m1 <- if (m1 == 0) matrix() else matrix(m[(c+1):(c+m1), (c+1):(c+m1)], ncol=m1, byrow=TRUE)
    s$m1m2 <- if (m2 == 0) matrix() else matrix(m[(c+1):(c+m1), (c+1+m1):(c+m1+m2)], ncol=m2, byrow=TRUE)
    s$m2m2 <- if (m2 == 0) matrix() else matrix(m[(c+1+m1):(c+m1+m2), (c+1+m1):(c+m1+m2)], ncol=m2, byrow=TRUE)
    return (s)
  }
  # Step 1.
  sliptFixed <- fixMissingValues(slipt)
  # Step 2.
  kernel     <- kernelMatrix(sliptFixed)
  # Step 3 & 4.
  laplacian  <- graphLaplacian(kernel)
  # Step 5-1.
  #A <- chol(kernel$Y$all)
  #B <- chol(kernel$X$all)
  #return (laplacian)
}
# Step 0 : Split two parts in seeds dataset and set the missing values.
sliptAndMissDataset <- function (dataset, missingRate) {
  # sd : slipted dataset
  sd      <- dataset
  divideX <- 1:floor(median(1:ncol(sd)))
  divideY <- (floor(median(1:ncol(sd)))+1):ncol(sd)
  # C, M1, M2
  M1M2    <- sample(1:nrow(sd), (floor(nrow(sd) * missingRate)))
  M1      <- sort(sample(M1M2, (length(M1M2) * runif(1))))
  M2      <- sort(setdiff(M1M2, M1))
  C       <- sort(setdiff(1:nrow(sd), M1M2))
  # Random set NA to X or Y part.
  sd[M1, divideX] <- NA
  sd[M2, divideY] <- NA
  # Order by missing status.
  sd <- rbind(sd[C,], sd[M1,], sd[M2,])
  rownames(sd) <- 1:nrow(sd)
  # Result list.
  result <- list(X=sd[,c(divideX)], Y=sd[,c(divideY)], c=length(C), m1=length(M1), m2=length(M2))
  return (result)
}

# Read the 'seeds' data from UCI datasets.
#seeds           <- read.table('seeds_dataset.txt', sep='\t', header=FALSE)
seeds           <- read.table('seeds_less_dataset.txt', sep='\t', header=FALSE)
colnames(seeds) <- c('area', 'perimeter', 'compactness', 'length_of_kernel', 'width_of_kernel', 'asymmetry_coefficient', 'length_of_kernel_groove', 'class_of_wheat')
# Step 0.
#missingRate       <- seq(0.1, 0.9, by=0.025)
missingRate       <- 0.5
seedSlipts        <- lapply(missingRate, function (rate) { return (sliptAndMissDataset(seeds[-ncol(seeds)], rate)) })
names(seedSlipts) <- paste("missingRate_", missingRate, sep="")
# Seeds uses CoKL algorithm.
results <- lapply(seedSlipts, function (slipt) { CoKL(slipt) })