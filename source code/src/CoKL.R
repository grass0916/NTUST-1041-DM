# CoKL algorithm.
CoKL <- function (slipt) {
  c  <- slipt$c
  m1 <- slipt$m1
  m2 <- slipt$m2
  mr <- slipt$missingRate
  or <- slipt$order
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
    kernelX <- as.matrix(s$X) %*% t(as.matrix(s$X))
    kernelY <- as.matrix(s$Y) %*% t(as.matrix(s$Y))
    kernel  <- list(X=matrixSeperator(kernelX, 2), Y=matrixSeperator(kernelY, 2))
    return (kernel)
  }
  # Step 3 & 4: Calculate Laplacian Lx and Ly.
  graphLaplacian <- function (k) {
    # k: kernel.
    diagonalX  <- diag(apply(k$X$all, 1, function (row) { sum(row) }), nrow(k$X$all), ncol(k$X$all))
    diagonalY  <- diag(apply(k$Y$all, 1, function (row) { sum(row) }), nrow(k$Y$all), ncol(k$Y$all))
    laplacianX <- diagonalX - k$X$all
    laplacianY <- diagonalY - k$Y$all
    laplacian  <- list(X=matrixSeperator(laplacianX, 2), Y=matrixSeperator(laplacianY, 2))
    return (laplacian)
  }
  matrixSeperator <- function (m, dim) {
    # s: set ; m: original matrix.
    s <- list()
    if (dim == 1) {
      s$all <- m
      s$c   <- matrix(as.matrix(m[1:c,]), ncol=ncol(m))
      s$m1  <- if (m1 == 0) matrix() else matrix(as.matrix(m[(c+1):(c+m1),]), ncol=ncol(m))
      s$m2  <- if (m2 == 0) matrix() else matrix(as.matrix(m[(c+m1+1):(c+m1+m2),]), ncol=ncol(m))
    } else if (dim == 2) {
      s$all  <- matrix(m, ncol=c+m1+m2, byrow=TRUE)
      s$cc   <- matrix(m[1:c, 1:c], ncol=c, byrow=TRUE)
      s$cm1  <- if (m1 == 0) matrix() else matrix(m[1:c, (c+1):(c+m1)], ncol=m1, byrow=TRUE)
      s$cm2  <- if (m2 == 0) matrix() else matrix(m[1:c, (c+1+m1):(c+m1+m2)], ncol=m2, byrow=TRUE)
      s$m1m1 <- if (m1 == 0) matrix() else matrix(m[(c+1):(c+m1), (c+1):(c+m1)], ncol=m1, byrow=TRUE)
      s$m1m2 <- if (m2 == 0) matrix() else matrix(m[(c+1):(c+m1), (c+1+m1):(c+m1+m2)], ncol=m2, byrow=TRUE)
      s$m2m2 <- if (m2 == 0) matrix() else matrix(m[(c+1+m1):(c+m1+m2), (c+1+m1):(c+m1+m2)], ncol=m2, byrow=TRUE)
    }
    return (s)
  }
  # Step 1.
  sliptFixed <- fixMissingValues(slipt)
  # Step 2.
  kernel     <- kernelMatrix(sliptFixed)
  # Step 3 & 4.
  laplacian  <- graphLaplacian(kernel)
  # Step 5 ~ 12.
  A <- matrixSeperator(sliptFixed$Y, 1)
  B <- matrixSeperator(sliptFixed$X, 1)
  convergenceLimit <- 0.001
  while (TRUE) {
    A$all[(c+1):(c+m1),]       <- A$m1 <- -solve(laplacian$X$m1m1) %*% (t(laplacian$X$cm1) %*% A$c - laplacian$X$m1m2 %*% A$m2)
    B$all[(c+m1+1):(c+m1+m2),] <- B$m2 <- -solve(laplacian$Y$m2m2) %*% (t(laplacian$Y$cm2) %*% B$c - t(laplacian$Y$m1m2) %*% B$m1)
    newKernel    <- kernelMatrix(list(X=B$all, Y=A$all))
    newLaplacian <- graphLaplacian(newKernel)
    # Is convergence or not (get value).
    convergenceValueX <- sqrt(sum((newKernel$X$all - kernel$X$all)^2) / (nrow(kernel$X$all) * ncol(kernel$X$all)))
    convergenceValueY <- sqrt(sum((newKernel$Y$all - kernel$Y$all)^2) / (nrow(kernel$Y$all) * ncol(kernel$Y$all)))
    # Next kernel and laplacian.
    kernel    <- newKernel
    laplacian <- newLaplacian
    # If convergence or not (process).
    if (convergenceValueX < convergenceLimit && convergenceValueY < convergenceLimit) break
  }
  kernel$missingRate <- mr
  kernel$order       <- or
  # kmeans_res <- kmeans(cbind(B$all, A$all), 3)
  # print(NMI(kmeans_res$cluster[or], kmeans(seeds[1:7], 3)$cluster))
  return (kernel)
}