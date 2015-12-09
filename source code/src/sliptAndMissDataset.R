# Step 0 : Split two parts in seeds dataset and set the missing values.
sliptAndMissDataset <- function (dataset, missingRate) {
  # sd : slipted dataset
  sd      <- dataset
  divideX <- 1:floor(median(1:ncol(sd)))
  divideY <- (floor(median(1:ncol(sd)))+1):ncol(sd)
  # C, M1, M2
  M1M2 <- sample(1:nrow(sd), max(floor(nrow(sd) * missingRate), 2))
  M1   <- sort(sample(M1M2, ceiling((length(M1M2) - 1) * runif(1))))
  M2   <- sort(setdiff(M1M2, M1))
  C    <- sort(setdiff(1:nrow(sd), M1M2))
  # Random set NA to X or Y part.
  sd[M2, divideX] <- NA
  sd[M1, divideY] <- NA
  # Order by missing status.
  order <- as.integer(c(rownames(sd[C,]), rownames(sd[M1,]), rownames(sd[M2,])))
  sd    <- rbind(sd[C,], sd[M1,], sd[M2,])
  rownames(sd) <- order
  # Result list.
  result <- list(X=sd[,c(divideX)], Y=sd[,c(divideY)], missingRate=missingRate, c=length(C), m1=length(M1), m2=length(M2), order=order)
  return (result)
}