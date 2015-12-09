NMI <- function (A, B) {
  if (length(A) != length(B)) {
    return ("ERROR")
  }
  total <- length(A)
  A_ids <- unique(A)
  B_ids <- unique(B)
  #
  MI <- 0
  for (A_id in A_ids) {
    for (B_id in B_ids) {
      #
      A_id_occur  <- which(A == A_id)
      B_id_occur  <- which(B == B_id)
      AB_id_occur <- intersect(A_id_occur, B_id_occur)
      #
      px  <- length(A_id_occur) / total
      py  <- length(B_id_occur) / total
      pxy <- length(AB_id_occur) / total
      #
      MI <- MI + pxy * log2(pxy / (px * py) + .Machine$double.eps)
    }
  }
  #
  Hx <- 0
  for (A_id in A_ids) {
    A_id_occur_count <- length(which(A == A_id))
    Hx <- Hx - (A_id_occur_count / total) * log2(A_id_occur_count / total + .Machine$double.eps)
  }
  Hy <- 0
  for (B_id in B_ids) {
    B_id_occur_count <- length(which(B == B_id))
    Hy <- Hy - (B_id_occur_count / total) * log2(B_id_occur_count / total + .Machine$double.eps)
  }
  #
  MIhat <- 2 * MI / (Hx + Hy)
  return (MIhat)
}