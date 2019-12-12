SNP.category <- function(bed, Z, threshold = 0.8) {
  m <- min(Z)
  M <- max(Z)
  Z <- (Z - m)/(M - m)
  # calcul de Z'G après standardisation de G (règle le cas des NA)
  standardize(bed) <- "p"
  ZG <- Z %*% bed
  n1 <- sum(Z)
  n0 <- sum(1-Z)
  # ...mais du coup il faut ruser pour retrouver q1 et q0
  q1 <- bed@p + sqrt(bed@p*(1-bed@p)/2)*(1/n1)*ZG
  q0 <- bed@p - sqrt(bed@p*(1-bed@p)/2)*(1/n0)*ZG
  var.ratio <- c(q1*(1-q1)/q0/(1-q0))
  cut( var.ratio, c(-Inf, threshold, 1/threshold, Inf) )
}
