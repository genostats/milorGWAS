# Check AND replaces by QR decomposition...
trans.X <- function(X, PCs = matrix(0, nrow=nrow(X), ncol=0), mean.y = 1) {
  if(any(is.na(X)))
    stop("Covariates can't be NA")

  PCs <- as.matrix(PCs) # cas où p = 1
  n.X  <- ncol(X)
  n.pc <- ncol(PCs)
  n <- n.X + n.pc

  qr.X <- qr( cbind(PCs, X) );
  if(qr.X$rank < n) {
    warning("Covariate matrix X is not full rank, removing col(s) ", paste(qr.X$pivot[ seq(qr.X$rank+1,n) ] - n.pc , collapse = ", "))
    X <- X[ , qr.X$pivot[seq(n.pc+1, qr.X$rank)] - n.pc]
    qr.X <- qr(X)
  }
  if(mean.y > 1e-4) {
    X1 <- cbind(1,X);
    qr.X1 <- qr(X1);
    if(qr.X1$rank == ncol(X1)) {
      warning("An intercept column was added to the covariate matrix X")
      X <- X1;
      qr.X <- qr.X1
    }
  }
  if( qr.X$rank == ncol(X) )
    qr.Q(qr.X)
  else
    qr.Q(qr(X))
}

