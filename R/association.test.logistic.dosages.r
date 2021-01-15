#' Mixed logistic regression for GWAS, using dosages
#' 
#' @param filename Name of a dosage file
#' @param Y phenotype vector. Default is column \code{pheno} of \code{x@ped}
#' @param X A matrix of covariates (defaults to a column of ones for the intercept)
#' @param K A genetic relationship matrix (or a list of such matrices)
#' @param beg Index of the first SNP tested for association
#' @param end Index of the last SNP tested for association
#' @param algorithm Algorithm to use
#' @param eigenK eigen decomposition of K (only if \code{p} > 0)
#' @param p Number of principal components to include in the model
#' @param ... Additional parameter for \code{gaston::logistic.mm.aireml}
#'
#' @details Dosage files can be VCF files with 'DS' or 'GP' fields. It is also possible to use a file
#' with columns 'id"', 'chr', 'pos', 'A1', 'A2', 'sample1', 'sample2', etc. These files should have a 
#' header with column names. 
#' @details For more details refer to \code{\link{association.test.logistic}} and \code{\link[gaston]{association.test}}.
#' 
#' @return A data frame giving for each SNP the association statistics.
#' 
#' @seealso  \code{\link{association.test.logistic}}, \code{\link[gaston]{association.test}}
#' 
#' @export
association.test.logistic.dosage <- function(filename, Y, X, K, beg, end,
                                             algorithm = c("amle", "offset"), eigenK, p = 0, 
                                             n.cores = 1L, ...) {
  filename <- path.expand(filename)
  dims <- dim.dosage.file(filename)
  nb.inds <- dims[1]
  nb.snps <- dims[2]

  if(missing(beg)) beg <- 1
  if(missing(end)) end <- nb.snps

  if(beg < 1 || end > nb.snps) stop("range too wide")
  if(length(Y) != nb.inds) stop("Dimension of Y and #individuals in ", filename, " mismatch")

  if(!is.list(K)) {
    if(nb.inds != nrow(K) | nb.inds != ncol(K))
      stop("K and x dimensions don't match")
  } else {
    if(any(nb.inds != sapply(K, nrow)) | any(nb.inds != sapply(K, ncol)))
      stop("K and x dimensions don't match")
  }

  if(missing(X)) X <- rep(1, nb.inds); # default = intercept 
  X <- as.matrix(X)

  if(nrow(X) != nb.inds) stop("Dimensions of Y and #individuals in ", filename, " mismatch")
  # preparation de X [Y COMPRIS DECOMPOSITION QR]
  if(p > 0)
    X <- cbind(X, eigenK$vectors[,seq_len(p)])
  X <- trans.X(X, mean.y = mean(Y))

  # c'est parti
  null.model <- logistic.mm.aireml(Y, X = X, K, get.P = TRUE, ... )
  omega <- null.model$BLUP_omega
  # ajout des effets fixes a l'offset
  if (!is.null(X)) omega <- omega + X %*% null.model$BLUP_beta

  arg <- list(filename = filename, Y = Y, X = X, K = K, algorithm = algorithm)
  arg$null.model <- null.model
  arg$omega <- omega

  if(n.cores == 1) {
    arg$beg = beg;
    arg$end = end;
    return(do.call(association.test.logistic.dosage.0, arg))
  }
  if(.Platform$OS.type != "unix") {
    stop("FORK cluster unavailable")
  }
  a <- round(seq(beg, end, length = n.cores + 1))
  BEG <- a[1:n.cores]
  END <- a[-1]-1;
  END[n.cores] <- end;

  ARG <- rep( list(arg), n.cores)
  ARG <- mapply(function(a, be, en) {a$beg = be; a$end = en; a}, ARG, BEG, END, SIMPLIFY=FALSE )
  cl <- parallel::makeForkCluster(n.cores)
  xx <- parallel::clusterApply(cl, ARG, function(a) do.call(association.test.logistic.dosage.0, a))
  parallel::stopCluster(cl)
  Reduce(rbind, xx)
}

