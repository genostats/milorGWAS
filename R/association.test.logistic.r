#' Mixed logistic regression for GWAS
#'
#' @param x a bedmatrix
#' @param Y phenotype vector. Default is column \code{pheno} of \code{x@ped}
#' @param X A matrix of covariates (defaults to a column of ones for the intercept)
#' @param K A genetic relationship matrix (or a list of such matrices)
#' @param beg Index of the first SNP tested for association
#' @param end Index of the last SNP tested for association
#' @param algorithm Algorithm to use
#' @param ... Additional parameter for \code{gaston::logistic.mm.aireml}
#'
#' @details Tests the association between the phenotype and requested SNPs in \code{x}.
#' The phenotype \code{Y} is a binary trait. A Wald test is performed using an approximate
#' method defined by the parameter \code{algorithm}.
#' 
#' @return A data frame giving for each SNP the association statistics.
#' @examples 
#' data(TTN)
#' x <- as.bed.matrix(TTN.gen, TTN.fam, TTN.bim)
#' standardize(x) <- "p"
#' ## Simulation data ##
#' set.seed(1)
#' # some covariables
#' X <- cbind(1, runif(nrow(x)))
#' # A GRM
#' ran <- random.pm( nrow(x))
#' # the random effects (tau = 1)
#' omega <- lmm.simu(1, 0, eigenK=ran$eigen)$omega
#' # linear term of the model
#' lin <- X %*% c(0.1,-0.2) + omega
#' # vector of probabilitues
#' pi <- 1/(1+exp( -lin ))
#' # vector of binary phenotypes
#' y <- rbinom(nrow(x), 1, pi)
#' # testing association with 1) the score test, 2) the offset algorithm, 3) the 'amle' algorithm
#' a1 <- association.test(x, y, X, K = ran$K, method = "lmm", response = "bin")
#' a2 <- association.test.logistic(x, y, X, K = ran$K, algorithm = "offset")
#' a3 <- association.test.logistic(x, y, X, K = ran$K, algorithm = "amle")
#' 
#' @export
association.test.logistic <- function(x, Y = x@ped$pheno, X = matrix(1, nrow(x)), K, beg = 1, end = ncol(x), 
                                      algorithm = c("offset", "amle"), ...) {

  if(beg < 1 | end > ncol(x)) stop("range too wide")
  if(is.null(x@mu) | is.null(x@p)) stop("Need mu and p to be set in x (use set.stats)")
  if(length(Y) != nrow(x)) stop("Dimensions of Y and x mismatch")
  
  X <- as.matrix(X)
  if(nrow(X) != nrow(x)) stop("Dimensions of Y and x mismatch")
  
  # check dimensions before anything
  n <- nrow(x)
  if(missing(K)) stop("argument K is mandatory")
 
  if(!is.list(K)) {
    if(n != nrow(K) | n != ncol(K))
      stop("K and x dimensions don't match")
  } else {
    if(any(n != sapply(K, nrow)) | any(n != sapply(K, ncol)))
      stop("K and x dimensions don't match")
  }

  # preparation de X [Y COMPRIS DECOMPOSITION QR]
  X <- trans.X(X, mean.y = mean(Y))

  # c'est parti
  model <- logistic.mm.aireml(Y, X = X, K, get.P = TRUE, ... )
  omega <- model$BLUP_omega   
  # ajout des effets fixes a l'offset
  if (!is.null(X)) omega <- omega + X%*%model$BLUP_beta 
  
  if(match.arg(algorithm) == "offset") { 
    t <- GWAS_logit_offset_bed(x@bed, x@p, Y, omega, X, beg-1, end-1, 1e-8, 25)
  } else {
    pi <- 1/(1+exp(-omega))
    t <- GWAS_approx_pql_bed(x@bed, Y-pi, model$P, x@p, beg-1, end-1)
  }
  
  t$p <- pchisq( (t$beta/t$sd)**2, df = 1, lower.tail = FALSE)

  # mise en forme
  L <- data.frame(chr = x@snps$chr, pos = x@snps$pos, id  = x@snps$id,  A1 = x@snps$A1, A2 = x@snps$A2, freqA2 = x@p, stringsAsFactors = FALSE)
  if(beg > 1 | end < ncol(x))  # avoid copy
    L <- L[beg:end,] 

  data.frame( c( L, t) )
}


