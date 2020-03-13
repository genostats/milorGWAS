#' SNP.category
#'
#' @param bed A bed matrix 
#' @param Z A vector of length \code{nrow(bed)} 
#' @param threshold Variance thresholds
#'
#' @details This function determines a SNP Category from a covariable \code{Z},
#' which can be for example an indicator variable for a population strata, 
#' or the first genomic principal component.
#' 
#' @return A factor giving the category of each SNP
#' 
#' @export
#' @seealso \code{\link{qqplot.pvalues}}
#' @examples
#' # a random vector of categories
#' ca <- sample(c("A","B","C"), 1e6, TRUE, c(0.05, 0.9, 0.05))
#' # a vector of p-values, with different distribution depending on the strata
#' p <- runif(1e6)**ifelse(ca == "A", .8, ifelse(ca == "B", 1, 1.2))
#' qqplot.pvalues(p, ca)
#' 

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
