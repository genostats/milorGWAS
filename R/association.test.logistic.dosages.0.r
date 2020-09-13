association.test.logistic.dosage.0 <- function(filename, Y, X, K, beg, end,
                                             algorithm = c("amle", "offset"), omega, null.model, ...) {
  if(match.arg(algorithm) == "offset") {
    t <- GWAS_logit_offset_dosages(filename, Y, omega, X, beg-1, end-1, 1e-8, 25)
  } else {
    pi <- 1/(1+exp(-omega))
    t <- GWAS_approx_pql_dosages(filename, Y-pi, null.model$P, beg-1, end-1, 1e-8)
  }

  t$p <- pchisq( (t$beta/t$sd)**2, df = 1, lower.tail = FALSE)
  data.frame( t, stringsAsFactors = FALSE )
}

