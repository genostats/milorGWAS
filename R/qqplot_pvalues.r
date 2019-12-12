qqplot.pvalues <- function(p, snp.cat, col.cat, col.abline = "red", CB = TRUE, col.CB = "gray80", CB.level = 0.95, thinning = TRUE, ...) {
  if(missing(snp.cat))
    return( gaston::qqplot.pvalues(p, col.abline, CB, col.CB, CB.level, thinning, ...) )

  if(!is.factor(snp.cat)) 
    snp.cat <- as.factor(snp.cat)
  
  if(is.list(p)) { # ok pour les data frame aussi
    if(is.null(p$p))
      stop("No p-values were found")
    p <- p$p
  }

  w <- !is.na(p) 
  p <- p[ w ] # on supprime ces valeurs silencieusement...
  snp.cat <- snp.cat[ w ]

  if(any(p>1) | any(p<0))
    stop("p-values should be in [0,1]\n")

  mi <- round(-log10(1/length(p)))*2
  w <- (p < 10**-mi)
  if(any(w)) { # mais ça on avertit
    warning("There are ", sum(w), " p-values lower than 1e-", mi, " that will be displayed as 1e-", mi)
    p[w] <- 10**-mi
  }

  args <- list(...)
  if(is.null(args$xlab))
    args$xlab <- expression(paste("expected ", -log[10](p)))
  if(is.null(args$ylab))
    args$ylab <- expression(paste("observed ", -log[10](p)))
  if(is.null(args$main))
    args$main <- "QQ plot of p-values"

  if(missing(snp.cat)) 
    n <- length(p)
  else 
    n <- max(table(snp.cat))
 
  args$type <- "n"
  args$x <- -log10( c(1,n) / (n+1) )
  args$y <- -log10( range(p) )
  do.call( plot, args )

  # confidence interval
  k <- n * 10**seq(-log10(n), 0, length = 200)
  hi <- -log10(qbeta( 0.5-CB.level/2, k, n+1-k))
  lo <- -log10(qbeta( 0.5+CB.level/2, k, n+1-k))
  polygon( -log10(c(k, rev(k))/(n+1)), c(lo, rev(hi)), col = col.CB, border= col.CB )

  # diag line with real limits
  segments(0, 0, -log10(1/n), -log10(1/n), col = col.abline) 

  args$type <- "p"
  if(missing(col.cat))
    col.cat <- seq_len(nlevels(snp.cat))

  for(i in seq_len(nlevels(snp.cat))) {
    le <- levels(snp.cat)[i]
    # extract
    pp <- p[ snp.cat == le ]
    nn <- length(pp)
    expected <- -log10( (nn:1)/(nn+1) )
    observed <- sort(-log10(pp))
    if(thinning) {
      w <- gaston:::manhattan.thinning(expected, observed, 10000, 10000)
      args$x <- expected[w]
      args$y <- observed[w] 
    } else {
      args$x <- expected
      args$y <- observed
    }
    args$col <- col.cat[i]
    do.call( points, args )
  }
}


