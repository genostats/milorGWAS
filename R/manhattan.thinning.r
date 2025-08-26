manhattan.thinning <- function(x, y, mx, my ) {
  if( is.unsorted(x) ) {
    o <- order(x);
    x1 <- x[o]
    y1 <- y[o]
    w <- .Call('_milorGWAS_manhattan_thinning', PACKAGE = "milorGWAS", x1, y1, mx, my)
    return(o[w]);
  } 
  return(.Call('_milorGWAS_manhattan_thinning', PACKAGE = "milorGWAS", x, y, mx, my))
}
