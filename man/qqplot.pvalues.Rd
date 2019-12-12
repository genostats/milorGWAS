\name{qqplot.pvalues}
\alias{qqplot.pvalues}

\title{ Stratified QQ plot of p-values }

\description{ Draws a QQ plot of p-values }
\usage{ qqplot.pvalues(p, snp.cat, col.cat, col.abline = "red", CB = TRUE, col.CB = "gray80", 
                       CB.level = 0.95, thinning = TRUE, ...) }

\arguments{
  \item{p}{ A vector of p-values, or a data.frame with a column named \code{p} }
  \item{snp.cat}{ SNP categories (optional) }
  \item{col.cat}{ Colors in which to plot SNP categories (optional) }
  \item{col.abline}{ Color of the line of slope 1. Set to \code{NA} to suppress. }
  \item{CB}{ \code{Logical}. If \code{TRUE}, a confidence band is included in the plot. }
  \item{col.CB}{ The color of the confidence band. }
  \item{CB.level}{ The level of the confidence band. }
  \item{thinning}{ \code{Logical}. If \code{TRUE}, not all points are displayed. }
  \item{...}{ Graphical parameters to be passed to \code{plot} and \code{points} }
}

\details{ This function draws a QQ plot of \eqn{p}-values, stratified by categories.
If the parameter \code{snp.cat} is missing, the function falls back on \code{gaston::qqplot.pvalues}.
}

\seealso{ \code{\link[gaston:qqplot.pvalues]{gaston:qqplot.pvalues}}, \code{\link{SNP.category}} } 

\examples{
# a random vector of categories
ca <- sample(c("A","B","C"), 1e6, TRUE, c(0.05, 0.9, 0.05))
# a vector of p-values, with different distribution depending on the strata
p <- runif(1e6)**ifelse(ca == "A", .8, ifelse(ca == "B", 1, 1.2))
qqplot.pvalues(p, ca)
}

