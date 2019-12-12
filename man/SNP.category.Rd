\name{SNP.category}
\alias{SNP.category}

\title{ SNP Category }

\description{ Determines a SNP Category }
\usage{ SNP.category(bed, Z, threshold = 0.8) }

\arguments{
  \item{bed}{ A bed matrix }
  \item{Z}{ A vector of length \code{nrow(bed)} }
  \item{threshold}{ Variance threshold }
}

\details{ This function determines a SNP Category from a covariable \code{Z},
which can be for example an indicator variable for a population strata, 
or the first genomic principal component.
}

\seealso{ \code{\link{qqplot.pvalues}} } 

