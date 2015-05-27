\name{overlap.grpreg-internal}
\title{Internal functions}
\alias{expandX}
\alias{gamma2beta}
\alias{indicatorMatrix}
\description{Internal functions in the package.}
\usage{
expandX(X, indicator.mat, grp.vec) 
gamma2beta(gamma, indicator.mat, grp.vec) 
indicatorMatrix(X, group) 
}
\author{
Yaohui Zeng <yaohui-zeng@uiowa.edu>
}

\details{These are not intended for use by users. \code{expandX} creates a new design matrix by duplicating the columns of the overlapped variables. \code{gamma2beta} transforms the latent coefficient matrix (or vector) into non-latent coefficients according to the group information.The function \code{indicatorMatrix} returns a matrix of indicators (1 or 0). "1" in the cell of i-th row and j-th column indicates group i contains variable j.
  }
\keyword{internal}