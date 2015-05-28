\name{overlap.grpreg-internal}
\title{Internal functions}
\alias{expandX}
\alias{gamma2beta}
\alias{incidenceMatrix}
\description{Internal functions in the package.}
\usage{
expandX(X, incidence.mat, grp.vec) 
gamma2beta(gamma, incidence.mat, grp.vec) 
incidenceMatrix(X, group) 
}
\author{
Yaohui Zeng <yaohui-zeng@uiowa.edu>
}

We can duplicate the columns of overlapped variables
in the design matrix to form a new design matrix,

\details{
These are not intended for use by users. \code{expandX} creates a new design matrix by duplicating the columns of the overlapped variables in design matrix 'X'. \code{gamma2beta} transforms the latent coefficient matrix (or vector) into non-latent form according to the grouping information. The function \code{incidenceMatrix} returns a matrix of indicators, 1 or 0. "1" in the cell of i-th row and j-th column indicates group i contains variable j.
  }
\keyword{internal}