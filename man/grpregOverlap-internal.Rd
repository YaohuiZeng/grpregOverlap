\name{grpregOverlap-internal}
\title{Internal functions}
\alias{gamma2beta}
\alias{incidenceMatrix}
\description{Internal functions in the package.}
\usage{
gamma2beta(gamma, incidence.mat, grp.vec) 
incidenceMatrix(X, group) 
}
\author{
Yaohui Zeng <yaohui-zeng@uiowa.edu>
}

\details{
These are not intended for use by users. \code{gamma2beta} transforms the latent coefficient matrix (or vector) into non-latent form according to the grouping information. The function \code{incidenceMatrix} returns a matrix of indicators, 1 or 0. "1" in the cell of i-th row and j-th column indicates group i contains variable j.
}

\keyword{internal}
