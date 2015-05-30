\name{plot.overlap.grpreg}
\alias{plot.overlap.grpreg}

\title{
Plot object "overlap.grpreg"
}
\description{
Plot the paths of estimated coefficients, latent coefficients (default), or L2 
norms of the groups for a fitted object "overlap.grpreg".
}
\usage{
\method{plot}{overlap.grpreg}(x, legend.loc, alpha = 1, latent = TRUE, log.l = FALSE, norm = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{The fitted \code{"overlap.grpreg"} object.}
  \item{legend.loc}{
  Where to put the legend?  If left unspecified, no legend is drawn.  
  See \code{legend} for details.}
  \item{alpha}{Controls alpha-blending.  Default is alpha=1.}
  \item{latent}{
  Should plot the paths of latent coefficients? Default is TRUE. Note \eqn{\hat{\gamma}}
  is used as the label of y axis to represent latent coefficient vector.
  }
  \item{log.l}{Should horizontal axis be on the log scale?  Default is FALSE.}
  \item{norm}{
  If \code{TRUE}, plot the norm of each group, rather than the individual 
  coefficients. Note that the norm is of latent coefficients.
  }
  \item{\dots}{
  Other graphical parameters to \code{plot}, \code{matlines}, or \code{legend}.
  }
}

\author{
Yaohui Zeng <yaohui-zeng@uiowa.edu>
}

\seealso{
\code{\link{overlap.grpreg}}
}

\examples{
## See examples in 'overlap.grpreg'
}

\keyword{graphics}