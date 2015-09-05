\name{select.grpregOverlap}
\alias{select.grpregOverlap}
\alias{select}

\title{
Select an value of lambda along a regularization path
}
\description{
The function selects a point along the regularization path of a fitted "grpregOverlap object" according to the AIC, BIC, or GCV criteria.
}
\usage{
select(obj,\dots)

\method{select}{grpregOverlap}(obj, criterion = c("BIC", "AIC", "GCV", "AICc", "EBIC"), 
    df.method = c("default", "active"), smooth = FALSE, ...)
}

\arguments{
  \item{obj}{A fitted \code{"grpregOverlap"} model object.}
  \item{criterion}{
  The criterion by which to select the regularization parameter. One of \code{"AIC"}, \code{"BIC"}, \code{"GCV"}, \code{"AICc"}, or \code{"EBIC"}; default is \code{"BIC"}.
}
  \item{df.method}{
  How should effective model parameters be calculated? One of: \code{"active"}, which counts the number of nonzero coefficients; or \code{"default"}, which uses the calculated \code{df} returned by \code{grpregOverlap}.  Default is \code{"default"}.
}
  \item{smooth}{
  Applies a smoother to the information criteria before selecting the optimal value.
  }
  \item{\dots}{
  For S3 method compatibility.
  }
}

\details{
See R Package \code{\link[grpreg]{grpreg}} for details about the selection criteria
as well as the methods for computing degrees of freedom.
}

\value{
A list containing:
  \item{lambda}{The selected value of the regularization parameter, \code{lambda}.}
  \item{beta}{The vector of coefficients at the chosen value of \code{lambda}.}
  \item{beta.latent}{The vector of latent coefficients at the chosen value of \code{lambda}.}
  \item{df}{The effective number of model parameters at the chosen value of \code{lambda}.}
  \item{IC}{A vector of the calculated model selection criteria for each point on the regularization path.}
}

\author{
Yaohui Zeng <yaohui-zeng@uiowa.edu>
}

\seealso{
\code{\link{grpregOverlap}}
}

\examples{
## See examples in "grpregOverlap".
}
