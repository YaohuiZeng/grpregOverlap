\name{cv.overlap.grpreg}
\alias{cv.overlap.grpreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Cross-validation for choosing regularization parameter \eqn{\lambda} for \code{overlap.grpreg}
}
\description{
Performs k-fold cross validation for penalized regression models with overlapping grouped covariates over a grid of values for the regularization parameter lambda.
}
\usage{
cv.overlap.grpreg(X, y, group, ..., nfolds = 10, seed, trace = FALSE)
}

\arguments{
  \item{X}{
The design matrix, without an intercept, as in \code{overlap.grpreg}.}
  \item{y}{
The response vector (or matrix), as in \code{overlap.grpreg}.}
  \item{group}{
A list of vectors containing group information, as in \code{overlap.grpreg}.
}
  \item{\dots}{
Additional arguments to \code{overlap.grpreg}.}
  \item{nfolds}{
The number of cross-validation folds.  Default is 10.
}
  \item{seed}{
Set the seed of the random number generator to obtain reproducible results.  
}
  \item{trace}{
If set to TRUE, print out the progress of the cross-validation. Default is FALSE.
  }
}

\details{
This function is built upon \code{cv.grpreg}.
}
\value{
An object with S3 class \code{"cv.overlap.grpreg"} and inherited from class \code{"cv.grpreg"}. 
}

\author{
Yaohui Zeng <yaohui-zeng@uiowa.edu>
}

\seealso{
\code{\link{overlap.grpreg}}, \code{\link{predict.overlap.grpreg}}, \code{\link[=summary.cv.overlap.grpreg]{summary}}
}

\examples{
## to be added
}

\concept{cross-validation}