\name{cv.grpregOverlap}
\alias{cv.grpregOverlap}

\title{
Cross-validation for choosing regularization parameter lambda
}

\description{
Performs k-fold cross validation for penalized regression models with overlapping grouped covariates over a grid of values for the regularization parameter lambda.
}
\usage{
cv.grpregOverlap(X, y, group, ..., nfolds = 10, seed, trace = FALSE)
}

\arguments{
  \item{X}{The design matrix, without an intercept, as in \code{grpregOverlap}.}
  \item{y}{The response vector (or matrix), as in \code{grpregOverlap}.}
  \item{group}{A list of vectors containing group information, as in \code{grpregOverlap}.}
  \item{\dots}{
Additional arguments to \code{grpregOverlap}.}
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
This function is built upon \code{\link[grpreg]{cv.grpreg}}. The class can directly call
\code{plot} function implemented for class \code{\link[grpreg]{cv.grpreg}}.
}

\value{
An object with S3 class \code{"cv.grpregOverlap"}, which inherits from class \code{"cv.grpreg"}. 
The following variables are contained in the class (adopted from \code{\link[grpreg]{cv.grpreg}}).
  \item{cve}{The error for each value of \code{lambda}, averaged across the cross-validation folds.}
  \item{cvse}{The estimated standard error associated with each value of for \code{cve}.}
  \item{lambda}{The sequence of regularization parameter values along which the cross-validation error was calculated.}
  \item{fit}{The fitted \code{grpreg} object for the whole data.}
  \item{min}{The index of \code{lambda} corresponding to \code{lambda.min}.}
  \item{lambda.min}{The value of \code{lambda} with the minimum cross-validation error.}
  \item{null.dev}{The deviance for the intercept-only model.}
  \item{pe}{If \code{family="binomial"}, the cross-validation prediction error for each value of \code{lambda}.}
}


\references{
  \itemize{
  \item Breheny P (2014). R package 'grpreg'. \url{http://cran.r-project.org/web/packages/grpreg/grpreg.pdf}
  }
}

\author{
  Yaohui Zeng and Patrick Breheny
  
  Maintainer: Yaohui Zeng <yaohui-zeng@uiowa.edu>
}


\seealso{
\code{\link{grpregOverlap}}, \code{\link{predict.grpregOverlap}}, \code{\link[=summary.cv.grpregOverlap]{summary}}, and \code{\link[grpreg]{cv.grpreg}}.
}

\examples{
## linear regression, a simulation demo.
set.seed(123)
group <- list(gr1 = c(1, 2, 3),
              gr2 = c(1, 4),
              gr3 = c(2, 4, 5),
              gr4 = c(3, 5),
              gr5 = c(6))
beta.latent.T <- c(5, 5, 5, 0, 0, 0, 0, 0, 5, 5, 0) # true latent coefficients.
# beta.T <- c(2, 3, 7, 0, 5, 0), true variables: 1, 2, 3, 5; true groups: 1, 4.
X <- matrix(rnorm(n = 6*100), ncol = 6)
X.latent <- expandX(X, group)
y <- X.latent \%*\% beta.latent.T + rnorm(100)

cvfit <- cv.grpregOverlap(X, y, group, penalty = 'grMCP')
summary(cvfit)
plot(cvfit)
par(mfrow=c(2,2))
plot(cvfit, type="all")
}

\concept{cross-validation}
