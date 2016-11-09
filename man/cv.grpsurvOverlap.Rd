\name{cv.grpsurvOverlap}
\alias{cv.grpsurvOverlap}

\title{
Cross-validation for choosing regularization parameter lambda for Cox models.
}

\description{
Performs k-fold cross validation for penalized regression models with overlapping grouped covariates over a grid of values for the regularization parameter lambda.
}
\usage{
cv.grpsurvOverlap(X, y, group, ..., nfolds = 10, seed, cv.ind, 
                  returnY = FALSE, trace = FALSE)
}

\arguments{
  \item{X}{The design matrix, without an intercept, as in \code{grpregOverlap}.}
  \item{y}{The time-to-event outcome matrix for survival analysis, as explained in \code{grpregOverlap}.}
  \item{group}{A list of vectors containing group information, as in \code{grpregOverlap}.}
  \item{\dots}{
Additional arguments to \code{grpregOverlap}.}
  \item{nfolds}{
The number of cross-validation folds.  Default is 10.
}
  \item{seed}{
Set the seed of the random number generator to obtain reproducible results.  
}
  \item{cv.ind}{User specified indices of which fold each observation belongs to. By default the observations are randomly assigned.}
  \item{returnY}{Should the linear predictors from the cross-validation folds be returned?  Default is FALSE; if TRUE, this will return a matrix in which the element for row i, column j is the fitted value for observation i from the fold in which observation i was excluded from the fit, at the jth value of lambda. See details in \code{\link[grpreg]{cv.grpsurv}}}
  \item{trace}{
If set to TRUE, print out the progress of the cross-validation. Default is FALSE.
  }
}

\details{
This function is built upon \code{\link[grpreg]{cv.grpsurv}}. The \code{plot}, \code{summary}, and \code{predict} functions are also supported. See details about the cross-validation approach for fitting survival models in \code{\link[grpreg]{cv.grpsurv}}.
}

\value{
An object with S3 class \code{"cv.grpsurvOverlap"}, which inherits from class \code{"cv.grpregOverlap"} and \code{"cv.grpsurv"}. 
The following variables are contained in the class (adopted from \code{\link[grpreg]{cv.grpsurv}}).
  \item{cve}{The error for each value of \code{lambda}, averaged across the cross-validation folds.}
  \item{lambda}{The sequence of regularization parameter values along which the cross-validation error was calculated.}
  \item{fit}{The fitted \code{grpreg} object for the whole data.}
  \item{min}{The index of \code{lambda} corresponding to \code{lambda.min}.}
  \item{lambda.min}{The value of \code{lambda} with the minimum cross-validation error.}
  \item{null.dev}{The cross-validated deviance for Cox model with max(lambda). See details in \code{\link[grpreg]{cv.grpsurv}}.}
}

\references{
  \itemize{
  \item Breheny P (2014). R package 'grpreg'. \url{https://CRAN.R-project.org/package=grpreg/grpreg.pdf}
  }
}

\author{
  Yaohui Zeng and Patrick Breheny
  
  Maintainer: Yaohui Zeng <yaohui-zeng@uiowa.edu>
}


\seealso{
\code{\link{grpregOverlap}}, \code{\link{predict.grpregOverlap}}, \code{\link[=summary.cv.grpregOverlap]{summary}}, and \code{\link[grpreg]{cv.grpreg}}.
}

\concept{cross-validation}
