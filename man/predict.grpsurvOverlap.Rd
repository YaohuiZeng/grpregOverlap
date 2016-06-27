\name{predict.grpsurvOverlap}
\alias{predict.grpsurvOverlap}

\title{
Model predictions based on a fitted \code{grpsurvOverlap} object.
}
\description{
Similar to other predict methods, this function returns predictions from a fitted 
\code{"grpsurvOverlap"} object.
}

\usage{
\method{predict}{grpsurvOverlap}(object, X, type=c("link", "response", "survival",
        "median", "norm", "coefficients", "vars", "nvars", "groups", "ngroups"),
        latent = FALSE, lambda, which=1:length(object$lambda), ...)
}

\arguments{
Adoptted from \code{\link[grpreg]{predict.grpsurv}}: 
  \item{object}{Fitted \code{"grpsurvOverlap"} model object.}
  \item{X}{Matrix of values at which predictions are to be made.  Not
    used for \code{type="coefficients"} or for some of the \code{type}
    settings in \code{predict}.}
  \item{latent}{Should plot the paths of latent coefficients? Default is TRUE. Note \eqn{\hat{\gamma}} is used as the label of y axis to represent latent coefficient vector.}
  \item{lambda}{Values of the regularization parameter \code{lambda} at
    which predictions are requested.  For values of \code{lambda} not in
    the sequence of fitted models, linear interpolation is used.}
  \item{which}{Indices of the penalty parameter \code{lambda} at which
    predictions are required.  By default, all indices are returned.  If
    \code{lambda} is specified, this will override \code{which}.}
  \item{type}{Type of prediction: \code{"link"} returns the
    linear predictors; \code{"response"} gives the risk (i.e., exp(link));
    \code{"survival"} returns the estimated survival function; \code{"median"}
    estimates median survival times.  The other options are all identical to
    their \code{grpreg} counterparts: \code{"coefficients"} returns the
    coefficients; \code{"vars"} returns the indices for the nonzero coefficients;
    \code{"groups"} returns the indices for the groups with at least one
    nonzero coefficient; \code{"nvars"} returns the number of nonzero
    coefficients; \code{"ngroups"} returns the number of groups with at
    least one nonzero coefficient; \code{"norm"} returns the L2 norm of
    the coefficients in each group.}
  \item{\dots}{Not used.}
}

\details{
See \code{\link[grpreg]{predict.grpsurv}} for more details.
}      


\author{
  Yaohui Zeng and Patrick Breheny
  
  Maintainer: Yaohui Zeng <yaohui-zeng@uiowa.edu>
}
