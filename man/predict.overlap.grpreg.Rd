\name{predict.overlap.grpreg}
\alias{predict.overlap.grpreg}
\alias{coef.overlap.grpreg}
\alias{predict.cv.overlap.grpreg}
\alias{coef.cv.overlap.grpreg}

\title{
Model predictions based on a fitted object
}
\description{
Similar to other predict methods, this function returns predictions from a fitted \code{"overlap.grpreg"} object.
}

\usage{
\method{predict}{overlap.grpreg}(object, X, type = c("link", "response", "class", "coefficients", "vars", "groups", "nvars", "ngroups", "norm"), latent = FALSE, lambda, which = 1:length(object$lambda), ...)
\method{predict}{cv.overlap.grpreg}(object, X, type = c("link", "response", "class", "coefficients", "vars", "groups", "nvars", "ngroups", "norm"), latent = FALSE, lambda = object$lambda.min, which=object$min, ...)
\method{coef}{overlap.grpreg}(object, lambda, latent = FALSE, which=1:length(object$lambda), drop=TRUE, ...)
\method{coef}{cv.overlap.grpreg}(object, latent = FALSE, lambda = object$lambda.min, which = object$min, ...)
}

\arguments{
  \item{object}{
A fitted \code{"overlap.grpreg"} or \code{"cv.overlap.grpreg"} model object.}
  \item{X}{
Matrix of values at which predictions are to be made.  Not used for \code{type="coefficients"}}
  \item{type}{
Type of prediction: \code{"link"} returns the linear predictors; \code{"response"} gives the fitted values; \code{"class"} returns the binomial outcome with the highest
probability; \code{"coefficients"} returns the coefficients; \code{"vars"} returns the indices for the nonzero coefficients; \code{"groups"} returns the indices for the groups with at least one nonzero coefficient; \code{"nvars"} returns the number of nonzero coefficients; \code{"ngroups"} returns the number of groups with at least one nonzero coefficient; \code{"norm"} returns the L2 norm of the coefficients in each group.
WARNING: if latent = 'TRUE' only for 'coefficients', 'vars', 'nvars', 'norm'...
}
  \item{latent}{
Should return prediction values at the latent level? Default is FALSE. ONLY set TRUE for \code{type} 'coefficients', 'vars', 'nvars' since for other types, there are no latent information computed. Otherwise, an warning message will be issued; If \code{type} = 'norm', this argument is overwritten to be 'TRUE' since the L2 norms are only computed at the latent level.
}
  \item{lambda}{
Values of the regularization parameter \code{lambda} at which predictions are requested.  For values of \code{lambda} not in the sequence of fitted models, linear interpolation is used.}
  \item{which}{
Indices of the penalty parameter \code{lambda} at which predictions are required.  By default, all indices are returned.  If \code{lambda} is specified, this will override \code{which}.
}
  \item{drop}{Drop the matrix to be a vector.}
  \item{\dots}{
Not used.
  }
}

\details{
\code{coef} and \code{predict} methods are provided for \code{"cv.overlap.grpreg"} options as a convenience.  They simply call \code{coef.overlap.grpreg} and \code{predict.overlap.grpreg} with \code{lambda} set to the value that minimizes the cross-validation error.
}

\value{
The object returned depends on type.
}

\author{
Yaohui Zeng <yaohui-zeng@uiowa.edu>
}

\seealso{
\code{\link{overlap.grpreg}}
}

\examples{
## To be added.
}

\keyword{predict}