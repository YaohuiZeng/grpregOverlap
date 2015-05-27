\name{overlap.grpreg}
\alias{overlap.grpreg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Fit penalized regression models with overlapping grouped covariates
}
\description{
Fit the regularization path of linear, logistic or poisson
  models with overlapping grouped covariates based on the latent group lasso
  approach (Jacob 2005).
}
\usage{
overlap.grpreg(X, y, group, penalty = c("grLasso", "grMCP", "grSCAD", "gel", 
"cMCP", "gBridge", "gLasso", "gMCP"), family = c("gaussian", "binomial", 
"poisson"), nlambda = 100, lambda, lambda.min = {
    if (nrow(X) > ncol(X)) 
        1e-04
    else 0.05
}, alpha = 1, eps = 0.001, max.iter = 1000, dfmax = p, gmax = J, gamma = 3, tau = 1/3, group.multiplier = {
    if (strtrim(penalty, 2) == "gr") 
        sqrt(table(group[group != 0]))
    else rep(1, J)
}, warn = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{X}{
The design matrix, without an intercept.  \code{overlap.grpreg}
    standardizes the data and includes an intercept by default.
}
  \item{y}{
The response vector, or a matrix in the case of multitask
    learning (see details).
}
  \item{group}{
Different from that in \code{grpreg}, \code{group} here must be a list of vectors, each containing indices (integers) or names (characters) of variables in a group. Variables that don't belong to any groups will be disgarded.
}
  \item{penalty}{
The penalty to be applied to the model.  For group selection, one of \code{grLasso}, \code{grMCP}, or \code{grSCAD}. For bi-level selection, one of \code{gel} or \code{cMCP}.  See Package 'grpreg' for details.
}
  \item{family}{
Either "gaussian", "binomial", or 'poisson', depending on the response.
}
  \item{nlambda}{
The number of \code{lambda} values.  Default is 100.
}
  \item{lambda}{
A user supplied sequence of \code{lambda} values. Typically, this is left unspecified, and the function automatically computes a grid of lambda values that ranges uniformly on the log scale over the relevant range of lambda values.
}
  \item{lambda.min}{The smallest value for \code{lambda}, as a fraction
    of \code{lambda.max}.  Default is .0001 if the number of observations
    is larger than the number of covariates and .05 otherwise.} 
  \item{alpha}{\code{grpreg} allows for both a group penalty and an L2
    (ridge) penalty; \code{alpha} controls the proportional weight of
    the regularization parameters of these two penalties.  The group
    penalties' regularization parameter is \code{lambda*alpha}, while
    the regularization parameter of the ridge penalty is
    \code{lambda*(1-alpha)}.  Default is 1: no ridge penalty.}
  \item{eps}{Convergence threshhold.  The algorithm iterates until the
    change (on the standardized scale) in any coefficient is less than
    \code{eps}.  Default is \code{.001}.  See details.}
  \item{max.iter}{Maximum number of iterations.  Default is 1000.  See
    details.}
  \item{dfmax}{Limit on the number of parameters allowed to be nonzero.
    If this limit is exceeded, the algorithm will exit early from the
    regularization path.}
  \item{gmax}{Limit on the number of groups allowed to have nonzero
    elements.  If this limit is exceeded, the algorithm will exit early
    from the regularization path.}
  \item{gamma}{Tuning parameter of the MCP penalty; defaults to 3.}
  \item{tau}{Tuning parameter for the group exponential lasso; defaults
    to 1/3.}
  \item{group.multiplier}{A vector of values representing multiplicative
    factors by which each group's penalty is to be multiplied.
    Often, this is a function (such as the square root) of the number of
    predictors in each group.  The default is to use the square root of
    group size for the group selection methods, and a vector of 1's
    (i.e., no adjustment for group size) for bi-level selection.}
  \item{warn}{Should the function give a warning if it fails to
    converge?  Default is TRUE.  See details.}
  \item{...}{Arguments passed to other functions (such as gBridge).}
}
\details{
Introduction to latent group lasso selection for overlapping grouped covariates. Explain latent coefficients.
}
\value{
An object with S3 class \code{"overlap.grpreg", "grpreg"} containing:
  \item{beta}{The fitted matrix of coefficients.  The number of rows is
    equal to the number of coefficients, and the number of columns is
    equal to \code{nlambda}.}
  \item{family}{Same as above.}
  \item{group}{Same as above.}
  \item{lambda}{The sequence of \code{lambda} values in the path.}
  \item{alpha}{Same as above.}
  \item{loss}{A vector containing either the residual sum of squares
  (\code{"gaussian"}) or negative log-likelihood (\code{"binomial"}) of
  the fitted model at each value of \code{lambda}.}
  \item{n}{Number of observations.}
  \item{penalty}{Same as above.}
  \item{df}{A vector of length \code{nlambda} containing estimates of
    effective number of model parameters all the points along the
    regularization path.  For details on how this is calculated, see
    Breheny and Huang (2009).}
  \item{iter}{A vector of length \code{nlambda} containing the number
    of iterations until convergence at each value of \code{lambda}.}
  \item{group.multiplier}{A named vector containing the multiplicative
    constant applied to each group's penalty.}
  \item{beta.latent}{The fitted matrix of latent coefficients.}
  \item{indicator.mat}{indicator matrix.}
  \item{incidence.mat}{incidence matrix.}
  \item{X.latent}{The expanded design matrix for latent formulation.}
}
\references{
\itemize{
    \item Jacob, L., Obozinski, G., & Vert, J. P. (2009, June). Group lasso with overlap and graph lasso. \emph{In Proceedings of the 26th annual international conference on machine learning, ACM}: 433-440.
    \item Breheny, P. and Huang, J. (2009) Penalized methods for bi-level variable selection.  \emph{Statistics and its interface}, \strong{2}: 369-380. \url{myweb.uiowa.edu/pbreheny/publications/Breheny2009.pdf}
    \item Huang J., Breheny, P. and Ma, S. (2012). A selective review of group selection in high dimensional models. \emph{Statistical Science}, \strong{27}: 481-499. \url{myweb.uiowa.edu/pbreheny/publications/Huang2012.pdf}
    \item Breheny, P. and Huang, J. (to appear) Group descent algorithms for nonconvex penalized linear and logistic regression models with grouped predictors. \emph{Statistics and Computing}, to appear. \url{www.springerlink.com/openurl.asp?genre=article&id=doi:10.1007/s11222-013-9424-2}
    \item Breheny, P. The group exponential lasso for bi-level variable selection. In submission. \url{myweb.uiowa.edu/pbreheny/publications/grpexp.pdf}
  }
}
\author{
Yaohui Zeng <yaohui-zeng@uiowa.edu>
}

\seealso{
\code{\link{cv.overlap.grpreg}}, \code{\link[=plot.overlap.grpreg]{plot}}.
}

\examples{
## to be added
}

\keyword{overlap.grpreg}
\keyword{models}
