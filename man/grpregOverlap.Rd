\name{grpregOverlap}
\alias{grpregOverlap}
\title{
Fit penalized regression models with overlapping grouped variables
}

\description{
Fit the regularization paths of linear, logistic, Poisson or Cox models with 
overlapping grouped covariates based on the latent group lasso approach 
(Jacob et al., 2009; Obozinski et al., 2011). Latent group MCP/SCAD as well as 
bi-level selection methods, namely the group exponential lasso (Breheny, 2015) 
and the composite MCP (Huang et al., 2012) are also available.

This function is a useful wrapper to the \code{grpreg} package's
\code{grpreg} and \code{grpsurv} (depending on the \code{family})
functions. Arguments can be passed through to these functions using
\code{...}, see \code{\link[grpreg]{grpreg}} and
\code{\link[grpreg]{grpsurv}} for usage and more details.
}

\usage{
grpregOverlap(X, y, group, 
    family=c("gaussian","binomial", "poisson", "cox"),
    returnX.latent = FALSE, returnOverlap = FALSE, ...)
}

\arguments{
  \item{X}{
  The design matrix, without an intercept.  \code{grpregOverlap} calls 
  \code{grpreg}, which standardizes the data and includes an intercept by default.
  }
  \item{y}{
  The response vector, or a matrix in the case of multitask learning. For survival analysis, \code{y} is the time-to-event outcome - a two-column matrix or
    \code{\link[survival]{Surv}} object.  The first column is the
    time on study (follow up time); the second column is a binary
    variable with 1 indicating that the event has occurred and 0
    indicating (right) censoring. See \code{\link[grpreg]{grpreg}} and  \code{\link[grpreg]{grpsurv}} for more details.
  }
  \item{group}{
  Different from that in \code{grpreg}, \code{group} here must be a list of vectors,
  each containing integer indices or  character names of variables in the group. 
  variables that not belong to any groups will be disgarded.
  }
  \item{family}{
  Either "gaussian", "binomial", or 'cox', depending on the response. If 
  \code{family} is missing, it is set to be 'gaussian'. Specify \code{family} = 'cox' for survival analysis (Cox models).
  }
  \item{returnX.latent}{
  Return the new expanded design matrix? Default is FALSE. Note the storage size 
  of this new matrix can be very large. Note: the name of this argument
  was recently changed so that returnX can be passed through to
  \code{\link[grpreg]{grpreg}} (in which case it will return the
  group-orthonormalized design.
  }
  \item{returnOverlap}{
  Return the matrix containing overlapps? Default is FALSE. It is a square matrix
  \eqn{C} such that \eqn{C[i, j]} is the number of overlapped variables between 
  group i and j. Diagonal value \eqn{C[i, i]} is therefore the number of 
  variables in group i.
  }
  \item{...}{
  Used to pass options (e.g., `group.multiplier`) to
  \code{\link[grpreg]{grpreg}}. Note: the \code{returnX} argument will
  not be passed through, since this will cause \code{grpregOverlap} to
  store X.latent in the fitted model object.
  }
}

\details{
The latent group lasso approach extends the group lasso to group variable selection 
with overlaps. The proposed \emph{latent group lasso} penalty is formulated in a 
way such that it's equivalent to a classical non-overlapping group lasso problem 
in an new space, which is expanded by duplicating the columns of overlapped variables.
For technical details, see (Jacob et al., 2009) and (Obozinski et al., 2011).

\code{grpregOverlap} takes input design matrix \code{X} and grouping information
\code{group}, and expands {X} to the new, non-overlapping space. It then calls
\code{grpreg} for modeling fitting based on group decent algorithm. Unlike
in \code{grpreg}, the interface for group bridge-penalized method is not implemented.

The expanded design matrix is named \code{X.latent}. It is a returned value in the fitted
object, provided \code{returnX.latent} is TRUE. The latent coeffecient (or norm) vector then 
corresponds to that. Note thaT when constructing \code{X.latent}, the columns in \code{X} 
corresponding to those variables not included in \code{group} will be removed automatically.

For more detailed explanation for the penalties and algorithm, see \code{\link[grpreg]{grpreg}}.
}

\value{
An object with S3 class \code{"grpregOverlap"} or \code{"grpsurvOverlap"} (for Cox models), which inherits \code{"grpreg"}, 
with following variables.
  \item{beta}{
  The fitted matrix of coefficients. The number of rows is equal to the number 
  of coefficients, and the number of columns is equal to \code{nlambda}.
  }
  \item{family}{Same as above.}
  \item{group}{Same as above.}
  \item{lambda}{
  The sequence of \code{lambda} values in the path.
  }
  \item{alpha}{Same as above.}
  \item{loss}{
  A vector containing either the residual sum of squares (\code{"gaussian"}) or 
  negative log-likelihood (\code{"binomial"}) or negative partial log-likelihood (\code{"cox"}) of the fitted model at each value of \code{lambda}.}
  \item{n}{Number of observations.}
  \item{penalty}{Same as above.}
  \item{df}{
  A vector of length \code{nlambda} containing estimates of effective 
  number of model parameters all the points along the regularization path.
  For details on how this is calculated, see Breheny and Huang (2009).
  }
  \item{iter}{
  A vector of length \code{nlambda} containing the number of iterations until 
  convergence at each value of \code{lambda}.
  }
  \item{group.multiplier}{
  A named vector containing the multiplicative constant applied to each group's 
  penalty.
  }
  \item{beta.latent}{
  The fitted matrix of latent coefficients. The number of rows is equal to the number 
  of coefficients, and the number of columns is equal to \code{nlambda}.
  }
  \item{incidence.mat}{
  Incidence matrix: I[i, j] = 1 if group i contains variable j; otherwise 0.
  }
  \item{grp.vec}{
  A vector of consecutive integers indicating grouping information of variables. This
  is equivalent to argument \code{group} in \code{\link[grpreg]{grpreg}}.
  }
  \item{overlap.mat}{
  A square matrix \eqn{C} where \eqn{C[i, j]} is the number of overlapped 
  variables between group i and j. Diagonal value \eqn{C[i, i]} is therefore the 
  number of variables in group i. Only returned if \code{returnOverlap} is TRUE.
  }
  \item{X.latent}{
  The new expanded design matrix for the latent group lasso formulation. The
  variables are reordered according to the order of groups. Only returned if
  \code{returnX.latent} is TRUE.
  }
  \item{W}{Matrix of \code{exp(beta)} values for each subject over all
    \code{lambda} values. (For Cox models only)}
  \item{time}{Times on study. (For Cox models only)}
  \item{fail}{Failure event indicator. (For Cox models only)}
}
\references{
  \itemize{
    \item Zeng, Y., and Breheny, P. (2016). Overlapping Group Logistic Regression with Applications to Genetic Pathway Selection. \emph{Cancer Informatics}, \strong{15}, 179-187. \url{http://doi.org/10.4137/CIN.S40043}. 
    \item Jacob, L., Obozinski, G., and Vert, J. P. (2009, June). Group lasso with overlap and graph lasso. \emph{In Proceedings of the 26th annual international conference on machine learning, ACM}: 433-440. \url{http://www.machinelearning.org/archive/icml2009/papers/471.pdf}
    \item Obozinski, G., Jacob, L., and Vert, J. P. (2011). Group lasso with overlaps: the latent group lasso approach. \url{http://arxiv.org/abs/1110.0413}.
    \item Breheny, P. and Huang, J. (2009) Penalized methods for bi-level variable selection.  \emph{Statistics and its interface}, \strong{2}: 369-380. \url{http://myweb.uiowa.edu/pbreheny/publications/Breheny2009.pdf}
    \item Huang J., Breheny, P. and Ma, S. (2012). A selective review of group selection in high dimensional models. \emph{Statistical Science}, \strong{27}: 481-499. \url{http://myweb.uiowa.edu/pbreheny/publications/Huang2012.pdf}
    \item Breheny P and Huang J (2015). Group descent algorithms for nonconvex penalized linear and logistic regression models with grouped predictors. \emph{Statistics and Computing}, \strong{25}: 173-187.\url{http://myweb.uiowa.edu/pbreheny/publications/group-computing.pdf}
    \item Breheny P and Huang J (2009). Penalized methods for bi-level variable selection. \emph{Statistics and Its Interface}, \strong{2}: 369-380. \url{http://myweb.uiowa.edu/pbreheny/publications/Breheny2009.pdf}
    \item Breheny P (2014). R package 'grpreg'. \url{https://CRAN.R-project.org/package=grpreg/grpreg.pdf}
  }
}

\author{
  Yaohui Zeng and Patrick Breheny
  
  Maintainer: Yaohui Zeng <yaohui-zeng@uiowa.edu>
}

\seealso{
\code{\link{cv.grpregOverlap}}, \code{\link{cv.grpsurvOverlap}}, \code{\link[=plot.grpregOverlap]{plot}}, 
\code{\link[=select.grpregOverlap]{select}}, \code{\link[grpreg]{grpreg}}, \code{\link[grpreg]{grpsurv}}.
}

\examples{
## linear regression, a simulation demo.
set.seed(123)
group <- list(gr1 = c(1, 2, 3), gr2 = c(1, 4), gr3 = c(2, 4, 5), 
              gr4 = c(3, 5), gr5 = c(6))
beta.latent.T <- c(5, 5, 5, 0, 0, 0, 0, 0, 5, 5, 0) # true latent coefficients.
# beta.T <- c(5, 5, 10, 0, 5, 0), true variables: 1, 2, 3, 5; true groups: 1, 4.
X <- matrix(rnorm(n = 6*100), ncol = 6)  
X.latent <- expandX(X, group)
y <- X.latent \%*\% beta.latent.T + rnorm(100)

fit <- grpregOverlap(X, y, group, penalty = 'grLasso')
# fit <- grpregOverlap(X, y, group, penalty = 'grMCP')
# fit <- grpregOverlap(X, y, group, penalty = 'grSCAD')
head(coef(fit, latent = TRUE)) # compare to beta.latent.T
plot(fit, latent = TRUE) 
head(coef(fit, latent = FALSE)) # compare to beta.T
plot(fit, latent = FALSE)

cvfit <- cv.grpregOverlap(X, y, group, penalty = 'grMCP')
plot(cvfit)
head(coef(cvfit))
summary(cvfit)

## logistic regression, real data, pathway selection
data(pathway.dat)
X <- pathway.dat$expression
group <- pathway.dat$pathways
y <- pathway.dat$mutation
fit <- grpregOverlap(X, y, group, penalty = 'grLasso', family = 'binomial')
plot(fit)
str(select(fit))
str(select(fit,criterion="AIC",df="active"))

\dontrun{
cvfit <- cv.grpregOverlap(X, y, group, penalty = 'grLasso', family = 'binomial')
coef(cvfit)
predict(cvfit, X, type='response')
predict(cvfit, X, type = 'class')
plot(cvfit)
plot(cvfit, type = 'all')
summary(cvfit)
}
}

\keyword{grpregOverlap}
\keyword{models}
