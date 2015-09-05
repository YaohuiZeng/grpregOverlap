\name{grpregOverlap}
\alias{grpregOverlap}
\title{
Fit penalized regression models with overlapping grouped variables
}

\description{
Fit the regularization path of linear, logistic or poisson models with 
overlapping grouped covariates based on the latent group lasso approach 
(Jacob et al., 2009; Obozinski et al., 2011). Latent group MCP/SCAD as well as 
bi-level selection methods, namely the group exponential lasso (Breheny, 2015) 
and the composite MCP (Huang et al., 2012) are also available.
}

\usage{
grpregOverlap(X, y, group, 
    penalty=c("grLasso", "grMCP", "grSCAD", "gel", "cMCP", "gLasso", "gMCP"), 
    family=c("gaussian","binomial", "poisson"), nlambda=100, lambda, 
    lambda.min={if (nrow(X) > ncol(X)) 1e-4 else .05}, alpha=1, eps=.001, 
    max.iter=1000, dfmax=ncol(X), gmax=length(group), gamma=3, tau=1/3, 
    group.multiplier={if (strtrim(penalty,2)=="gr")  sqrt(sapply(group, length)) 
                         else rep(1, length(group))}, 
    warn=TRUE, ...) 
}

\arguments{
  \item{X}{
  The design matrix, without an intercept.  \code{grpregOverlap} calls 
  \code{grpreg}, which standardizes the data and includes an intercept by default.
  }
  \item{y}{
  The response vector, or a matrix in the case of multitask learning. See 
  \code{\link[grpreg]{grpreg}} for more details.
  }
  \item{group}{
  Different from that in \code{grpreg}, \code{group} here must be a list of vectors,
  each containing integer indices or  character names of variables in the group. 
  variables that not belong to any groups will be disgarded.
  }
  \item{penalty}{
  The penalty to be applied to the model.  Specify \code{grLasso}, \code{grMCP}, or
  \code{grSCAD} for group selection. Or specify \code{gel} or \code{cMCP} for 
  bi-level selection, i.e., selecting important groups as well as important 
  variables in those groups. See \code{\link[grpreg]{grpreg}} for more details.
  }
  \item{family}{
  Either "gaussian", "binomial", or 'poisson', depending on the response. If 
  \code{family} is missing, it is set to be 'gaussian'.
  }
  \item{nlambda}{
  The number of \code{lambda} values. Default is 100.
  }
  \item{lambda}{
  A user supplied sequence of \code{lambda} values. Typically, this is left 
  unspecified, and the function automatically computes a grid of lambda values 
  that ranges uniformly on the log scale over the relevant range of lambda values.
  }
  \item{lambda.min}{
  The smallest value for \code{lambda}, as a fraction of \code{lambda.max}.
  Default is .0001 if the number of observations is larger than the number of 
  covariates and .05 otherwise.
  } 
  \item{alpha}{
  Adopted from \code{grpreg}, the L2 (ridge) penalty is also allowed along with 
  the group penalty. \code{alpha} controls the proportional weight of the 
  regularization parameters of these two penalties. The regularization parameter 
  of the group penalty is \code{lambda*alpha}, while that of the ridge penalty is
  \code{lambda*(1-alpha)}. Default is 1: no L2 penalty.
  }
  \item{eps}{Convergence threshhold.  The algorithm iterates until the change (on the standardized scale) in any coefficient is less than \code{eps}. Default is \code{.001}.
  }
  \item{max.iter}{
  The maximum number of iterations. Default is 1000. See \code{\link[grpreg]{grpreg}} 
  for more details.
  }
  \item{dfmax}{
  Limit on the number of parameters allowed to be nonzero. If this limit 
  is exceeded, the algorithm will exit early from the regularization path.
  Default is the total number of covariates.
  }
  \item{gmax}{
  Limit on the number of groups allowed to have nonzero elements. If this limit 
  is exceeded, the algorithm will exit early from the regularization path.
  Default is the total number of groups.
  }
  \item{gamma}{
  Tuning parameter of the MCP penalty; defaults to 3.
  }
  \item{tau}{
  Tuning parameter for the group exponential lasso; defaults to 1/3.
  }
  \item{group.multiplier}{
  A vector of values representing multiplicative factors by which each group's 
  penalty is to be multiplied. Often, this is a function (such as the square root) 
  of the number of predictors in each group. The default is to use the square root 
  of group size for the group selection methods, and a vector of 1's (i.e., no 
  adjustment for group size) for bi-level selection.
  }
  \item{warn}{
  Should the function give a warning if it fails to converge?  Default is TRUE.  
  See \code{\link[grpreg]{grpreg}} for more details.
  }
  \item{...}{
  Not used currently. 
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
\code{grpreg} for modeling fitting based on group coordinate decent algorithm. Unlike
in \code{grpreg}, the interface for group bridge-penalized method is not implemented
yet, and will be adopted in later version.

The expanded design matrix is named \code{X.latent} as a returned value in the fitted
object. The latent coeffecient (or norm) vector then corresponds to that. Note that
When forming \code{X.latent}, the columns in \code{X} corresponding to those variables 
not included in \code{group} will be removed automatically.

For more detailed explanation to the penalties and algorithm, see \code{\link[grpreg]{grpreg}}.
}

\value{
An object with S3 class \code{"grpregOverlap"}, which inherits \code{"grpreg"}, 
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
  negative log-likelihood (\code{"binomial"}) of the fitted model at each value 
  of \code{lambda}.}
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
  number of variables in group i.
  }
  \item{X.latent}{
  The new expanded design matrix for the latent group lasso formulation. The
  variables are reordered according to the order of groups.
  }
}
\references{
  \itemize{
    \item Jacob, L., Obozinski, G., and Vert, J. P. (2009, June). Group lasso with overlap and graph lasso. \emph{In Proceedings of the 26th annual international conference on machine learning, ACM}: 433-440. \url{http://dl.acm.org/citation.cfm?id=1553431}
    \item Obozinski, G., Jacob, L., and Vert, J. P. (2011). Group lasso with overlaps: the latent group lasso approach. \url{http://arxiv.org/abs/1110.0413}.
    \item Breheny, P. and Huang, J. (2009) Penalized methods for bi-level variable selection.  \emph{Statistics and its interface}, \strong{2}: 369-380. \url{http://myweb.uiowa.edu/pbreheny/publications/Breheny2009.pdf}
    \item Huang J., Breheny, P. and Ma, S. (2012). A selective review of group selection in high dimensional models. \emph{Statistical Science}, \strong{27}: 481-499. \url{http://myweb.uiowa.edu/pbreheny/publications/Huang2012.pdf}
    \item Breheny P and Huang J (2015). Group descent algorithms for nonconvex penalized linear and logistic regression models with grouped predictors. \emph{Statistics and Computing}, \strong{25}: 173-187.\url{http://myweb.uiowa.edu/pbreheny/publications/group-computing.pdf}
    \item Breheny P and Huang J (2009). Penalized methods for bi-level variable selection. \emph{Statistics and Its Interface}, \strong{2}: 369-380. \url{http://myweb.uiowa.edu/pbreheny/publications/Breheny2009.pdf}
    \item Breheny P (2014). R package 'grpreg'. \url{http://cran.r-project.org/web/packages/grpreg/grpreg.pdf}
  }
}

\author{
Yaohui Zeng <yaohui-zeng@uiowa.edu>
}

\seealso{
\code{\link{cv.grpregOverlap}}, \code{\link[=plot.grpregOverlap]{plot}}, 
\code{\link[=select.grpregOverlap]{select}}, and \code{\link[grpreg]{grpreg}}.
}

\examples{
## linear regression, a simulation demo.
set.seed(123)
group <- list(gr1 = c(1, 2, 3), gr2 = c(1, 4), gr3 = c(2, 4, 5), 
              gr4 = c(3, 5), gr5 = c(6))
beta.latent.T <- c(5, 5, 5, 0, 0, 0, 0, 0, 5, 5, 0) # true latent coefficients.
# beta.T <- c(2, 3, 7, 0, 5, 0), true variables: 1, 2, 3, 5; true groups: 1, 4.
X <- matrix(rnorm(n = 6*100), ncol = 6)  
X.latent <- expandX(X, group)
y <- X.latent \%*\% beta.latent.T + rnorm(100)

fit <- grpregOverlap(X, y, group, penalty = 'grLasso')
# fit <- grpregOverlap(X, y, group, penalty = 'grMCP')
# fit <- grpregOverlap(X, y, group, penalty = 'grSCAD')
coef(fit, latent = TRUE) # compare to beta.latent.T
plot(fit, latent = TRUE) 
coef(fit, latent = FALSE) # compare to beta.T
plot(fit, latent = FALSE)

cvfit <- cv.grpregOverlap(X, y, group, penalty = 'grMCP')
plot(cvfit)
coef(cvfit)
summary(cvfit)

## logistic regression, real data, pathway selection
data(pathway.dat)
X <- pathway.dat$expression
group <- pathway.dat$pathways
y <- pathway.dat$mutation
fit <- grpregOverlap(X, y, group, penalty = 'grLasso', family = 'binomial')
plot(fit)
select(fit)
select(fit,criterion="AIC",df="active")

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
