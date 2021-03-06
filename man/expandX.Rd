\name{expandX}
\alias{expandX}

\title{
Expand design matrix according to grouping information
}
\description{
\code{expandX} creates a new design matrix by duplicating the columns of the 
overlapped variables in design matrix \code{X}.
}
\usage{
expandX(X, group)
}

\arguments{
  \item{X}{
The design matrix, without an intercept, as in \code{grpregOverlap}.
}
  \item{group}{
A list of vectors containing group information, as in \code{grpregOverlap}.
}
}

\value{
A matrix expanded based on \code{X}, with duplicated columns corresponding to
variables being overlapped between groups.
}

\author{
Yaohui Zeng <yaohui-zeng@uiowa.edu>
}
\note{
When expanding \code{X},  The columns of \code{X} corresponding to the variables 
not included in \code{group} will be removed automatically.
}

\seealso{
\code{\link{grpregOverlap}}, \code{\link{overlapMatrix}}.
}

\examples{
set.seed(123)
group <- list(gr1 = c(1, 2, 3), gr2 = c(1, 4), gr3 = c(2, 4, 5), 
              gr4 = c(3, 5), gr5 = c(6))
beta.latent.T <- c(5, 5, 5, 0, 0, 0, 0, 0, 5, 5, 0) # true latent coefficients.
# beta.T <- c(2, 3, 7, 0, 5, 0), true variables: 1, 2, 3, 5; true groups: 1, 4.
X <- matrix(rnorm(n = 6*100), ncol = 6)  
X.latent <- expandX(X, group)
}
