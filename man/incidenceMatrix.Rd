\name{incidenceMatrix}
\alias{incidenceMatrix}

\title{
Compute the incidence matrix indicating group memebership
}
\description{
Compute a \code{p-by-J} sparse matrix indicating which group(s) the variable(s) belong(s) to.
}
\usage{
incidenceMatrix(X, group)
}

\arguments{
  \item{X}{
  The design matrix, without an intercept, as in \code{grpregOverlap}.
}
  \item{group}{
  A list of vectors containing group information, as in \code{grpregOverlap}.
}
}

\details{
This function is intended to provide an interface for users to construct the so-called 
incidence matrix, denoted as \code{M}, with dimension p-by-J, where \code{p} is the number of variables, and \code{J} is the number of groups. If variable \code{i} is contained by group
\code{j}, then \code{M[i, j] = 1}; otherwise, \code{M[i, j] = 0}.

Note that since a variable can be contained by multiple groups, so the sum of row \code{M[i, ]}
can be larger than 1.
}

\value{
A sparse \code{p-by-J} matrix as described above.
}

\author{
Yaohui Zeng <yaohui-zeng@uiowa.edu>
}

\seealso{
\code{\link{grpregOverlap}}, \code{\link{overlapMatrix}}
}

\examples{
data(pathway.dat)
X <- pathway.dat$expression
group <- pathway.dat$pathways
incidence.mat <- incidenceMatrix(X, group)
dim(incidence.mat)
}
