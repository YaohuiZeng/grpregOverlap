\name{select.overlap.grpreg}
\alias{select.overlap.grpreg}
\alias{select}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Select an value of lambda along a regularization path
}
\description{
The function selects a point along the regularization path of a fitted 'overlap.grpreg object' according to the AIC, BIC, or GCV criteria. (Experimental!)
}
\usage{
select(obj,\dots)
\method{select}{overlap.grpreg}(obj, criterion = c("BIC", "AIC", "GCV", "AICc", "EBIC"), df.method = c("default", "active"), smooth = FALSE, ...)
}

\arguments{
  \item{obj}{
A fitted \code{"overlap.grpreg"} model object.
}
  \item{criterion}{
The criterion by which to select the regularization parameter. One of \code{"AIC"}, \code{"BIC"}, \code{"GCV"}, \code{"AICc"}, or \code{"EBIC"}; default is \code{"BIC"}.
}
  \item{df.method}{
How should effective model parameters be calculated? One of: \code{"active"}, which counts the number of nonzero coefficients; or \code{"default"}, which uses the calculated \code{df} returned by \code{grpreg}.  Default is \code{"default"}.
}
  \item{smooth}{
Applies a smoother to the information criteria before selecting the optimal value.
  }
  \item{\dots}{
For S3 method compatibility.
  }
}

\details{
See R Package 'grpreg' for details.
%%The criteria are defined as follows, where \eqn{\ell}{L} is the log
%%  likelihood, \eqn{\nu}{df} is the degrees of freedom, and \eqn{n}{n} is
%%  the sample size: \cr
%%  \deqn{AIC = 2L + 2df}{AIC = 2*L + 2*df}
%%  \deqn{BIC = 2*L + log(n)*df}{BIC = 2*L + log(n)*df}
%%  \deqn{GCV = 2*L/((1-df/n)^2)}{GCV= 2*L/((1-df/n)^2)}
%%  \deqn{AICc = AIC + 2*df*(df+1)/(n-df-1)}{AICc = AIC + 2*df*(df+1)/(n-df-1)}
%%  \deqn{EBIC = BIC + 2*log(p choose df)}{EBIC = BIC + 2*log(p choose df)}
}

\value{
A list containing:
  \item{lambda}{The selected value of the regularization parameter, \code{lambda}.}
  \item{beta}{The vector of coefficients at the chosen value of \code{lambda}.}
  \item{beta.latent}{The vector of latent coefficients at the chosen value of \code{lambda}.}
  \item{df}{The effective number of model parameters at the chosen value of \code{lambda}.}
  \item{IC}{A vector of the calculated model selection criteria for each point on the regularization path.}
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
