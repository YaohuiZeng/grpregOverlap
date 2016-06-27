\name{plot.cv.grpregOverlap}
\alias{plot.cv.grpregOverlap}

\title{Plots the cross-validation curve from cross-validated object}
\description{Plots the cross-validation curve from a a \code{cv.grpregOverlap} or \code{cv.grpsurvOverlap} object.}

\usage{
\method{plot}{cv.grpregOverlap}(x, log.l=TRUE, type=c("cve", "rsq", "scale",
"snr", "pred", "all"),  selected=TRUE, vertical.line=TRUE, col="red",
...)}
\arguments{
  \item{x}{A \code{cv.grpregOverlap} or \code{cv.grpsurvOverlap} object.}
  \item{log.l}{Should horizontal axis be on the log scale?  Default is
    TRUE.}
  \item{type}{What to plot on the vertical axis.  \code{cve} plots the
    cross-validation error (deviance); \code{rsq} plots an estimate of
    the fraction of the deviance explained by the model (R-squared);
    \code{snr} plots an estimate of the signal-to-noise ratio;
    \code{scale} plots, for \code{family="gaussian"}, an estimate of the
    scale parameter (standard deviation); \code{pred} plots, for
    \code{family="binomial"}, the estimated prediction error; \code{all}
    produces all of the above.}
  \item{selected}{If \code{TRUE} (the default), places an axis on top of
    the plot denoting the number of groups in the model (i.e., that
    contain a nonzero regression coefficient) at that value of
    \code{lambda}.}  
  \item{vertical.line}{If \code{TRUE} (the default), draws a vertical
    line at the value where cross-validaton error is minimized.}
  \item{col}{Controls the color of the dots (CV estimates).}
  \item{\dots}{Other graphical parameters to \code{plot}}
}

\author{
  Yaohui Zeng and Patrick Breheny
  
  Maintainer: Yaohui Zeng <yaohui-zeng@uiowa.edu>
}

