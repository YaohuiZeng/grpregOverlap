
## function: plot grpregOverlap
# ------------------------------------------------------------------------------
plot.cv.grpregOverlap <- function(x, log.l=TRUE, type=c("cve", "rsq", "scale", "snr", "pred", "all"), selected=TRUE, vertical.line=TRUE, col="red", ...) {
  if (x$fit$family == 'cox') {
    class(x) <- c('cv.grpsurv', "cv.grpreg")
    plot(x, log.l = log.l, type = type, selected = selected,
         vertical.line = vertical.line, col = col, ...)
  } else {
    plot(x, log.l = log.l, type = type, selected = selected,
         vertical.line = vertical.line, col = col, ...)
  }
} 