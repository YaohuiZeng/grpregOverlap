## function: cross-validation
# ------------------------------------------------------------------------------
cv.grpregOverlap <- function(X, y, group, ...) {
  fit <- grpregOverlap(X=X, y=y, group=group, returnX = TRUE, ...)
  cvfit <- cv.grpreg(X = fit$X.latent, y = y, group = fit$grp.vec, ...)
  cvfit$fit <- fit
  val <- structure(cvfit, class = c('cv.grpregOverlap', 'cv.grpreg'))
  val
}
# ------------------------------------------------------------------------------
