## function: cross-validation for survival
# ------------------------------------------------------------------------------
cv.grpsurvOverlap <- function(X, y, group, ...) {
  fit <- grpregOverlap(X=X, y=y, group=group, returnX = TRUE, family = 'cox', ...)
  cvfit <- cv.grpsurv(X = fit$X.latent, y = y, group = fit$grp.vec, ...)
  cvfit$fit <- fit
  val <- structure(cvfit, class = c('cv.grpregOverlap', 'cv.grpsurv', 'cv.grpreg'))
  val
}
# ------------------------------------------------------------------------------
