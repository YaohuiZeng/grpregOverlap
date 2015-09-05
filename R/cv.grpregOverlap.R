## function: cross-validation
# ------------------------------------------------------------------------------
cv.overlap.grpreg <- function(X, y, group, ..., nfolds=10, seed, trace=FALSE) {
  fit <- overlap.grpreg(X=X, y=y, group=group, ...)
  cvfit <- cv.grpreg(X = fit$X.latent, y = y, group = fit$grp.vec, ...,
                     nfolds = nfolds, seed = seed, 
                     trace = trace)
  cvfit$fit <- fit
  val <- structure(cvfit, class = c('cv.overlap.grpreg', 'cv.grpreg'))
  val
}
# ------------------------------------------------------------------------------
