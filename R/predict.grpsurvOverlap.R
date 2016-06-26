## function: predict grpsurvOverlap
# -------------------------------------------------------------------------------
predict.grpsurvOverlap <- function(object, X, 
                                   type=c("link", "response", "survival", "median",
                                          "norm", "coefficients", "vars", "nvars",
                                          "groups", "ngroups"),
                                   latent = FALSE, lambda, 
                                   which=1:length(object$lambda), ...) {
  
  type <- match.arg(type)
  if (type %in% c("norm", "coefficients", "vars", "nvars", "groups", "ngroups")) {
    class(object) <- 'grpregOverlap'
    return(predict(object=object, X=X, type=type, latent=latent,
                   lambda=lambda, which=which, ...))
  } else {
    class(object) <- 'grpsurv'
    return(predict(object = object, X = X, type = type,
                   lambda = lambda, which = which))
  }
}
# -------------------------------------------------------------------------------
