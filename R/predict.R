## function: predict grpregOverlap
# -------------------------------------------------------------------------------
predict.grpregOverlap <- function(object, X, 
                                   type=c("link", "response", "class", 
                                          "coefficients", "vars", "groups", 
                                          "nvars", "ngroups", "norm"), 
                                   latent = FALSE, lambda, 
                                   which=1:length(object$lambda), ...) {
  family <- object$family
  if (!missing(X) && is.character(X)) {
    type <- X
    X <- NULL
  }
  type <- match.arg(type)
  beta <- coef(object=object, lambda=lambda, latent = latent, 
               which=which, drop = FALSE, ...)
  if (type == 'coefficients') return(beta)
  if (length(dim(object$beta)) == 2) {
    if (type=="vars") {
      if (family != 'cox') {
        return(drop(apply(beta[-1, , drop=FALSE]!=0, 2, FUN=which)))
      } else {
        return(drop(apply(beta != 0, 2, FUN=which)))
      }
    }
    if (type=="nvars") {
      if (family != 'cox') {
        v <- drop(apply(beta[-1, , drop=FALSE]!=0, 2, FUN=which))
      } else {
        v <- drop(apply(beta != 0, 2, FUN=which))
      }
      if (is.list(v)) {
        res <- sapply(v, length)
      } else {
        res <- length(v)
      }
      return(res)
    }
    if (type == 'norm') {
      beta <- coef(object=object, lambda=lambda, latent = TRUE, 
                   which=which, drop = FALSE, ...)
      if (!latent) {
        cat("The returned is the L2 norm of the latent coefficients! Set latent = 'TRUE' to avoid this warning message.\n\n")
      } 
      if (family != 'cox') {
        return(drop(apply(beta[-1, , drop=FALSE], 2, 
                          function(x) tapply(x, object$grp.vec, function(x){sqrt(sum(x^2))}))))
      } else {
        return(drop(apply(beta, 2, function(x) tapply(x, object$grp.vec, function(x){sqrt(sum(x^2))}))))
      }
    }
  } else {
    if (type=="vars") 
      stop("Predicting type 'vars' not implemented with multivariate outcomes\n")
    if (type=="nvars") {
      return(drop(apply(beta[,-1, , drop=FALSE]!=0, 3, FUN=sum)))
    }
    if (type == 'norm') {
      beta <- coef(object=object, lambda=lambda, latent = TRUE, 
                   which=which, drop = FALSE, ...)
      if (!latent) {
        cat("The returned is the L2 norm of the latent coefficients. Set latent = 'TRUE' to avoid this warning message.\n\n")
      } 
      if (family != 'cox') {
        return(drop(apply(beta[, -1, , drop=FALSE], 3, function(x) apply(x, 2, function(x){sqrt(sum(x^2))})))) 
      } else {
        return(drop(apply(beta, 3, function(x) apply(x, 2, function(x){sqrt(sum(x^2))})))) 
      }
    }
  }
  if (!missing(X) && !is.null(X)) {
    X <- expandX(X, object$group)
  }
  if (latent) {
    cat("Only latent 'coefficients', 'vars', 'nvars', 'norm' can be returned! Set latent = 'FALSE' to suppress this message.\n\n")
  }
  obj.new <- object
  obj.new$group <- object$grp.vec
  obj.new$beta <- object$beta.latent
  if (family != 'cox') {
    class(obj.new) <- 'grpreg'
  } else {
    class(obj.new) <- 'grpsurv'
  }
  return(predict(obj.new, X=X, type=type, lambda=lambda, which=which, ...))
}
# -------------------------------------------------------------------------------

## function: coef.grpregOverlap, coef for grpregOverlap
# -------------------------------------------------------------------------------
coef.grpregOverlap <- function(object, lambda, latent = FALSE, 
                               which=1:length(object$lambda), drop=TRUE, ...) {
  family <- object$family
  obj.new <- object
  obj.new$beta <- object$beta.latent
  class(obj.new) <- 'grpreg'

  ## latent beta
  beta <- coef(object = obj.new, lambda = lambda, 
               which = which, drop=FALSE, ...)  
  if (!latent) {
    beta <- gamma2beta(beta, incidence.mat = object$incidence.mat, 
                       grp.vec = object$grp.vec, family = family)
  }
  if (drop) return(drop(beta)) else return(beta)
}
# -------------------------------------------------------------------------------
