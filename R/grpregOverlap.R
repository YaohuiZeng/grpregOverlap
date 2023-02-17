## function: overlapping group selection based on R Package 'grpreg' 
## update (6/21/2016): adapt for cox model
# ------------------------------------------------------------------------------
grpregOverlap <- function(X, y, group, 
                          family=c("gaussian","binomial", "poisson", 'cox'), 
                          returnX.latent = FALSE,
                          returnOverlap = FALSE,
                          ...) {

  # Error checking
  if (is.matrix(X)) {
    tmp <- try(X <- as.matrix(X), silent=TRUE)
    if (class(tmp)[1] == "try-error")  {
      stop("X must be a matrix or able to be coerced to a matrix")
    }   
  }
  if (storage.mode(X)=="integer") X <- 1.0*X

  incid.mat <- incidenceMatrix(X, group) # group membership incidence matrix
  over.mat <- over.temp <- Matrix(incid.mat %*% t(incid.mat)) # overlap matrix
  grp.vec <- rep(1:nrow(over.mat), times = diag(over.mat)) # group index vector
  X.latent <- expandX(X, group)
 
  diag(over.temp) <- 0
  if (all(over.temp == 0)) {
    cat("Note: There are NO overlaps between groups at all!", "\n") 
    cat("      Now conducting non-overlapping group selection ...")
  }
  
  family <- match.arg(family)
  if (family != 'cox') {
    fit <- grpreg(X = X.latent, y = y, group = grp.vec, family = family, ...)
  } else {
    ## survival analysis
    fit <- grpsurv(X = X.latent, y = y, group = grp.vec, ...)
  }
 
  fit$beta.latent <- fit$beta # fit$beta from grpreg is latent beta
  fit$beta <- gamma2beta(gamma = fit$beta, incid.mat, grp.vec, family = family)
  fit$incidence.mat <- incid.mat
  fit$group <- group
  fit$grp.vec <- grp.vec # this is 'group' argument in Package 'grpreg'
  fit$family <- family
  if (returnX.latent) {
    fit$X.latent <- X.latent
  } 
  if (returnOverlap) {
    fit$overlap.mat <- over.mat
  }
  
  if (family != 'cox') {
    # get results, store in new class 'grpregOverlap', and inherited from 'grpreg'
    val <- structure(fit,
                     class = c('grpregOverlap', 'grpreg'))
  } else {
    val <- structure(fit, 
                     class = c("grpsurvOverlap", "grpregOverlap"))
  }
  val
}
# -------------------------------------------------------------------------------

## function: convert latent beta coefficients (gamma's) to non-latent beta's
## update (6/21/2016): adapt for cox model
# -------------------------------------------------------------------------------
gamma2beta<- function(gamma, incidence.mat, grp.vec, family) {
  # gamma: matrix, ncol = length(lambda), nrow = # of latent vars.
  p <- ncol(incidence.mat)
  J <- nrow(incidence.mat)
  beta <- matrix(0, ncol = ncol(gamma), nrow = p)
  
  if (family != 'cox') {
    intercept <- gamma[1, , drop = FALSE]
    gamma <- gamma[-1, , drop = FALSE]
  } else {
    # Cox model doesn't have an intercept
    gamma <- gamma
  }
  
  for (i in 1:J) {
    ind <- which(incidence.mat[i, ] == 1)
    beta[ind, ] <- beta[ind, ] + gamma[which(grp.vec == i), , drop = FALSE]
  }
  if (family != 'cox') {
    beta <- rbind(intercept, beta)         
    rownames(beta) <- c("(Intercept)", colnames(incidence.mat))
  } else {
    rownames(beta) <- colnames(incidence.mat)
  }
  beta
}
# -------------------------------------------------------------------------------


## function: expand design matrix X to overlapping design matrix (X.latent)
# -------------------------------------------------------------------------------
expandX <- function(X, group) {
  incidence.mat <- incidenceMatrix(X, group) # group membership incidence matrix
  over.mat <- Matrix(incidence.mat %*% t(incidence.mat), sparse = TRUE) 
                     #dimnames = dimnames(incidence.mat)) # overlap matrix
  grp.vec <- rep(1:nrow(over.mat), times = diag(over.mat)) # group index vector
  
  # expand X to X.latent
  X.latent <- NULL
  names <- NULL

  ## the following code will automatically remove variables not included in 'group'
  for(i in 1:nrow(incidence.mat)) {
    idx <- incidence.mat[i,]==1
    X.latent <- cbind(X.latent, X[, idx, drop=FALSE])
    names <- c(names, colnames(incidence.mat)[idx])
#     colnames(X.latent) <- c(colnames(X.latent), colnames(X)[incidence.mat[i,]==1])
  }
  colnames(X.latent) <- paste('grp', grp.vec, '_', names, sep = "")
  X.latent
}
# -------------------------------------------------------------------------------
