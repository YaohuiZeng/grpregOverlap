## function: overlapping group selection based on R Package 'grpreg' 
##           by Dr. Patrick Breheny <patrick-breheny@uiowa.edu>
# ------------------------------------------------------------------------------
overlap.grpreg <- function(X, y, group, 
                           penalty=c("grLasso", "grMCP", "grSCAD", "gel", 
                                     "cMCP", "gBridge", "gLasso", "gMCP"), 
                           family=c("gaussian","binomial", "poisson"), 
                           nlambda=100, lambda, 
                           lambda.min={if (nrow(X) > ncol(X)) 1e-4 else .05},
                           alpha=1, eps=.001, max.iter=1000, dfmax=p, gmax=J,
                           gamma=3, tau=1/3, 
                           group.multiplier={if (strtrim(penalty,2)=="gr") 
                             sqrt(table(group[group!=0])) else rep(1,J)}, 
                           warn=TRUE, ...) {
  # Error checking
  if (class(X) != "matrix") {
    tmp <- try(X <- as.matrix(X), silent=TRUE)
    if (class(tmp)[1] == "try-error")  {
      stop("X must be a matrix or able to be coerced to a matrix")
    }   
  }
  if (storage.mode(X)=="integer") X <- 1.0*X

  incid.mat <- incidenceMatrix(X, group) # group membership incidence matrix
  over.mat <- over.temp <- Matrix(incid.mat %*% t(incid.mat)) # overlap matrix
  grp.vec <- rep(1:nrow(over.mat), times = diag(over.mat)) # group index vector
  
  overlap <- TRUE
  diag(over.temp) <- 0
  if (all(over.temp == 0)) {
    overlap <- FALSE
    cat("Note: There are NO overlaps among groups at all! Now conducting non-overlapping group selection!")
    val <- grpreg(X = X, y = y, group = grp.vec, ...)
  } else {
    X.latent <- expandX(X, incid.mat, grp.vec)
    fit <- grpreg(X = X.latent, y = y, group = grp.vec, ...)
    fit$beta.latent <- fit$beta # fit$beta from grpreg is latent beta
    fit$beta <- gamma2beta(gamma=fit$beta, incid.mat, grp.vec)
    fit$incidence.mat <- incid.mat
    fit$overlap.mat <- over.mat
    fit$group <- group
    fit$grp.vec <- grp.vec
    fit$X.latent <- X.latent
    # get results, store in new class 'overlap.grpreg'
    val <- structure(fit,
                     class = c('overlap.grpreg', 'grpreg'))
    val
  }
}
# -------------------------------------------------------------------------------

## function: convert latent beta coefficients (gamma's) to non-latent beta's
# -------------------------------------------------------------------------------
gamma2beta<- function(gamma, incidence.mat, grp.vec) {
  # gamma: matrix, ncol = length(lambda), nrow = # of latent vars.
  p <- ncol(incidence.mat)
  J <- nrow(incidence.mat)
  beta <- matrix(0, ncol = ncol(gamma), nrow = p)
  
  intercept <- gamma[1, , drop = FALSE]
  gamma <- gamma[-1, , drop = FALSE]
  
  for (i in 1:J) {
    ind <- which(incidence.mat[i, ] == 1)
    beta[ind, ] <- beta[ind, ] + gamma[which(grp.vec == i), , drop = FALSE]
  }
  beta <- rbind(intercept, beta)         
  rownames(beta) <- c("(Intercept)", colnames(incidence.mat))
  beta
}
# -------------------------------------------------------------------------------

## function: expand design matrix X to overlapping design matrix (X.latent)
# -------------------------------------------------------------------------------
expandX <- function(X, incidence.mat, grp.vec) {
  # expand X to X.latent
  X.latent <- NULL
  for(i in 1:nrow(incidence.mat)) {
    X.latent <- cbind(X.latent, X[, incidence.mat[i,]==1])
  }
  colnames(X.latent) <- paste('grp', grp.vec, '_', 
                              colnames(X.latent), sep = "")
  X.latent
}
# -------------------------------------------------------------------------------
