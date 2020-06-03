
## function: overlap matrix c[i,j] = # of overlaps between group i and j.
##           Diagonals are group size. Provide interface for user to check 
##           overlapping structure.
# ------------------------------------------------------------------------------
overlapMatrix <- function(X, group) {
  inc.mat <- incidenceMatrix(X, group)
  over.mat <- Matrix(inc.mat %*% t(inc.mat), sparse = TRUE, dimnames = dimnames(inc.mat))
  over.mat
}
# ------------------------------------------------------------------------------

## function: incidence matrix: I[i, j] = 1 if group i contains variable j.
# ------------------------------------------------------------------------------
incidenceMatrix <- function(X, group) {
  n <- nrow(X)
  p <- ncol(X)
  if (! is.list(group)) {
    stop("Argument 'group' must be a list of integer indices or character names of variables!")
  }
  J <- length(group)
  grp.mat <- Matrix(0, nrow = J, ncol = p, sparse = TRUE, 
                    dimnames=list(as.character(rep(NA, J)),
                                  as.character(rep(NA, p))))    
  if(is.null(colnames(X))) {
    colnames(X) <- paste("V", 1:ncol(X), sep="")    
  }
  if (is.null(names(group))) {
    names(group) <- paste("grp", 1:J, sep="")
  }
  
  if (class(group[[1]]) == 'numeric') {
    for (i in 1:J) {
      ind <- group[[i]]
      grp.mat[i, ind] <- 1
      colnames(grp.mat)[ind] <- colnames(X)[ind]
    }
  } else { ## character, names of variables
    for (i in 1:J) {
      grp.i <- as.character(group[[i]])
      ind <- colnames(X) %in% grp.i
      grp.mat[i, ] <- 1*ind
      colnames(grp.mat)[ind] <- colnames(X)[ind]
    }
  }
  rownames(grp.mat) <- as.character(names(group))
  # check grp.mat
  if (all(grp.mat == 0)) {
    stop("The names of variables in X don't match with names in group!")
  }

  grp.mat
}
# ------------------------------------------------------------------------------
