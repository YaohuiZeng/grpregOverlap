library(testthat)
library(overlap.grpreg)
library(grpreg)
# test_file("test_incidence.R")

context("Testing overlap.grpreg()")

test_that("Non-overlapping fit againt grpreg:", {
  data(birthwt.grpreg)
  X <- as.matrix(birthwt.grpreg[,-1:-2])
  ## no overlap, should be the same as from 'grpreg'
  group <- list(c(1, 2, 3), c(4, 5, 6), c(7, 8), c(9), c(10, 11), 
                c(12), c(13), c(14, 15, 16))
  group2 <- rep(1:length(group), sapply(group, length))
  
  ## linear regression
  y <- birthwt.grpreg$bwt
  invisible(capture.output({
    fit <- overlap.grpreg(X, y, group, penalty = 'grLasso')
  }))
  # x == x.latent, test expandX()
  expect_equal(all(X == fit$X.latent), TRUE) 
  # beta = beta.latent, test gamma2beta()
  # Pass test ONLY IF variable indices in ascending order, otherwise, 
  # variables are reordered in beta.latent according to group order.
  expect_equal(all(fit$beta == fit$beta.latent), TRUE)
  
  # equivalent to fit 'grpreg'
  fit2 <- grpreg(X, y, group2, penalty = 'grLasso')
  expect_identical(fit$beta, fit2$beta)
  expect_equal(all(fit2$beta == fit$beta.latent), TRUE)
  
  ## logistic regression
  y <- birthwt.grpreg$low
  invisible(capture.output({
    fit <- overlap.grpreg(X, y, group, penalty = 'grLasso', family = 'binomial')
  }))
  # x == x.latent, test expandX()
  expect_equal(all(X == fit$X.latent), TRUE) 
  # beta = beta.latent, test gamma2beta()
  # Pass test ONLY IF variable indices in ascending order, otherwise, 
  # variables are reordered in beta.latent according to group order.
  expect_equal(all(fit$beta == fit$beta.latent), TRUE)
  
  # equivalent to fit 'grpreg'
  fit2 <- grpreg(X, y, group2, penalty = 'grLasso', family = 'binomial')
  expect_identical(fit$beta, fit2$beta)
  expect_equal(all(fit2$beta == fit$beta.latent), TRUE)
  
  ## logistic regression
  y <- birthwt.grpreg$low
  invisible(capture.output({
    fit <- overlap.grpreg(X, y, group, penalty = 'grLasso', family = 'binomial')
  }))
  # x == x.latent, test expandX()
  expect_equal(all(X == fit$X.latent), TRUE) 
  # beta = beta.latent, test gamma2beta()
  # Pass test ONLY IF variable indices in ascending order, otherwise, 
  # variables are reordered in beta.latent according to group order.
  expect_equal(all(fit$beta == fit$beta.latent), TRUE)
  
  # equivalent to fit 'grpreg'
  fit2 <- grpreg(X, y, group2, penalty = 'grLasso', family = 'binomial')
  expect_identical(fit$beta, fit2$beta)
  expect_equal(all(fit2$beta == fit$beta.latent), TRUE)
})


test_that("predict, coef, select, cv, against grpreg: ", {
  data(birthwt.grpreg)
  X <- as.matrix(birthwt.grpreg[,-1:-2])
  ## no overlap, should be the same as from 'grpreg'
  group <- list(c(1, 2, 3), c(4, 5, 6), c(7, 8), c(9), c(10, 11), 
                c(12), c(13), c(14, 15, 16))
  group2 <- rep(1:length(group), sapply(group, length))
  ## logistic regression
  y <- birthwt.grpreg$low
  invisible(capture.output({
    fit <- overlap.grpreg(X, y, group, penalty = 'grMCP', family = 'binomial')
  }))
  
  ## test predict, coef againt grpreg
  fit2 <- grpreg(X, y, group2, penalty = 'grMCP', family = 'binomial')
  
  expect_equal(all(coef(fit, lambda=.001) == coef(fit2, lambda=.001)), TRUE)
  expect_equal(all(predict(fit, X, type="link", lambda=.001) == 
                     predict(fit2, X, type="link", lambda=.001)
    ), TRUE)
  expect_equal(all(predict(fit, X, type="response", lambda=.001) ==
                     predict(fit2, X, type="response", lambda=.001)
    ), TRUE)
  expect_equal(all(
    ), TRUE)
  
  expect_equal(all(predict(fit, X, type="class", lambda=.001) ==
                     predict(fit2, X, type="class", lambda=.001)
    ), TRUE)
  expect_equal(all(predict(fit, type="vars", lambda=.07) ==
                     predict(fit, type='vars', latent=T, lambda=.07)
    ), TRUE)
  expect_equal(all(predict(fit, type="nvars", lambda=.07) ==
                     predict(fit, type='nvars', latent=T, lambda=.07)
  ), TRUE)
  expect_equal(all(predict(fit, type="coefficients", lambda=.07) ==
                     predict(fit, type="coefficients", latent = T, lambda=.07)
    ), TRUE)

  expect_equal(all(predict(fit, type="vars", lambda=.07) ==
                     predict(fit2, type="vars", lambda=.07)
    ), TRUE)
  
  expect_equal(all(predict(fit, type="groups", lambda=.07) ==
                     predict(fit2, type="groups", lambda=.07)
    ), TRUE)
  
  invisible(capture.output({
    expect_equal(all(predict(fit, type="norm", lambda=.07) ==
                       predict(fit2, type="norm", lambda=.07)
    ), TRUE)
    }))
  
  invisible(capture.output({
    cvfit <- cv.overlap.grpreg(X, y, group, family="binomial", penalty="grMCP", 
                               seed = 1234)
  }))
  
  cvfit2 <- cv.grpreg(X, y, group2, family="binomial", penalty="grMCP", 
                      seed = 1234)
  expect_equal(all(coef(cvfit) == coef(cvfit2)), TRUE)
  expect_equal(all(predict(cvfit, X) == predict(cvfit2, X)), TRUE)
  expect_equal(all(predict(cvfit, X, type="response") ==
                     predict(cvfit2, X, type="response")
    ), TRUE)
  expect_equal(all(predict(cvfit, type="groups") == 
                     predict(cvfit2, type="groups")
    ), TRUE)
  
  ## test select
  y <- birthwt.grpreg$bwt
  invisible(capture.output({
    fit <- overlap.grpreg(X, y, group, penalty="grLasso")
  }))
  fit2 <- grpreg(X, y, group2, penalty="grLasso")
  expect_equal(all(select(fit)$lambda == select(fit2)$lambda), TRUE)
  expect_equal(all(select(fit)$beta == select(fit2)$beta), TRUE)
  suppressWarnings(
    expect_equal(all(select(fit,crit="AIC",df="active")$beta.latent ==
                     select(fit2,crit="AIC",df="active")$beta
     ),TRUE)  
  )  
})

# test_that("Overlapping fit: ", {
#   data(pathway.dat)
#   ## logistic regression, pathway selection
#   X <- pathway.dat$expression
#   group <- pathway.dat$pathways
#   y <- pathway.dat$mutation
#   fit <- overlap.grpreg(X, y, group, penalty = 'grLasso', family = 'binomial')
# #   fit <- overlap.grpreg(X, y, group, penalty = 'grMCP', family = 'binomial')
# #   fit <- overlap.grpreg(X, y, group, penalty = 'grSCAD', family = 'binomial')
# #   fit <- overlap.grpreg(X, y, group, penalty = 'gel', family = 'binomial')
# #   fit <- overlap.grpreg(X, y, group, penalty = 'cMCP', family = 'binomial')
# #   print(object.size(cvfit), units = 'Mb')
#   plot(fit)
#   plot(fit, latent = FALSE)
#   predict(fit, type = 'ngroups')
#   plot(fit, norm = T)
#   str(select(fit, "AIC"))
#   
#   predict(fit, X, type="class", lambda=0.04)
# 
#   cvfit <- cv.overlap.grpreg(X, y, group, penalty = 'grLasso', family = 'binomial')
#   coef(cvfit)
#   predict(cvfit, X)
#   predict(cvfit, X, type='response')
#   predict(cvfit, X, type = 'link')
#   predict(cvfit, X, type = 'class')
#   plot(cvfit)
#   plot(cvfit, type = 'all')
#   summary(cvfit)
# })













test_that('group information, correct?', {
  X <- as.matrix(birthwt.grpreg[,-1:-2])
  group <- list(c('a', 'b', 'c'), c('d', 'e', 'f'), c('g', 'h'), c('i'), 
                c('j', 'k'), c('l'), c('m'), c('n', 'o', 'p'))
  expect_that(incidenceMatrix(X, group), 
              throws_error("The names of variables in X don't match with names in group"))
  
  group <- list(c(1, 2, 3), c(4, 5, 6), c(7, 8), c(10, 11), c(14, 15, 16))
  ## create new groups for variables 9, 12, 13, put them at bottom
  expect_equal(rownames(incidenceMatrix(X, group)), 
               paste(paste("grp", 1:8, sep="")))
  
  colnames(X) <- letters[1:ncol(X)]
  group <- list(c('a', 'b', 'c'), c('d', 'e', 'f'), c('g', 'h'),  
                c('j', 'k'), c('n', 'o', 'p'))
  expect_equal(rownames(incidenceMatrix(X, group)), 
               paste(paste("grp", 1:8, sep="")))
})

test_that("incidence matrix, correct?", {
  X <- as.matrix(birthwt.grpreg[,-1:-2])
  group <- list(c(1, 2, 3), c(4, 5, 6), c(7, 8), c(9), c(10, 11), 
                c(12), c(13), c(14, 15, 16))
  expect_is(incidenceMatrix(X, group), 'dgCMatrix')
  
  X <- pathway.dat$expression
  group <- pathway.dat$pathways
  incid.mat <- incidenceMatrix(X, group)
  expect_is(incid.mat, 'dgCMatrix')
  expect_equal(dim(incid.mat), c(308, 308))
})





