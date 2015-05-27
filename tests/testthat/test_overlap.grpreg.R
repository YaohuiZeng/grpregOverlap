library(testthat)
library(overlap.grpreg)
library(grpreg)
# test_file("test_incidence.R")

context("Testing overlap.grpreg()")
data(birthwt.grpreg)
data(pathway.dat)

X <- pathway.dat$expression
group <- pathway.dat$pathways
# Linear regression
fit <- overlap.grpreg(X, y, group, penalty = 'grLasso')
fit <- overlap.grpreg(X, y, group, penalty = 'grMCP')
fit <- overlap.grpreg(X, y, group, penalty = 'grSCAD')
fit <- overlap.grpreg(X, y, group, penalty = 'gel')
fit <- overlap.grpreg(X, y, group, penalty = 'cMCP')
str(fit)
plot(fit)
plot(fit, latent = FALSE)
plot(fit, norm = T)
select(fit, "AIC")

test_that("group argument, a list?", {
  X <- as.matrix(birthwt.grpreg[,-1:-2])
  group <- c(1,1,1,2,2,2,3,3,4,5,5,6,7,8,8,8)
  expect_that(incidenceMatrix(X, group), 
              throws_error("Argument 'group' must be a list"))
  
  group <- list(c(1, 2, 3), c(4, 5, 6), c(7, 8), c(9), c(10, 11), 
                c(12), c(13), c(14, 15, 16))
  expect_that(incidenceMatrix(X, group), 
              not(throws_error("Argument 'group' must be a list")))
  
  group <- lapply(group, function(x) colnames(X)[x])
  expect_that(incidenceMatrix(X, group), 
              not(throws_error("Argument 'group' must be a list")))
  
  X <- pathway.dat$expression
  group <- pathway.dat$pathways
  expect_that(incidenceMatrix(X, group), 
              not(throws_error("Argument 'group' must be a list")))  
})

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





