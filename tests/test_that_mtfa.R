rm(list=ls())
require(testthat)
require(EMMIXmfa)

set.seed(1)
Y <- iris[, -5]
g <- 3
q <- 2
p <- ncol(Y)
n <- nrow(Y)

context("mtfa cov")
model <- mtfa(Y, g, q, nkmeans = 2, nrandom = 2, tol = 1.e-1, 
              sigma_type = "unique", D_type = "unique")
test_that("cov uu", {
  expect_that(model, is_a("emmix"))
})
model <- mtfa(Y, g, q, nkmeans = 2, nrandom = 2, tol = 1.e-1, 
              sigma_type = "unique", D_type = "common")

test_that("cov uc", {
  expect_that(model, is_a("emmix"))
})
context("mtfa")
model <- mtfa(Y, g, q, nkmeans = 2, nrandom = 2, tol = 1.e-1)
mod2 <- mtfa (Y, g, q, init_para = model)
test_that("mtfa with init_para as past model", {
  expect_that(mod2, is_a("mtfa"))
})
test_that("mtfa with init_para as past model", {
  expect_that(mod2, is_a("emmix"))
})
X <- Y
X[1, 4] <- NA
expect_that(mtfa(X, g, q),
            throws_error("`Y' has missing value"))
X[1, 4] <- "a"
expect_that(mtfa(X, g, q),
            throws_error("`Y' has a non-numeric element"))
X <- X[,1]
expect_that(mtfa(X, g, q),
            throws_error("The data must have more than one variable."))
X <- Y
expect_that(mtfa(X, g, ncol(X)),
            throws_error("The number of factors must be less than the number of variables."))
expect_that(mtfa(X, g, ncol(X)+1),
            throws_error("The number of factors must be less than the number of variables."))
expect_that(mtfa(X, g = g, q = -1),
            throws_error("q must be a positive integer."))
expect_that(mtfa(X, g = -2, q = q),
            throws_error("g must be a positive integer."))
expect_that(mtfa(X, g = 2.5, q = q),
            throws_error("g must be a positive integer."))
expect_that(mtfa(X, g = g, q = 1.5),
            throws_error("q must be a positive integer."))
test_that("mtfa with itmax NULL", {
  expect_that(mtfa (Y, g, q, itmax=NULL), throws_error())
})
test_that("mtfa with itmax NULL", {
  expect_that(mtfa (Y, g, q, itmax=NULL), throws_error())
})
test_that("mtfa with itmax neg", {
  expect_that(mtfa (Y, g, q, itmax=-1), throws_error())
})

context("q = 1")
q <- 1
context("mtfa")
model <- mtfa(Y, g, q, nkmeans = 2, nrandom = 2, tol = 1.e-5)

X <- Y
X[1, 4] <- NA
expect_that(mtfa(X, g, q),
            throws_error("`Y' has missing value"))

X[1, 4] <- "a"
expect_that(mtfa(X, g, q),
            throws_error("`Y' has a non-numeric element"))

X <- X[,1]
expect_that(mtfa(X, g, q),
            throws_error("The data must have more than one variable."))

X <- Y
expect_that(mtfa(X, g, ncol(X)),
            throws_error("The number of factors must be less than the number of variables."))

expect_that(mtfa(X, g, ncol(X)+1),
            throws_error("The number of factors must be less than the number of variables."))

expect_that(mtfa(X, g = g, q = -1),
            throws_error("q must be a positive integer."))

expect_that(mtfa(X, g = -2, q = q),
            throws_error("g must be a positive integer."))

expect_that(mtfa(X, g = 2.5, q = q),
            throws_error("g must be a positive integer."))

expect_that(mtfa(X, g = g, q = 1.5),
            throws_error("q must be a positive integer."))

test_that("mtfa with itmax NULL", {
  expect_that(mtfa (Y, g, q, itmax=NULL), throws_error())
})

test_that("mtfa with itmax NULL", {
  expect_that(mtfa (Y, g, q, itmax=NULL), throws_error())
})

test_that("mtfa with itmax neg", {
  expect_that(mtfa (Y, g, q, itmax=-1), throws_error())
})

test_that("mtfa with init_para as past model", {
  expect_that(mtfa (Y, g, q, init_para = model), is_a("mtfa"))
})

dat <- rmix(100, model)

test_that("mtfa with init_para as past model", {
  expect_that(dat, is_a("matrix"))
})

test_that("mtfa with init_para as past model", {
  expect_that(nrow(dat), equals(100))
})

test_that("mtfa with init_para as past model", {
  expect_that(ncol(dat), equals(p))
})
