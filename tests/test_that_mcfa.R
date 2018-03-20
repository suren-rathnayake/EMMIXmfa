rm(list=ls())
require(testthat)
require(EMMIXmcfa)

set.seed(1)
Y <- iris[, -5]
g <- 3
q <- 2
p <- ncol(Y)
n <- nrow(Y)

context("mcfa")
model <- mcfa(Y, g, q, nkmeans = 2, nrandom = 2, tol = 1.e-5)

expect_that(model, is_a("mcfa"))
expect_that(model, is_a("emmix"))

expect_that(g,   equals(model$g))
expect_that(q,   equals(model$q))

expect_that(g,  equals(length(model$pivec)))
expect_that(1,   equals(sum(model$pivec)))

expect_that(p,   equals(nrow(model$A)))
expect_that(q,   equals(ncol(model$A)))

expect_that(q,   equals(nrow(model$xi)))
expect_that(g,   equals(ncol(model$xi)))

dim_omega <- dim(model$omega)
expect_that(q,   equals(dim_omega[1]))
expect_that(q,   equals(dim_omega[2]))
expect_that(g,   equals(dim_omega[3]))

expect_that(p,   equals(nrow(model$D)))
expect_that(p,   equals(ncol(model$D)))

expect_that(n,   equals(nrow(model$tau)))
expect_that(g,   equals(ncol(model$tau)))

expect_that(n,   equals(nrow(model$Fmat)))
expect_that(q,   equals(ncol(model$Fmat)))

expect_that(n,   equals(nrow(model$UC)))
expect_that(q,   equals(ncol(model$UC)))

dim_U <- dim(model$U)
expect_that(n,   equals(dim_U[1]))
expect_that(q,   equals(dim_U[2]))
expect_that(g,   equals(dim_U[3]))

expect_that(n,   equals(length(model$clust)))


X <- Y
X[1, 4] <- NA
expect_that(mcfa(X, g, q),
            throws_error("`Y' has missing value"))

X[1, 4] <- "a"
expect_that(mcfa(X, g, q),
            throws_error("`Y' has a non-numeric element"))

X <- X[,1]
expect_that(mcfa(X, g, q=1),
            throws_error("The number of factors must be less than the number of variables."))

X <- Y
expect_that(mcfa(X, g, ncol(X)),
            throws_error("The number of factors must be less than the number of variables."))


expect_that(mcfa(X, g, ncol(X)+1),
            throws_error("The number of factors must be less than the number of variables."))

expect_that(mcfa(X, g = g, q = -1),
            throws_error("q must be a positive integer."))


expect_that(mcfa(X, g = -2, q = q),
            throws_error("g must be a positive integer."))


expect_that(mcfa(X, g = 2.5, q = q),
            throws_error("g must be a positive integer."))


expect_that(mcfa(X, g = g, q = 1.5),
            throws_error("q must be a positive integer."))

test_that("mcfa with itmax NULL", {
  expect_that(mcfa (Y, g, q, itmax=NULL), throws_error())
})

test_that("mcfa with itmax NULL", {
  expect_that(mcfa (Y, g, q, itmax=NULL), throws_error())
})

test_that("mcfa with itmax neg", {
  expect_that(mcfa (Y, g, q, itmax=-1), throws_error())
})

test_that("mcfa with init_para as past model", {
  expect_that(mcfa (Y, g, q, init_para = model), is_a("mcfa"))
})


context("q = 1")

q <- 1

context("mcfa")
model <- mcfa(Y, g, q, nkmeans = 2, nrandom = 2, tol = 1.e-5)

expect_that(model, is_a("mcfa"))
expect_that(model, is_a("emmix"))

expect_that(g,   equals(model$g))
expect_that(q,   equals(model$q))

expect_that(g,  equals(length(model$pivec)))
expect_that(1,   equals(sum(model$pivec)))

expect_that(p,   equals(nrow(model$A)))
expect_that(q,   equals(ncol(model$A)))

expect_that(q,   equals(nrow(model$xi)))
expect_that(g,   equals(ncol(model$xi)))

dim_omega <- dim(model$omega)
expect_that(q,   equals(dim_omega[1]))
expect_that(q,   equals(dim_omega[2]))
expect_that(g,   equals(dim_omega[3]))

expect_that(p,   equals(nrow(model$D)))
expect_that(p,   equals(ncol(model$D)))

expect_that(n,   equals(nrow(model$tau)))
expect_that(g,   equals(ncol(model$tau)))

expect_that(n,   equals(nrow(model$Fmat)))
expect_that(q,   equals(ncol(model$Fmat)))

expect_that(n,   equals(nrow(model$UC)))
expect_that(q,   equals(ncol(model$UC)))

dim_U <- dim(model$U)
expect_that(n,   equals(dim_U[1]))
expect_that(q,   equals(dim_U[2]))
expect_that(g,   equals(dim_U[3]))

expect_that(n,   equals(length(model$clust)))


context("missing")

X <- Y
X[1, 4] <- NA
expect_that(mcfa(X, g, q),
            throws_error("`Y' has missing value"))

X[1, 4] <- "a"
expect_that(mcfa(X, g, q),
            throws_error("`Y' has a non-numeric element"))

X <- X[,1]
expect_that(mcfa(X, g, q=1),
            throws_error("The number of factors must be less than the number of variables."))

X <- Y
expect_that(mcfa(X, g, ncol(X)),
            throws_error("The number of factors must be less than the number of variables."))


expect_that(mcfa(X, g, ncol(X)+1),
            throws_error("The number of factors must be less than the number of variables."))

expect_that(mcfa(X, g = g, q = -1),
            throws_error("q must be a positive integer."))


expect_that(mcfa(X, g = -2, q = q),
            throws_error("g must be a positive integer."))


expect_that(mcfa(X, g = 2.5, q = q),
            throws_error("g must be a positive integer."))


expect_that(mcfa(X, g = g, q = 1.5),
            throws_error("q must be a positive integer."))

test_that("mcfa with itmax NULL", {
  expect_that(mcfa (Y, g, q, itmax=NULL), throws_error())
})

test_that("mcfa with itmax NULL", {
  expect_that(mcfa (Y, g, q, itmax=NULL), throws_error())
})

test_that("mcfa with itmax neg", {
  expect_that(mcfa (Y, g, q, itmax=-1), throws_error())
})

a <- mcfa (Y, g, q, init_para = model)
test_that("mcfa with init_para as past model", {
  expect_that(a, is_a("mcfa"))
})
