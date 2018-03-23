rm(list=ls())
require(testthat)
require(EMMIXmfa)

set.seed(1)
Y <- iris[, -5]
g <- 3
q <- 2
p <- ncol(Y)
n <- nrow(Y)
context("mfa")
model <- mfa(Y, g, q, nkmeans = 2, nrandom = 2, tol = 1.e-5, 
             sigma_type = "unique", D_type = "common")
expect_that(model, is_a("mfa"))
expect_that(model, is_a("emmix"))
expect_that(g,   equals(model$g))
expect_that(q,   equals(model$q))
expect_that(g,  equals(length(model$pivec)))
expect_that(1,   equals(sum(model$pivec)))
expect_that(p,   equals(nrow(model$mu)))
expect_that(g,   equals(ncol(model$mu)))
dim_D <- dim(model$D)
if (model$D_type == "unique") {
  expect_that(p,   equals(dim_D[1]))
  expect_that(p,   equals(dim_D[2]))
  expect_that(g,   equals(dim_D[3]))
} else {
  expect_that(p,   equals(dim_D[1]))
  expect_that(p,   equals(dim_D[2]))  
}
dim_sigma <- dim(model$B)
if (model$sigma_type == "unique") {
  expect_that(p,   equals(dim_sigma[1]))
  expect_that(q,   equals(dim_sigma[2]))
  expect_that(g,   equals(dim_sigma[3]))
} else {
  expect_that(p,   equals(dim_sigma[1]))
  expect_that(q,   equals(dim_sigma[2]))  
}

expect_that(n, equals(nrow(model$tau)))
expect_that(g, equals(ncol(model$tau)))
expect_that(n, equals(nrow(model$Fmat)))
expect_that(q, equals(ncol(model$Fmat)))
expect_that(n, equals(nrow(model$UC)))
expect_that(q, equals(ncol(model$UC)))
dim_U <- dim(model$U)
expect_that(n, equals(dim_U[1]))
expect_that(q, equals(dim_U[2]))
expect_that(g, equals(dim_U[3]))
expect_that(n, equals(length(model$clust)))

model <- mfa(Y, g, q, nkmeans = 2, nrandom = 2, tol = 1.e-5,
                          sigma_type = "unique", D_type = "unique")
expect_that(model, is_a("mfa"))
expect_that(model, is_a("emmix"))
expect_that(g,   equals(model$g))
expect_that(q,   equals(model$q))
expect_that(g,  equals(length(model$pivec)))
expect_that(1,   equals(sum(model$pivec)))
expect_that(p,   equals(nrow(model$mu)))
expect_that(g,   equals(ncol(model$mu)))
dim_D <- dim(model$D)
if (model$D_type == "unique") {
    expect_that(p,   equals(dim_D[1]))
  expect_that(p,   equals(dim_D[2]))
    expect_that(g,   equals(dim_D[3]))
} else {
    expect_that(p,   equals(dim_D[1]))
  expect_that(p,   equals(dim_D[2]))
}
dim_sigma <- dim(model$B)
if (model$sigma_type == "unique") {
    expect_that(p,   equals(dim_sigma[1]))
  expect_that(q,   equals(dim_sigma[2]))
    expect_that(g,   equals(dim_sigma[3]))
} else {
    expect_that(p,   equals(dim_sigma[1]))
  expect_that(q,   equals(dim_sigma[2]))
}

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




test_that("mfa with init_para as past model", {
  expect_that(model, is_a("mfa"))
})

test_that("mfa with init_para as past model", {
  expect_that(model, is_a("emmix"))
})

X <- Y
X[1, 4] <- NA
expect_that(mfa(X, g, q),
            throws_error("`Y' has missing value"))
X[1, 4] <- "a"
expect_that(mfa(X, g, q),
            throws_error("`Y' has a non-numeric element"))
X <- X[,1]
expect_that(mfa(X, g, q=1),
            throws_error("The data must have more than one variable."))
X <- Y
expect_that(mfa(X, g, ncol(X)),
            throws_error("The number of factors must be less than the number of variables."))
expect_that(mfa(X, g, ncol(X)+1),
            throws_error("The number of factors must be less than the number of variables."))
expect_that(mfa(X, g = g, q = -1),
            throws_error("q must be a positive integer."))
expect_that(mfa(X, g = -2, q = q),
            throws_error("g must be a positive integer."))
expect_that(mfa(X, g = 2.5, q = q),
            throws_error("g must be a positive integer."))
expect_that(mfa(X, g = g, q = 1.5),
            throws_error("q must be a positive integer."))
test_that("mfa with itmax NULL", {
  expect_that(mfa (Y, g, q, itmax=NULL), throws_error())
})
test_that("mfa with itmax NULL", {
  expect_that(mfa (Y, g, q, itmax=NULL), throws_error())
})
test_that("mfa with itmax neg", {
  expect_that(mfa (Y, g, q, itmax=-1), throws_error())
})

library(mvtnorm)
p1 <- 10
p2 <- 100
p = p1+p2
q = 2
g = 5
pi1=0.15; pi2=0.2; pi3=0.15; pi4=0.2; pi5=0.3
A1 <- rbind( c(0.5, 0),   c(-0.9, 0),    c(0.3, 0),
             c(0.6, 0.8), c(0.2, -0.7), c(-0.7, 0.5),
             c(0, 0.6), c(0, -0.4), c(0, 0.3), c(0, -0.5))
A2 <- matrix(rep(0, p2*q), c(p2, q))
xi1 <- c( 0,   2.5)
xi2 <- c(-2.5, 0)
xi3 <- c( 2.5, 0)
xi4 <- c( 0,  -2.5)
xi5 <- c( 0,   0)
Om1 <- rbind( c(0.1, 0), c(0, 0.45))
Om2 <- rbind( c(0.45,0), c(0, 0.1))
Om3 <- rbind( c(0.45,0), c(0, 0.1))
Om4 <- rbind( c(0.1, 0), c(0, 0.45))
Om5 <- rbind( c(1, 0.9), c(0.9, 1)) 	
A <- rbind(A1, A2)
D <- diag(c(runif(p1, 0.1, 0.3), runif(p2, 0.3, 0.8)))
pivec <- c(pi1, pi2, pi3, pi4, pi5)

mu    <- array(NA, c(p, g))
mu <- cbind(A%*%xi1, A%*%xi2, A%*%xi3, A%*%xi4, A%*%xi5)

sigma <- array(NA, c(p, p, g))
sigma[,,1] <- A %*% Om1 %*% t(A) + D
sigma[,,2] <- A %*% Om2 %*% t(A) + D
sigma[,,3] <- A %*% Om3 %*% t(A) + D
sigma[,,4] <- A %*% Om4 %*% t(A) + D
sigma[,,5] <- A %*% Om5 %*% t(A) + D

Y <- array(NA, c(n, p))
cls <- array(NA, c(n, 1))

for( i in 1:n) 
{
  Pi <- 0
  r <- runif(1)
  for( j in 1:g)
  {
    Pi <- Pi + pivec[j]
    if( r < Pi )
    {
      comp <- j
      cls[i] <- j
      break
    }
  }
  Y[i,] <-  rmvnorm(1, mu[,comp], sigma[,,comp])
}


context("q = 1")
q <- 1
g <- 2
p <- ncol(Y)
n <- nrow(Y)

context("mfa")
model <- mfa(Y, g, q, nkmeans = 2, nrandom = 2, tol = 1.e-5,
        sigma_type = "unique", D_type = "unique")

expect_that(model, is_a("mfa"))
expect_that(model, is_a("emmix"))
expect_that(g,   equals(model$g))
expect_that(q,   equals(model$q))
expect_that(g,  equals(length(model$pivec)))
expect_that(1,   equals(sum(model$pivec)))
expect_that(p,   equals(nrow(model$mu)))
expect_that(g,   equals(ncol(model$mu)))

dim_D <- dim(model$D)
if (model$D_type == "unique") {
  expect_that(p,   equals(dim_D[1]))
  expect_that(p,   equals(dim_D[2]))
  expect_that(g,   equals(dim_D[3]))
} else {
  expect_that(p,   equals(dim_D[1]))
  expect_that(p,   equals(dim_D[2]))  
}
dim_sigma <- dim(model$B)
if (model$sigma_type == "unique") {
  expect_that(p,   equals(dim_sigma[1]))
  expect_that(q,   equals(dim_sigma[2]))
  expect_that(g,   equals(dim_sigma[3]))
} else {
  expect_that(p,   equals(dim_sigma[1]))
  expect_that(q,   equals(dim_sigma[2]))  
}
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

context("q = 2")
q <- 2
context("mfa")
model <- mfa(Y, g, q, nkmeans = 2, nrandom = 2, tol = 1.e-5,
             sigma_type = "unique", D_type = "unique")
expect_that(model, is_a("mfa"))
expect_that(model, is_a("emmix"))
expect_that(g,   equals(model$g))
expect_that(q,   equals(model$q))
expect_that(g,  equals(length(model$pivec)))
expect_that(1,   equals(sum(model$pivec)))


expect_that(p,   equals(nrow(model$mu)))
expect_that(g,   equals(ncol(model$mu)))

dim_D <- dim(model$D)
if (model$D_type == "unique") {
  expect_that(p,   equals(dim_D[1]))
  expect_that(p,   equals(dim_D[2]))
  expect_that(g,   equals(dim_D[3]))
} else {
  expect_that(p,   equals(dim_D[1]))
  expect_that(p,   equals(dim_D[2]))  
}
dim_sigma <- dim(model$B)
if (model$sigma_type == "unique") {
  expect_that(p,   equals(dim_sigma[1]))
  expect_that(q,   equals(dim_sigma[2]))
  expect_that(g,   equals(dim_sigma[3]))
} else {
  expect_that(p,   equals(dim_sigma[1]))
  expect_that(q,   equals(dim_sigma[2]))  
}
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
expect_that(mfa(X, g, q),
            throws_error("`Y' has missing value"))
X[1, 4] <- "a"
expect_that(mfa(X, g, q),
            throws_error("`Y' has a non-numeric element"))
X <- Y[,1]
expect_that(mfa(X, g, q=p),
            throws_error("The data must have more than one variable."))
X <- Y
expect_that(mfa(X, g, ncol(X)),
            throws_error("The number of factors must be less than the number of variables."))
expect_that(mfa(X, g, ncol(X)+1),
            throws_error("The number of factors must be less than the number of variables."))
expect_that(mfa(X, g = g, q = -1),
            throws_error("q must be a positive integer."))
expect_that(mfa(X, g = -2, q = q),
            throws_error("g must be a positive integer."))
expect_that(mfa(X, g = 2.5, q = q),
            throws_error("g must be a positive integer."))
expect_that(mfa(X, g = g, q = 1.5),
            throws_error("q must be a positive integer."))
test_that("mfa with itmax NULL", {
  expect_that(mfa (Y, g, q, itmax=NULL), throws_error())
})
test_that("mfa with itmax NULL", {
  expect_that(mfa (Y, g, q, itmax=NULL), throws_error())
})
test_that("mfa with itmax neg", {
  expect_that(mfa (Y, g, q, itmax=-1), throws_error())
})
a <- mfa (Y, g, q, init_para = model, sigma_type = model$sigma_type, D_type = model$D_type)
test_that("mfa with init_para as past model", {
  expect_that(a, is_a("mfa"))
})
