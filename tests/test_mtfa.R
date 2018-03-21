rm(list=ls())
require(testthat)
require(EMMIXmfa)

library(mvtnorm)
p1 <- 10
p2 <- 100
p = p1 + p2
q = 2
g = 5
n <- 50
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

for( i in 1 : n) 
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

g <- 3
q <- 2
p <- ncol(Y)
n <- nrow(Y)
context("mtfa")
model <- mtfa(Y, g, q, nkmeans = 2, nrandom = 2, tol = 1.e-5, 
             sigma_type = "c", D_type = "common")
expect_that(model, is_a("mtfa"))
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



expect_that(n,   equals(length(model$clust)))

test_that("mtfa with init_para as past model", {
  expect_that(model, is_a("mtfa"))
})

test_that("mtfa with init_para as past model", {
  expect_that(model, is_a("emmix"))
})

expect_that(mtfa(Y, g, q, itmax = 0),
            throws_error("Maximum number of iterations, itmax, must be greather than one."))

expect_that(mtfa(Y, g, q),
            throws_error("The data must have more than one variable."))


expect_that(mtfa(Y, g, q = p),
            throws_error("The number of factors must be less than the number of variables."))

expect_that(mtfa(Y, g = g, q = -1),
            throws_error("q must be a positive integer."))

expect_that(mtfa(Y, g = -2, q = q),
            throws_error("g must be a positive integer."))


expect_that(mtfa(Y, g = 2.5, q = q),
            throws_error("g must be a positive integer."))


expect_that(mtfa(Y, g = g, q = 1.5),
            throws_error("q must be a positive integer."))

expect_that(mtfa(Y, g = g, q = q, sigma_type == 'common', D_type == 'unique'),
            throws_error("D_type = 'unique' not available with sigma_type = 'common'."))

expect_that(mtfa(Y, g = g, q = q, sigma_type == 'd', D_type == 'unique'),
            throws_error("sigma_type needs to be either 'unique' or 'common'"))

expect_that(mtfa(Y, g = g, q = q, sigma_type == 'common', D_type == 'a'),
            throws_error("D_type needs to be either 'unique' or 'common'."))

expect_that(mtfa(Y, g = g, q = q, conv_measure = "adiff"),
            throws_error("conv_measure needs to be either 'diff' or 'ratio'."))

expect_that(mtfa(Y, g = g, q = q, warn_messages = 1),
            throws_error("warn_messages must either be TRUE or FALSE."))

expect_that(mtfa(Y, g = g, q = q, init_para = model, sigma_type = "c"),
            throws_error("`init_para$sigma_type` is not same as `sigma_type`."))

expect_that(mtfa(Y, g = g, q = q, init_para = model, D_type = "u"),
            throws_error("`init_para$D_type` is not same as `D_type`."))

expect_that(mtfa(Y, g = g, q = q, df_update = "u"),
            throws_error("df_update is either TRUE or FALSE'."))

expect_that(mtfa(Y, g = g, q = q, df_update = TRUE, df_init = letters[1:g]),
            throws_error("df_init has a non-numeric element."))



model <- mtfa(Y, g, q, nkmeans = 20, nrandom = 20, tol = 1.e-5,
             sigma_type = "unique", D_type = "unique", df_update = FALSE)
expect_that(model, is_a("mtfa"))
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

model <- mtfa(Y, g, q, nkmeans = 2, nrandom = 2, tol = 1.e-5,
             sigma_type = "common", D_type = "common")
expect_that(model, is_a("mtfa"))
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




context("mtfa")
model <- mtfa(Y, g, q, nkmeans = 2, nrandom = 2, tol = 1.e-5, 
             sigma_type = "unique", D_type = "common")
expect_that(model, is_a("mtfa"))
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

model <- mtfa(Y, g, q, nkmeans = 2, nrandom = 2, tol = 1.e-5,
             sigma_type = "unique", D_type = "unique")
expect_that(model, is_a("mtfa"))
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

model <- mtfa(Y, g, q, nkmeans = 2, nrandom = 2, tol = 1.e-5,
             sigma_type = "common", D_type = "common")
expect_that(model, is_a("mtfa"))
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

fac <- factor_scores(model, Y)

expect_that(fac, is_a("list"))
expect_that(model, is_a("emmix"))
expect_that(dim(fac$U)[1],   equals(n))
expect_that(dim(fac$U)[2],   equals(q))
expect_that(dim(fac$U)[3],  equals(g))

expect_that(dim(fac$UC)[1],   equals(n))
expect_that(dim(fac$UC)[2],   equals(q))

expect_that(dim(fac$Fmat)[1],   equals(n))
expect_that(dim(fac$Fmat)[2],   equals(q))


