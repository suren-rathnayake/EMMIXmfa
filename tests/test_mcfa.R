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
A1 <- rbind(c(0.5, 0),   c(-0.9, 0),    c(0.3, 0),
             c(0.6, 0.8), c(0.2, -0.7), c(-0.7, 0.5),
             c(0, 0.6), c(0, -0.4), c(0, 0.3), c(0, -0.5))
A2 <- matrix(rep(0, p2*q), c(p2, q))
xi1 <- c(0, 2.5)
xi2 <- c(-2.5, 0)
xi3 <- c(2.5, 0)
xi4 <- c(0, -2.5)
xi5 <- c(0, 0)
Om1 <- rbind(c(0.1, 0), c(0, 0.45))
Om2 <- rbind(c(0.45, 0), c(0, 0.1))
Om3 <- rbind(c(0.45, 0), c(0, 0.1))
Om4 <- rbind(c(0.1, 0), c(0, 0.45))
Om5 <- rbind(c(1, 0.9), c(0.9, 1)) 	
A <- rbind(A1, A2)
D <- diag(c(runif(p1, 0.1, 0.3), runif(p2, 0.3, 0.8)))
pivec <- c(pi1, pi2, pi3, pi4, pi5)
mu <- array(NA, c(p, g))
mu <- cbind(A%*%xi1, A%*%xi2, A%*%xi3, A%*%xi4, A%*%xi5)

sigma <- array(NA, c(p, p, g))
sigma[,, 1] <- A %*% Om1 %*% t(A) + D
sigma[,, 2] <- A %*% Om2 %*% t(A) + D
sigma[,, 3] <- A %*% Om3 %*% t(A) + D
sigma[,, 4] <- A %*% Om4 %*% t(A) + D
sigma[,, 5] <- A %*% Om5 %*% t(A) + D

Y <- array(NA, c(n, p))
cls <- array(NA, c(n, 1))

for( i in 1 : n)  {
  Pi <- 0
  r <- runif(1) 
  for( j in 1 : g) {
    Pi <- Pi + pivec[j]
    if(r < Pi) {
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

context("mcfa")

model <- mcfa(Y, g, q, nkmeans = 2, nrandom = 2, tol = 1.e-5)
expect_that(model, is_a("mcfa"))
expect_that(model, is_a("emmix"))

expect_that(g, equals(model$g))
expect_that(q, equals(model$q))

expect_that(g, equals(length(model$pivec)))
expect_that(1, equals(sum(model$pivec)))

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

expect_that(n, equals(nrow(model$tau)))
expect_that(g, equals(ncol(model$tau)))

expect_that(n, equals(nrow(model$Umean)))
expect_that(q, equals(ncol(model$Umean)))

expect_that(n, equals(nrow(model$Uassign)))
expect_that(q, equals(ncol(model$Uassign)))

dim_U <- dim(model$Uscores)
expect_that(n, equals(dim_U[1]))
expect_that(q, equals(dim_U[2]))
expect_that(g, equals(dim_U[3]))

expect_that(n, equals(length(model$clust)))

g <- 1
q <- 1

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

expect_that(n,   equals(nrow(model$Umean)))
expect_that(q,   equals(ncol(model$Umean)))

expect_that(n,   equals(nrow(model$Uassign)))
expect_that(q,   equals(ncol(model$Uassign)))

dim_U <- dim(model$Uscores)
expect_that(n,   equals(dim_U[1]))
expect_that(q,   equals(dim_U[2]))
expect_that(g,   equals(dim_U[3]))

expect_that(n,   equals(length(model$clust)))

context("mctfa")

q <- 1
model <- mctfa(Y, g, q, nkmeans = 2, nrandom = 2, tol = 1.e-5)
expect_that(model, is_a("mctfa"))
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

expect_that(n, equals(nrow(model$tau)))
expect_that(g, equals(ncol(model$tau)))

expect_that(n, equals(nrow(model$Umean)))
expect_that(q, equals(ncol(model$Umean)))

expect_that(n, equals(nrow(model$Uassign)))
expect_that(q, equals(ncol(model$Uassign)))

dim_U <- dim(model$Uscores)
expect_that(n, equals(dim_U[1]))
expect_that(q, equals(dim_U[2]))
expect_that(g, equals(dim_U[3]))

expect_that(n, equals(length(model$clust)))

expect_that(g, equals(length(model$v)))

q <- 2
model <- mctfa(Y, g, q, nkmeans = 2, nrandom = 2, tol = 1.e-5,
               clust = model$clust, conv_measure = 'diff', df_update = FALSE)

expect_that(model, is_a("mctfa"))
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

expect_that(n,   equals(nrow(model$Umean)))
expect_that(q,   equals(ncol(model$Umean)))

expect_that(n,   equals(nrow(model$Uassign)))
expect_that(q,   equals(ncol(model$Uassign)))

dim_U <- dim(model$Uscores)
expect_that(n,   equals(dim_U[1]))
expect_that(q,   equals(dim_U[2]))
expect_that(g,   equals(dim_U[3]))

expect_that(n,   equals(length(model$clust)))

expect_that(g, equals(length(model$v)))

expect_that(rep(30, g), is_equivalent_to(model$v))


fac <- factor_scores(model, Y)

expect_that(fac, is_a("list"))
expect_that(model, is_a("emmix"))

expect_that(dim(fac$Uscores)[1], equals(n))
expect_that(dim(fac$Uscores)[2], equals(q))
expect_that(dim(fac$Uscores)[3], equals(g))

expect_that(dim(fac$Uassign)[1], equals(n))
expect_that(dim(fac$Uassign)[2], equals(q))

expect_that(dim(fac$Umean)[1], equals(n))
expect_that(dim(fac$Umean)[2], equals(q))


