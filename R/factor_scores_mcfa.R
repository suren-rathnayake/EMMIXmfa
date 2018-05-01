factor_scores_mcfa <- function(Y, g, q, pivec, A, xi, omega, D,
                               tau = NULL, clust = NULL, ...) {

p <- ncol(Y)
if(is.null(p))
  p <- 1
n <- nrow(Y)

U <- array(0, c(n, q, g))
gamma <- array(0, c(p, q, g))

invD <- diag(1 / diag(D))
for (i in 1 : g) {
  gamma[,, i] <- (invD - invD %*% A %*%
                   chol.inv(chol.inv(omega[,, i]) +
                   t(A) %*% invD %*% A) %*%
                   t(A) %*% invD) %*% A %*% omega[,, i]

  U[,, i] <- t(replicate(n , xi[, i, drop = FALSE], 'matrix'))

  U[,, i] <- U[,, i] +
              (Y - t(replicate(n, A%*% xi[, i, drop = FALSE], "matrix"))) %*%
              as.matrix(gamma[,, i])
}

if(is.null(tau))
  tau <- tau.mcfa(Y, g, q, pivec, A, xi, omega, D)

if(is.null(clust))
  clust <- apply(tau, 1, which.max)

UC <- array(0, c(n, q))
Umean <- array(0, c(n, q))

for (i in 1 : n) {
  UC[i, ] <- U[i,, clust[i]]
  Umean[i, ] <- tau[i, ] %*% t(matrix(U[i,, ], c(q, g)))
}

return(list(Uscores = U, Uclust = UC, Umean = Umean))
}
