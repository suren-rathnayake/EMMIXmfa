factor_scores_mctfa <- function(Y, g, q, pivec, A, xi, omega, D, v,
                          tau = NULL, clust = NULL, ...) {

if(!is.matrix(Y))
  Y <- as.matrix(Y)
p <- ncol(Y)
n <- nrow(Y)
U <- array(0, c(n, q, g))
gamma <- array(0, c(p, q))
inv_D <- diag(1/diag(D))

for (i in 1 : g) {
  gamma <- (inv_D - sweep(A, 1, diag(inv_D), '*') %*%
            chol.inv(chol.inv(omega[,, i]) +
            t(A * diag(inv_D)) %*% A) %*%
            sweep(t(A), 2, diag(inv_D), '*')) %*% A %*% omega[,, i]
  U[,, i] <- sweep(sweep(Y, 2, A %*% xi[, i, drop = F], '-') %*% gamma,
                 2, xi[, i, drop = F], "+")
}
if (is.null(tau)) {
  tau <- tau.mctfa(Y, g, q, pivec, A, xi, omega, D, v)
}
if (is.null(clust)) {
  clust <- apply(tau, 1, which.max)
}
UC <- array(0, c(n, q))
Umean <- array(0, c(n, q))
for (i in 1 : n) {
  UC[i, ] <- U[i,, clust[i]]
  Umean[i, ] <- tau[i, ] %*% t(matrix(U[i,, ], c(q, g)))
}
return(list(Uscores = U, Uclust = UC, Umean = Umean))
}
