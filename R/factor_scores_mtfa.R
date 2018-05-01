factor_scores_mtfa <- function (Y, g, q, pivec, B, mu, D, sigma_type, D_type,
                                v, tau = NULL, clust = NULL, ...) {

if (!is.matrix(Y)) {
 Y <- as.matrix(Y)
}
p <- ncol(Y)
n <- nrow(Y)
U <- array(0, c(n, q, g))

if (sigma_type == "common") {
  gamma <- matrix(0, nrow = p, ncol = q)
  for (i in 1 : g) {
    inv_D <- diag(1/diag(D))
    B_inv_D <- B *diag(inv_D)
    gamma <- (inv_D - B_inv_D %*%
                  (chol.inv(diag(q) +  t(B_inv_D) %*% B)) %*%
                  t(B_inv_D)) %*% B
    U[,, i] <- sweep(Y, 2, mu[, i, drop=FALSE], '-') %*% gamma
  }
} else {
  gamma <- array(0, c(p, q))
  for (i in 1 : g) {

    if (D_type == 'common') {
      inv_D <- diag(1/diag(D))
    } else {
      inv_D <- diag(1/diag(D[,, i]))
    }
    B_inv_D <- B[,, i] * diag(inv_D)
    gamma <- (inv_D - B_inv_D %*%
                  (chol.inv(diag(q) +  t(B_inv_D) %*% B[,, i])) %*%
                  t(B_inv_D)) %*% B[,, i]
    U[,, i] <- sweep(Y, 2, mu[, i, drop=FALSE], '-') %*% gamma
  }
}
if (is.null(tau)) {
  tau <- tau.mtfa(Y, g, q, pivec, B, mu, D, sigma_type, D_type, v)$tau
}
if (is.null(clust)) {
  clust <- apply(tau, 1, which.max)
}
UC <- array(0, c(n, q))
Umean <- array(0, c(n, q))
for (j in 1 : n) {
  UC[j, ] <- U[j,, clust[j]]
  Umean[j, ] <- tau[j,] %*% t(matrix(U[j,,], c(q, g)))
}

return(list(Uscores = U, Uclust = UC, Umean = Umean))
}
