tau.mctfa <- function(Y, g, q, pivec, A, xi, omega, D, v, ...) {
  
if(!is.matrix(Y)) 
  Y <- as.matrix(Y)
p <- ncol(Y)
n <- nrow(Y)

Fji <- array(NA, c(n,g))
inv_D <- diag(1/diag(D))
#  inv_D %*% A
inv_D_A <- sweep(A, 1, diag(inv_D), '*')

for(i in 1 : g) {
  inv_O <- try(chol.inv(chol.inv(omega[,, i]) +
                         t(inv_D_A) %*% A), silent = TRUE)
  if (class(inv_O) == "try-error")
    return(loglike <-
      paste('Ill-conditioned or singular Sigma[,', i, ']'))

  inv_S <- try(inv_D - inv_D_A %*% inv_O %*% t(inv_D_A), silent = TRUE)
  if (class(inv_S) == "try-error") {
      return(loglike <- paste('Ill-conditioned or singular Sigma_',
                                                  i, sep = ''))
  }  

  logdetS <- log(det(as.matrix(omega[,,i]))) + sum(log(diag(D))) -
              log(det(inv_O))

  mahal_dist <- mahalanobis(Y, t(A %*% xi[, i, drop = FALSE]), inv_S, TRUE)
  
  Fji[,i] <- lgamma(v[i] / 2 + p / 2) - 0.5 * logdetS -
              0.5 * p* log(pi * v[i]) - lgamma(v[i] / 2) - 
              0.5 * (v[i] + p) * log(1 + mahal_dist / v[i])
}

Fji <- sweep(Fji, 2, log(pivec), '+')
Fjmax <- apply(Fji, 1, max)
Fji <- sweep(Fji, 1, Fjmax, '-')
Fji <- exp(Fji)
tau <- sweep(Fji, 1, rowSums(Fji), '/')

return(tau)
}
