Mstep.mctfa <- function(Y, g, q, pivec, A, xi, omega, D, v, df_update,
                        tau, W, ...) {

args <- list(...)
if (!("v.min" %in% names(args)))
  v.min <- 0.0001
if (!("v.max" %in% names(args)))
  v.max <- 200

p <- ncol(Y)
n <- nrow(Y)

tauW <- tau * W
inv_D <- diag(1 / diag(D))

A1 <- array(0, c(p, q))
A2 <- array(0, c(q, q))
Di <- array(0, c(p))
#  inv_D %*% A
inv_D_A <- sweep(A, 1, diag(inv_D), '*')

for (i in 1 : g) {

  # gamma = (A Omega_i A^T + D) ^{-1} A Omega_i
  # using Woodbury Identity
  gamma <- (inv_D - inv_D_A %*%
               chol.inv(chol.inv(omega[,, i]) +  t(A) %*% inv_D_A) %*%
                 t(inv_D_A)) %*% A %*% omega[,, i]

#  gamma <- (inv_D - inv_D_A %*% #sweep(A, 1, diag(inv_D), '*') %*%
#             chol.inv(chol.inv(omega[,, i]) +
#             t(A * diag(inv_D)) %*% A) %*%
#             sweep(t(A), 2, diag(inv_D), '*')) %*% A %*% omega[,, i]

  ti <- sum(tau[, i])
  twi <- sum(tauW[, i])
  xi_i <- xi[, i, drop = F]
  
  # tau_ij w_ij y_j
  twY <- sweep(Y, 1, tauW[,i], '*')

  # y_j - A xi_i
  Y_A_xi_i <- sweep(Y, 2, A %*% xi_i, '-')

  # tau_ij w_ij (y_j - A xi_i)
  twY_A_xi_i <- sweep(Y_A_xi_i, 1, tauW[, i], '*')

  # xi_i + gamma^T \sum{ tau_ij w _ij(y_j - A xi_i) } / sum {tau_ij w_ij}
  xi[, i] <- xi_i + t(gamma) %*% as.matrix(colSums(twY_A_xi_i)) / twi

  #
  omega[,, i] <- (diag(q) - t(gamma) %*% A) %*% omega[,, i] +
                    t(gamma) %*% t(Y_A_xi_i) %*% twY_A_xi_i %*% gamma / ti -
                     (xi_i - xi[, i]) %*% t(xi_i - xi[, i])
  
  # tau_ij w_ij yj xi_i^T + w_ij tau_ij yj (yj - A xi)^T gamma
  A1 <- A1 + colSums(twY) %*% t(xi_i) + t(Y) %*% twY_A_xi_i %*% gamma
  #  
  A2 <- A2 + omega[,, i]* ti + xi[, i] %*% t(xi[, i]) * twi
  
  Di <- Di + colSums(Y * twY) 
  
  pivec[i] <- ti / n
  
  if (df_update) {
    f <- function(x) {
          -digamma(x / 2) + log(x / 2) + 1 +
          sum(tau[, i] * (log(W[, i]) -
          W[, i])) / sum(tau[, i]) +
          digamma((v[i] + p) / 2) - log((v[i] + p) / 2) }

    fmin <- try(f(v.min))
    fmax <- try(f(v.max))
    if (class(fmin) == class(fmax)) {
      if ((class(fmin) == "numeric") && is.finite(fmin) && is.finite(fmax)) {
        if (sign(fmin) == sign(fmax)) {
          if (abs(fmin) < abs(fmax)) {v[i] <- v.min} else {v[i] <- v.max}
        } else {
          v[i] <- uniroot(f, c(v.min, v.max))$root
        }
      } else {
        model <- "Error estimating DOF"
        class(model) <- "error"
        return(model)
      }
    } else {
      model <- "Error estimating DOF"
      class(model) <- "error"
      return(model)
    }
  }
}

A <- try(A1 %*% chol.inv(A2))
if (class(A) == "try-error") {
  model <- "ill-conditioned or a singular matrix"
  class(model) <- "error"
  return(model)
}
D <- diag(Di - rowSums((A %*% A2) * A )) / n
model <- list(g = g, q = q, pivec = pivec, A = A, xi = xi, omega = omega,
              D = D, v = v, df_update = df_update)
class(model) <- "mctfa"
return(model)
}
