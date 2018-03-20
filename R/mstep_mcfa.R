Mstep.mcfa <- function (Y, g, q, pivec, A, xi, omega, D, ...) {

p <- ncol(Y)
if (is.null(p))
  p <- 1
n <- nrow(Y)
p <- ncol(Y)

inv_D <- diag(1 / diag(D))
A1 <- array(0, c(p, q))
A2 <- array(0, c(q, q))
Di <- array(0, p)

inv_D <- diag(1 / diag(D))
tau <- tau.mcfa(Y, g, q, pivec, A, xi, omega, D)
#  inv_D %*% A
inv_D_A <- sweep(A, 1, diag(inv_D), '*')
for (i in 1 : g) {

  # gamma = (A Omega_i A^T + D) ^{-1} A Omega_i
  # using Woodbury Identity
  gamma <- (inv_D - inv_D_A %*% 
              chol.inv(chol.inv(omega[,, i]) +  t(A) %*% inv_D_A) %*%
                t(inv_D_A)) %*% A %*% omega[,, i]

  ti <- sum(tau[, i])
  xi_i <- xi[, i, drop = FALSE]
  
  # tau_ij * y_j  
  tY <- sweep(Y, 1, tau[, i], '*')
  
  # y_j - A xi_i
  Y_Axi_i <- sweep(Y, 2, A %*% xi_i, '-')

  # tau_ij * (y_j - A xi_i)
  tY_Axi_i <- sweep(Y_Axi_i, 1, tau[, i], '*')
  
  # xi_i = xi_i + gamma^T \sum{ tau_ij (y_j - A xi_i) } / sum {tau_ij}
  xi[, i] <- xi_i + t(gamma) %*% as.matrix(colSums(tY_Axi_i)) / ti
  
  # 
  omega[,, i] <- (diag(q) - t(gamma) %*% A) %*% omega[,, i] +
                  t(gamma) %*% t(Y_Axi_i) %*% tY_Axi_i %*% gamma / ti -
                  (xi_i - xi[,i]) %*% t(xi_i - xi[, i])
  # tau_ij yj xi_i^T +   tau_ij yj (yj - A xi)^T gamma
  A1 <- A1 + colSums(tY) %*% t(xi_i) + t(Y) %*% tY_Axi_i %*% gamma
  #
  A2 <- A2 + (omega[,, i] + xi[, i] %*% t(xi[, i])) * ti

  #Di <- Di + diag(t(Y) %*% tY)
  Di <- Di + colSums(Y * tY)

  pivec[i] <- ti / n
}

A <- try(A1 %*% chol.inv(A2))
if (class(A) == "try-error") {

  model <- "tried to invert an ill-conditioned or a singular matrix"
  class(model) <- "error"
  return(model)
}

D <- diag(Di - rowSums((A %*% A2) * A )) / n

model <- list(g = g, q = q, pivec = pivec, A = A,
              xi = xi, omega = omega, D = D)
return(model)
}
