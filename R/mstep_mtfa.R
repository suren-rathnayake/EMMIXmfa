Mstep.mtfa <- function(Y, g, q, pivec, B, mu, D, sigma_type, D_type,
                                            v, df_update, tau, W, ...) {

args <- list(...)
if(!("v.min" %in% names(args)))
  v.min <- 0.0001
if(!("v.max" %in% names(args)))
  v.max <- 200
p <- ncol(Y)
n <- nrow(Y)
pivec <- t(colSums(tau)/n)
tauW <- tau* W

for (i in 1 : g) {
  mu[, i] <-
           matrix(colSums(sweep(Y,  MARGIN=1, tauW[, i], '*')))/sum(tauW[, i])

  if (df_update) {

    f <- function(x) {
          -digamma(0.5* x) + log(0.5* x) + 1 +
            sum(tau[,i]* (log(W[,i]) - W[,i]))/sum(tau[,i]) +
            digamma(0.5* (v[i]+p)) - log(0.5* (v[i]+p))
    }

    fmin <- try(f(v.min))
    fmax <- try(f(v.max))

    if (class(fmin) == class(fmax)) {
      if ((class(fmin) == "numeric") && is.finite(fmin)
                        && is.finite(fmax)) {

        if (sign(fmin) == sign(fmax)) {

          if(abs(fmin) < abs(fmax)) {v[i] <- v.min} else {v[i] <- v.max}
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

TW <- tau.mtfa(Y, g, q, pivec, B, mu, D, sigma_type, D_type, v)
tau <- TW$tau
W <- TW$W
tauW <- tau * W

if (sigma_type == 'common') {

  inv_D <- diag(1/diag(D))
  B_inv_D <- B *diag(inv_D)
  beta <- (inv_D - B_inv_D %*%
            (chol.inv(diag(q) +  t(B_inv_D) %*% B)) %*%
            t(B_inv_D)) %*% B
  H <- (diag(q) - t(beta) %*% B)
  V <- matrix(0, nrow = p, ncol = p)

  for (i in 1 : g) {
    Ymu <- sweep(Y, 2,  mu[, i, drop=FALSE], '-')
    V   <- V + t(Ymu) %*% sweep(Ymu, 1, tauW[, i], '*')
  }
  B <- try(V %*% beta %*% solve(n*H + t(beta) %*% V %*% beta))

  if (class(B) == 'try-error') {
    ERR <- "inversion of a singular matrix"
    class(ERR) <- "error"
    return(ERR)
  }

  D <- diag(diag(V) +
           rowSums(B %*% (n * H + t(beta) %*% V %*% beta) * B) -
            2 * diag(B %*% t(beta) %*% V)) / n
}

if (sigma_type == "unique") {

  if (D_type == "unique") {

    for(i in 1 : g) {
			inv_D <- diag(1 / diag(D[,, i]))
      B_inv_D <- B[,, i] *diag(inv_D)
      beta <- (inv_D - B_inv_D %*%
                  (chol.inv(diag(q) +  t(B_inv_D) %*% B[,, i])) %*%
                  t(B_inv_D)) %*% B[,, i]
      Ymu <- sweep(Y, 2, mu[, i, drop = FALSE], '-')
      V <- t(Ymu) %*% sweep(Ymu, 1, tauW[, i], '*')
      H <- diag(q) - t(beta) %*% B[,, i]
      tau_i <- sum(tau[,i])
      B2 <- try(solve(H * tau_i +  t(beta) %*% V %*% beta))
      if (class(B2) == 'try-error') {
        ERR <- "inversion of a singular matrix"
        class(ERR) <- "error"
        return(ERR)
      }

      B[,, i] <- V %*% beta %*% B2
		  D[,, i] <-  diag((diag(V) - 2*diag(B[,, i] %*% t(beta) %*% V) +
                          rowSums((B[,, i] %*% (H*tau_i +
                          t(beta) %*% V %*% beta)) * B[,, i])) / tau_i)
    }
  }

  if (D_type == "common") {

    inv_D <- diag(1 / diag(D))
    D <- array(NA, c(p,  p, g))

    for (i in 1 : g) {
      B_inv_D <- sweep(as.matrix(B[,, i]), 1, diag(inv_D), "*")
      beta <- (inv_D - B_inv_D %*%
               (chol.inv(diag(q) +  t(B_inv_D) %*% B[,, i])) %*%
               t(B_inv_D)) %*% B[,, i]
      Ymu <- sweep(Y, 2, mu[, i, drop=FALSE], '-')
      V <- t(Ymu) %*% sweep(Ymu, 1, tauW[, i], '*')
      H <- diag(q) - t(beta) %*% B[,, i]
      tau_i <- sum(tau[,i])
      B2 <- try(solve(H * tau_i +  t(beta) %*% V %*% beta))
      if (class(B2) == 'try-error') {
        ERR <- "inversion of a singular matrix"
        class(ERR) <- "error"
        return(ERR)
      }
      B[,, i] <- V %*% beta %*% B2
      D[,, i] <- diag(diag(V) - 2*diag(B[,, i] %*% t(beta) %*% V) +
                  rowSums((B[,, i] %*% (H * tau_i
                    + t(beta) %*% V %*% beta)) * B[,, i]))
   }
   D <- diag(diag(apply(D, c(1,2), sum))/n)
 }
}

model <- list(g = g, q = q,  pivec = pivec, B = B, mu = mu, D = D,
              sigma_type = sigma_type, D_type = D_type, v = v,
              df_update = df_update)
class(model) <- "mtfa"
return(model)
}
