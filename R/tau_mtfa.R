tau.mtfa <- function(Y, g, q, pivec, B, mu, D, sigma_type, D_type, v, ...) {

if (!is.matrix(Y))
	Y <- as.matrix(Y)
p <- ncol(Y)
n <- nrow(Y)
Fji <- W <- array(NA, c(n,g))

if (sigma_type == 'common') {

	if (D_type == 'common') {
		inv_D <- diag(1 / diag(D))
    B_inv_D <- B * diag(inv_D)
    inv_O <- try(chol.inv(diag(q) +  t(B_inv_D) %*% B))
    if (class(inv_O) == "try-error")
    	return(loglike <-
                 paste('Ill-conditioned or singular sigma'))
  	inv_S <- try(inv_D - B_inv_D %*% inv_O %*% t(B_inv_D))
    if (class(inv_S) == "try-error")
    	return(loglike <-
                paste('Ill-conditioned or singular sigma'))

    for (i in 1 : g) {

			logdetS <- sum(log(diag(D))) -  log(det(inv_O))
      mahal_dist <- mahalanobis(Y, mu[, i, drop=FALSE], inv_S, TRUE)
      Fji[,i] <- lgamma(0.5 * (v[i] + p)) -
                    0.5 * logdetS -  0.5 * p * log(pi * v[i]) -
                    lgamma(0.5 * v[i]) -
                    0.5 * (v[i] + p) * log(1 + mahal_dist / v[i])
      W[, i] <- (v[i] + p) / (v[i] + mahal_dist)
		}
  }
} else {

	if (D_type == 'common') {

		inv_D <- diag(1 / diag(D))
    for (i in 1 : g) {

      B_inv_D <- B[,, i] * diag(inv_D)
      inv_O <- try(chol.inv(diag(q) +  t(B_inv_D) %*% B[,, i]))
      if (class(inv_O) == "try-error")
      		return(loglike <-
                 paste('Ill-conditioned or singular Sigma[,', i, ']'))
      inv_S <- try(inv_D - B_inv_D %*% inv_O %*% t(B_inv_D))
      if (class(inv_S) == "try-error")
          return(loglike <-
                   paste('Ill-conditioned or singular Sigma[,', i, ']'))
      logdetS <- sum(log(diag(D))) -  log(det(inv_O))
      mahal_dist <- mahalanobis(Y, mu[, i, drop = FALSE], inv_S, TRUE)
      Fji[,i] <- lgamma(0.5 * (v[i] + p)) -
                    0.5 * logdetS -  0.5 * p * log(pi * v[i]) -
                    lgamma(0.5 * v[i]) -
                    0.5 * (v[i] + p) * log(1 + mahal_dist / v[i])
      W[, i] <- (v[i] + p) / (v[i] + mahal_dist)
    }
  }

  if (D_type == 'unique') {
		for (i in 1 : g) {

        inv_D <- diag(1 / diag(D[,, i]))
        B_inv_D <- B[,, i] * diag(inv_D)
        inv_O <- try(chol.inv(diag(q) +  t(B_inv_D) %*% B[,, i]))
        if (class(inv_O) == "try-error")
					return(loglike <-
                 paste('Ill-conditioned or singular Sigma[,', i, ']'))
        inv_S <- try(inv_D - B_inv_D %*% inv_O %*% t(B_inv_D))
        if (class(inv_S) == "try-error")
         	return(loglike
							<- paste('Ill-conditioned or singular Sigma[,', i, ']'))
        logdetS <- sum(log(diag(D[,, i]))) - log(det(inv_O))
        mahal_dist <- mahalanobis(Y, mu[, i, drop=FALSE], inv_S, TRUE)
        Fji[,i] <- lgamma(0.5 * (v[i] + p)) -  0.5 * logdetS -
                    0.5 * p * log(pi * v[i]) -
                    lgamma(0.5 * v[i]) - 0.5 * (v[i] + p) *
										log(1 + mahal_dist / v[i])
        W[, i]  <- (v[i] + p) / (v[i] + mahal_dist)
    }
  }
}

Fji <- sweep(Fji, 2, log(pivec), '+')
Fjmax <- apply(Fji, 1, max)
Fji <- sweep(Fji, 1, Fjmax, '-')
Fji <- exp(Fji)
tau <- sweep(Fji, 1, rowSums(Fji), '/')
return (list(tau = tau, W = W))
}
