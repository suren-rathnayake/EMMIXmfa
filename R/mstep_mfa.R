Mstep.mfa <- function(Y, g, q, pivec, B, mu, D, 
                      sigma_type, D_type, tau, ...) {
                      
p <- ncol(Y) 
n <- nrow(Y)
pivec <- t(colSums(tau) / n)

for (i in 1:g)
  mu[, i] <- matrix(colSums(sweep(Y,  MARGIN=1, tau[, i], '*'))) / 
                    sum(tau[, i])

tau <- tau.mfa(Y = Y, g = g, q = q, pivec = pivec, B = B, mu = mu, 
                D = D, sigma_type = sigma_type, D_type = D_type)

if (sigma_type == 'common') {

  inv_D <- diag(1 / diag(D))
  B_inv_D <- B * diag(inv_D)
  beta <- (inv_D - B_inv_D %*%
             (chol.inv(diag(q) +  t(B_inv_D) %*% B)) %*%
             t(B_inv_D)) %*% B
  W <- (diag(q) - t(beta) %*% B)
  V <- matrix(0, nrow = p, ncol = p)
  
  for(i in 1 : g) {
     Ymu <- sweep(Y, 2, mu[, i, drop=FALSE], '-')
     V   <- V + t(Ymu) %*% sweep(Ymu, 1, tau[, i], '*')
  }
  B <- try(V %*% beta %*% solve(n * W + t(beta) %*% V %*% beta))

  if(class(B) == 'try-error') {
    ERR <- "inversion of a singular matrix"
    class(ERR) <- "error"
    return(ERR)
  }
  D <- diag(diag(V) +
       rowSums(B %*% (n * W + t(beta) %*% V %*% beta) * B) -
    	 2 * diag(B %*% t(beta) %*% V)) / n
}

if(sigma_type == "unique") {

	if (D_type == "unique") {

		for(i in 1 : g) {
    
			inv_D <- diag(1 / diag(D[,, i]))
  		B_inv_D <- B[,, i] * diag(inv_D)
  		beta <- (inv_D - B_inv_D %*%
               (chol.inv(diag(q) +  t(B_inv_D) %*% B[,, i])) %*%
               t(B_inv_D)) %*% B[,, i]
    	Ymu <- sweep(Y, 2, mu[, i, drop=FALSE], '-')
     	V <- t(Ymu) %*% sweep(Ymu, 1, tau[, i], '*')
			W <- diag(q) - t(beta) %*% B[,, i]
			tau_i <- sum(tau[,i])
      B2 <- try(solve(W * tau_i + t(beta) %*% V %*% beta))
      
      if(class(B2) == 'try-error') {
        ERR <- "inversion of a singular matrix"
        class(ERR) <- "error"
        return(ERR)
      }
      
      B[,, i] <- V %*% beta %*% B2
			D[,, i] <- diag((diag(V) - 2 * diag(B[,, i] %*% t(beta) %*% V) +
                  rowSums((B[,, i] %*% (W * tau_i + 
                  t(beta) %*% V %*% beta)) * B[,, i])) / tau_i)
  	}
  }

  if (D_type == "common") {
 
    inv_D <- diag(1 / diag(D))
    D <- array(NA, c(p, p, g))
    for (i in 1 : g) {
    
      B_inv_D <- sweep(as.matrix(B[,, i]), 1, diag(inv_D), "*")
      beta <- (inv_D - B_inv_D %*%
               (chol.inv(diag(q) +  t(B_inv_D) %*% B[,, i])) %*%
               t(B_inv_D)) %*% B[,, i]
      Ymu <- sweep(Y, 2, mu[, i, drop = FALSE], '-')
      V <- t(Ymu) %*% sweep(Ymu, 1, tau[, i], '*')
      W <- diag(q) - t(beta) %*% B[,, i]
      tau_i <- sum(tau[, i])
      B2 <- try(solve(W * tau_i + t(beta) %*% V %*% beta))
      if(class(B2) == 'try-error') {
        ERR <- "inversion of a singular matrix"
        class(ERR) <- "error"
        return(ERR)
      }
      
      B[,, i] <- V %*% beta %*% B2
      D[,, i] <- diag(diag(V) - 2 * diag(B[,, i] %*% t(beta) %*% V) +
                  rowSums((B[,, i] %*% (W * tau_i + t(beta) %*% V %*% beta)) * 
                  B[,, i]))
    }
    D <- diag(diag(apply(D, c(1, 2), sum))/n)
 }
}

model <- list(g = g, q = q,  pivec = pivec, B = B, mu = mu, D = D, 
              sigma_type = sigma_type, D_type = D_type)
class(model) <- "mfa"
return(model)
}
