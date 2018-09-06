tau.mcfa <- function(Y, g, q, pivec, A, xi, omega, D, ...) {

	p <- ncol(Y) 
	if (is.null(p)) 
	  p <- 1
	n <- nrow(Y)

	Fji <- array(NA, c(n, g))
	inv_D <- diag(1 / diag(D))
	#  inv_D %*% A
	inv_D_A <- sweep(A, 1, diag(inv_D), '*')
	for (i in 1 : g) {

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

	  # |D + A \Om A^T| using matrix determinant lemma
	  # |D| + |\Om| + |\Om^{-1} + A^T D^{-1} A|
	  logdetD <- log(det(as.matrix(omega[,, i]))) +
	                sum(log(diag(D))) -
	                log(det(inv_O))

	  mahal_dist <- mahalanobis(Y, t(A %*% xi[, i, drop = FALSE]), inv_S, TRUE)

	  Fji[,i] <- -0.5 * mahal_dist - (p / 2) * log(2 * pi) - 0.5 * logdetD
	}

	Fji <- sweep(Fji, 2, log(pivec), '+')
	Fjmax <- apply(Fji, 1, max)
	Fji <- sweep(Fji, 1, Fjmax, '-')
	Fji <- exp(Fji)
	tau <- sweep(Fji, 1, rowSums(Fji), '/')

	return(tau)
}
