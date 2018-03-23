rmfa <- function(n, model, ...) {
  
  g <- model$g
  q <- model$q 
  
  n_mix <- sample.int(g, n, replace = TRUE, prob = model$pivec)
  tabn <- table(n_mix)
  dat <- matrix(NA, nrow = n, ncol = ncol(model$D))
  
  if (model$sigma_type == "common") {

    sigma <- with(model, B %*% t(B) + D)
    for (i in 1 : g) {
      mu <- model$mu[, i] 
      dat[which(n_mix == i), ] <- mvtnorm::rmvnorm(tabn[i], mean = mu, 
                                                  sigma = sigma, ...)
    } 
  } else {

    for (i in 1 : g) {

      if (model$D_type == "common") {
        sigma <- with(model, B[,, i] %*%t(B[,, i]) + D)
      } else {
        sigma <- with(model, B[,, i] %*%t(B[,, i]) + D[,, i])
      }

      mu <- model$mu[, i] 

      dat[which(n_mix == i), ] <- mvtnorm::rmvnorm(tabn[i], mean = mu, 
                                                 sigma = sigma, ...)
    } 
  }
  dat
}

rmtfa <- function(n, model, ...) {
  
  g <- model$g
  q <- model$q 
  
  n_mix <- sample.int(g, n, replace = TRUE, prob = model$pivec)
  tabn <- table(n_mix)
  dat <- matrix(NA, nrow = n, ncol = ncol(model$D))
  
  if (model$sigma_type == "common") {

    sigma <- with(model, B %*% t(B) + D)
    for (i in 1 : g) {
      delta <- model$mu[, i] 
      dat[which(n_mix == i), ] <- mvtnorm::rmvt(tabn[i], delta = delta, 
                                                  sigma = sigma, 
                                                  df = model$v[i],
                                                  type = "shifted", ...)
    } 
  } else {

    for (i in 1 : g) {

      if (model$D_type == "common") {
        sigma <- with(model, B[,, i] %*%t(B[,, i]) + D)
      } else {
        sigma <- with(model, B[,, i] %*%t(B[,, i]) + D[,, i])
      }

      delta <- model$mu[, i] 

      dat[which(n_mix == i), ] <- mvtnorm::rmvt(tabn[i], delta = delta, 
                                                 sigma = sigma, 
                                                 df = model$v[i],
                                                 type = "shifted", ...)
    } 
  }
  dat
}

rmcfa <- function(n, model, ...) {
  
  g <- model$g
  q <- model$q 
  
  n_mix <- sample.int(g, n, replace = TRUE, prob = model$pivec)
  tabn <- table(n_mix)
  dat <- matrix(NA, nrow = n, ncol = ncol(model$D))
  
  for (i in 1 : g) {
    
    mu <- with(model, A %*% xi[, i])
    sigma <- with(model, A %*% omega[,, i] %*% t(A) + D)
    
    dat[which(n_mix == i), ] <- mvtnorm::rmvnorm(tabn[i], mean = mu, 
                                                 sigma = sigma, ...)
  }
  
  dat
}

rmctfa <- function(n, model, ...) {
  
  g <- model$g
  q <- model$q
  
  n_mix <- sample.int(g, n, replace = TRUE, prob = model$pivec)
  tabn <- table(n_mix)
  dat <- matrix(NA, nrow = n, ncol = ncol(model$D))
  
  for (i in 1 : g) {
    
    delta <- with(model, A %*% xi[, i])
    sigma <- with(model, A %*% omega[,, i] %*% t(A) + D)
    
    dat[which(n_mix == i), ] <- mvtnorm::rmvt(tabn[i], delta = delta, 
                                              sigma = sigma, 
                                              df = model$v[i], 
                                              type = "shifted", ...)
  }
  
 dat
}


rmix <- function(n, model, ...) {

if (!requireNamespace("mvtnorm", quietly = TRUE)) {
 
  stop("rmix require mvtnorm package. Please `install.packages(mvtnorm)`")
}

if (!(any(class(model) == "mcfa") || any(class(model) == "mctfa") ||
      any(class(model) == "mfa") || any(class(model) == "mtfa"))) {

  stop("STOP:model must be of class mcfa, mctfa, mfa or mtfa")
}

ncomp <- length(model$pivec)
p <- dim(model$D)[1]
sigma <- array(NA, c(p, p, ncomp))
dat <- matrix(NA, nrow = n, ncol = p)

if (any(class(model) == "mcfa")) {
  dat <- rmcfa(n, model, ...)
} 

if (any(class(model) == "mctfa")) {
  dat <- rmctfa(n, model, ...)
} 

if (any(class(model) == "mfa")) {
  dat <- rmfa(n, model, ...)
} 

if (any(class(model) == "mtfa")) {
  dat <- rmtfa(n, model, ...)
} 

dat 
}
