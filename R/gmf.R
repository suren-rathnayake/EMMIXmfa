gmf <- function(Y, q, maxit = 1000, lambda = 0.01, cor_rate = 0.9) {

#Ls <- 1
Y <- t(Y)

p <- nrow(Y)
n <- ncol(Y)

Ls <- sum(Y^2)
#A <- matrix(rnorm(p * q), nrow = p, ncol = q)
#B <- matrix(rnorm(q * n), nrow = q, ncol = n)
A <- rnorm(p * q)
B <- rnorm(q * n)

verr <- rep(0, maxit)
grr  <-  rep(0, maxit)

ab <- .C(cgmf,
         as.double(Y),
         as.double(A),
         as.double(B),
         as.integer(n),
         as.integer(p),
         as.integer(q),
         as.double(lambda),
         as.integer(maxit),
         as.double(Ls),
         as.double(cor_rate),
         as.double(verr),
         as.double(grr),
         PACKAGE='EMMIXmfa')

#Ls <- -Inf
#Y <- t(Y)
#  
#p <- nrow(Y)
#n <- ncol(Y)
#  
#A <- matrix(rnorm(p * q), nrow = p, ncol = q)
#B <- matrix(rnorm(q * n), nrow = q, ncol = n)
#  
#for (m in 1 : maxit) {
#  for (id_p in 1 : p) {
#    for (id_n in 1 : n) {
#      S <- A[id_p, , drop = FALSE] %*% B[, id_n, drop = FALSE]
#      E <- Y[id_p, id_n] - S
#        
#      for (f in 1 : q) {
#        alpha <- A[id_p, f] * B[f, id_n]
#        A[id_p, f] <- A[id_p, f] + lambda * E * B[f, id_n]
#        E <- E + alpha - A[id_p, f] * B[f, id_n]
#        alpha <- A[id_p, f] * B[f, id_n]
#        B[f, id_n] <- B[f, id_n] + lambda * E * A[id_p, f]
#        E <- E + alpha - A[id_p, f] * B[f, id_n]
#      }
#    }
#  }
#    
#  L <- sum((Y - A %*% B)^2) / (n * p)
#
#  if (L > Ls) {
#    Ls <- L
#    lambda <- cor_rate * lambda
#  }
#}
 
A <- matrix(ab[[2]], c(p, q))
B <- matrix(ab[[3]], c(n, q))

return(list(A = A, B = B))
}
