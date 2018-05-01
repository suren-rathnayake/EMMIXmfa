print.emmix <- function(x, ...) {
  cat("Call:\n")
  print(x$call)

  cat("\nCoefficients: \n")
  cat("pi_i : ", round(x$pivec, 3), "\n")

  if ((any(class(x) == "mfa") || any(class(x) == "mtfa"))) {
    for(j in 1 : x$g) {
      cat("mu_", j, ":\n")
      print(x$mu[, j])
    }

    if(x$sigma_type == 'common') {
      cat("B: \n ")
      print(x$B)
      cat("diag D: \n")
      print(diag(x$D))
    } else {
      for(j in 1 : x$g) {
        cat("B_",j, ":\n")
        print(x$B[,, j])
      }

      if(x$D_type == 'common')
        cat("diag D: \n", diag(x$D), "\n")

      if (x$D_type == 'unique') {
        for (i in 1 : x$g)
          cat("diag D_", i, ":\n", diag(x$D[,, i]), "\n")
      }
    }

    if(any((class(x) == "mtfa"))) {
      cat("nu: \n", x$v, '\n')
    }
  }

  if((any(class(x) == "mcfa") || any(class(x) == "mctfa"))) {

    cat("A: \n")
    print(x$A)

    for(j in 1 : x$g) {
      cat("xi_",j,":\n")
      print(x$xi[, j])
    }
    for(j in 1:x$g) {
      cat("omega_",j, ":\n")
      print(x$omega[,, j])
    }
    cat("diag D: \n", diag(x$D), "\n")
    if(any((class(x) == "mctfa"))) {
      cat("nu: \n", x$v, "\n")
    }
  }
}

print.mcfa <- function(x, ...) {
  
  class(x) <- c("emmix", class(x))
  print(x)

  # g <- x$g
  # q <- x$q
  # cat("Call:\n")
  # print(x$call)
  # cat("\n", paste0("Coefficients: (for i in 1 to ", g, ")", "\n"))
  # cat("\nMixing Proportions\n")
  # cat("pi_i : ", round(x$pivec, 3), "\n")
  # cat("\nLoading Matrix\n")
  # colnames(x$A) <- paste0("q_",  1 : q)
  # cat("A: \n")
  # print(x$A)
  # cat("\nFactor Means (in columns)\n")
  # colnames(x$xi) <- paste0("xi_", 1 : g)
  # rownames(x$xi) <- paste0("q_",  1 : q)
  # print(x$xi)
  # cat("\nFactor Covariance Matrices\n")
  # colnames(x$omega) <- paste0("q_",  1 : q)
  # rownames(x$omega) <- paste0("q_",  1 : q)
  # for(j in 1 : g) {
  #   cat(paste0("omega_", j, ":"), "\n")
  #   print(x$omega[,, j])
  #   cat("\n")
  # }
  # cat("Diagonal of the Error Covariance Matrix\n")
  # cat("diag D: \n", diag(x$D), "\n")
}

print.mctfa <- function(x, ...) {

  class(x) <- c("emmix", class(x))
  print(x)
  
  #   g <- x$g
  #   q <- x$q
  #   cat("Call:\n")
  #   print(x$call)
  #   cat("\n", paste0("Coefficients: (for i in 1 to ", g, ")", "\n"))
  #   cat("\nMixing Proportions\n")
  #   cat("pi_i : ", round(x$pivec, 3), "\n")
  #   cat("\nLoading Matrix\n")
  #   colnames(x$A) <- paste0("q_",  1 : q)
  #   cat("A: \n")
  #   print(x$A)
  #   cat("\nFactor Means (in columns)\n")
  #   colnames(x$xi) <- paste0("xi_", 1 : g)
  #   rownames(x$xi) <- paste0("q_",  1 : q)
  #   print(x$xi)
  #   cat("\nFactor Covariance Matrices\n")
  #   colnames(x$omega) <- paste0("q_",  1 : q)
  #   rownames(x$omega) <- paste0("q_",  1 : q)
  #   for(j in 1 : g) {
  #     cat(paste0("omega_", j, ":"), "\n")
  #     print(x$omega[,, j])
  #     cat("\n")
  #   }
  #   cat("Diagonal of the Error Covariance Matrix\n")
  #   cat("diag D: \n", diag(x$D), "\n")
  #   cat("nu: \n", x$v, "\n")
}

print.mfa <- function(x, ...) {
  class(x) <- c("emmix", class(x))
  print(x)
}

print.mtfa <- function(x, ...) {
  class(x) <- c("emmix", class(x))
  print(x)
}
