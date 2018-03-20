est.mcfa <- function(init_para, Y, itmax, tol, conv_measure, ...) {

p <- ncol(Y)
n <- nrow(Y)

init_para$logL <- do.call("loglike.mcfa", c(list(Y=Y), init_para))

if ((class(init_para$logL) == "try-error") ||
      (class(init_para$logL) == 'character')) {

  FIT <- paste('in computing the log-likelihood before EM-steps')
  class(FIT) <- "error"
  return(FIT)
}


for (niter in 1 : itmax) {

  FIT <- do.call('Mstep.mcfa', c(list(Y = Y), init_para))

  if (class(FIT) == 'error') {

    FIT <- paste('in ', niter,
                   'iteration of the M-step:',
                   FIT)
    class(FIT) <- "error"
    return(FIT)
  }

  FIT$logL <- try(do.call('loglike.mcfa', c(list(Y = Y), FIT)))

  if ((class(FIT$logL) == "try-error") ||
        (class(FIT$logL) == 'character')) {

    FIT <- paste('in computing the log-likelihood after the ', niter,
                   'th the M-step', FIT$logL, sep = '')
    class(FIT) <- "error"
    return(FIT)
  }

  if ((FIT$logL == -Inf) | is.na(FIT$logL)){

    FIT <- paste('Log likelihood computed after the', niter,
      'th iteration of the M-step is not finite', sep='')
    class(FIT) <- "error"
    return(FIT)
  }

  if ((conv_measure == "diff") && (abs(FIT$logL - init_para$logL) < tol))
    break

  if ((conv_measure == "ratio") && 
      (abs((FIT$logL - init_para$logL)/FIT$logL) < tol))
    break

  init_para <- FIT
}

class(FIT) <- "mcfa"
return(FIT)
}
