est.mctfa <- function(init_para, Y, itmax, tol, conv_measure, ...) {

p <- ncol(Y)
n <- nrow(Y)

loglike_and_tau <- try(do.call('logL_tau.mctfa',
                            c(list(Y = Y), init_para)), silent=TRUE)

if ((class(loglike_and_tau) == "try-error") ||
          (class(loglike_and_tau) == 'character')) {
  FIT <- paste('in computing the log-likelihood before EM-steps')
  class(FIT) <- "error"
  return(FIT)
}

# append the list  loglike_and_tau to init_para
init_para <- append(init_para, loglike_and_tau)

if (class(init_para$logL) == 'character') {

  FIT <- paste('in computing the log-likelihood before the EM-steps,',
                 init_para$logL)
  class(FIT) <- "error"
  return(FIT)
}

for (niter in 1 : itmax) {

  FIT <- do.call('Mstep.mctfa', c(list(Y=Y), init_para))
  if (class(FIT) == 'error') {
     FIT <- paste('in ', niter,
                   'iteration of the M-step:',
                   FIT)
     class(FIT) <- "error"
     return(FIT)
  }

  loglike_and_tau <- try(do.call('logL_tau.mctfa', c(list(Y = Y), FIT)),
                      silent = TRUE)

  if ((class(loglike_and_tau) == "try-error") ||
                    (class(loglike_and_tau) == 'character')) {

    FIT <- paste('in computing the log-likelihood after the ', niter,
                   'th the M-step:', FIT$logL, sep='')
    class(FIT) <- "error"
    return(FIT)
  }

  # append the list  loglike_and_tau 
  FIT <- append(FIT, loglike_and_tau)

  if ((class(FIT$logL)=="NULL")|| (class(FIT$logL) == 'character')) {

    FIT <- paste('in computing the log-likelihood after the ', niter,
                   'th the M-step:', FIT$logL, sep='')
    class(FIT) <- "error"
    return(FIT)
  } else {
    if ((FIT$logL == -Inf) || is.na(FIT$logL)) {

      FIT <- paste('Log-likelihood computed after the ', niter,
                       'th iteration of the M-step is not finite', sep='')
      class(FIT) <- "error"
      return(FIT)
    }
  }

  if ((conv_measure == "diff") && (abs(FIT$logL - init_para$logL) < tol))
     break

  if ((conv_measure == "ratio") &&
      (abs((FIT$logL - init_para$logL) / FIT$logL) < tol))
  break

init_para <- FIT
}

class(FIT) <- "mctfa"
return(FIT)
}
