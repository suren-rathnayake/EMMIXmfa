# est.mtfa <- function(Y, g, q, itmax, tol, pivec, B, mu, D, sigma_type,
# 						D_type, v, df_update, conv_measure, ...) {
est.mtfa <- function (init_para, Y, itmax, tol,
                                        v, df_update, conv_measure, ...) {

p <- ncol(Y)
n <- nrow(Y)
fit <-  c(init_para, list(v = v, df_update = df_update))

# fit <- list(g = g, q = q, pivec = pivec, B = B, mu = mu, D = D,
# 						sigma_type = sigma_type, D_type = D_type, v = v,
# 						df_update = df_update)

# loglikeNtau <- try(do.call('logL_tau.mtfa', c(list(Y=Y), fit)),
                                    # silent = TRUE)

 loglikeNtau <- do.call('logL_tau.mtfa', c(list(Y=Y), fit))

if ((class(loglikeNtau) == "try-error")||
              (class(loglikeNtau) == 'character')) {
	FIT <- paste('in computing the log-likelihood before EM-steps')
    class(FIT) <- "error"
    return(FIT)
}
if (class(loglikeNtau$logL) == 'character') {
	FIT <- paste('in computing the log-likelihood before EM-steps',
                   loglikeNtau$logL)
    class(FIT) <- "error"
    return(FIT)
}

fit <- append(fit, loglikeNtau)

for (niter in 1 : itmax) {

  FIT <- do.call('Mstep.mtfa', c(list(Y=Y), fit))

   if (class(FIT) == 'error') {
        FIT <- paste('in ', niter,
                   'iteration of the M-step', FIT)
        class(FIT) <- "error"
        return(FIT)
   }

	loglikeNtau <- try(do.call('logL_tau.mtfa', c(list(Y=Y), FIT)),
        											silent=TRUE)

    if ((class(loglikeNtau) == "try-error") ||
                    (class(loglikeNtau) == 'character')) {

        FIT <- paste('in computing the log-likelihood after the ', niter,
                   'th the M-step', FIT$logL, sep='')
        class(FIT) <- "error"
        return(FIT)
    }

	FIT <- append(FIT, loglikeNtau)

	if ((class(FIT$logL)=="NULL") || (class(FIT$logL) == 'character')) {

		FIT <- paste('in computing the log-likelihood after the ', niter,
                     'th the M-step', FIT$logL, sep='')
        class(FIT) <- "error"
        return(FIT)
    } else {
        if ((FIT$logL == -Inf) || is.na(FIT$logL)) {

            FIT <- paste('the log-likelihood computed after the ', niter,
                       'th iteration of the M-step is not finite', sep='')
            class(FIT) <- "error"
            return(FIT)
        }
	}

	if ((conv_measure == "diff") && (abs(FIT$logL-fit$logL) < tol))
        break

	if ((conv_measure == "ratio") && (abs((FIT$logL-fit$logL) / FIT$logL) < tol))
        break

	fit <- FIT
}

class(FIT) <- "mtfa"
return(FIT)
}
