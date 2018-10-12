mcfa <- function (Y, g, q, itmax = 500, nkmeans = 5, nrandom = 20,
                          tol = 1.e-5,  init_clust = NULL, init_para = NULL,
                          init_method = NULL, conv_measure = 'diff',
                          warn_messages = TRUE, ...) {

if (!is.matrix(Y))
  Y <- as.matrix(Y)

# check and fix arguments
ERR <- is_valid_args.mcfa(Y, g, q, itmax, nkmeans, nrandom, tol,
                          init_clust, init_para, init_method,
                          conv_measure, warn_messages)
# return error is the arguments are mis-specified
if (class(ERR) == "error") {

  stop(unclass(ERR))
} else {

  init_method <- ERR$init_method
  conv_measure <- ERR$conv_measure
}

p <- ncol(Y)
n <- nrow(Y)
warn_msg <- NULL
maxLOGL <- -Inf

pb <- txtProgressBar(style = 3, char = '.')
prog <- 0

# Fit a model using given initial parameter estimates
if (!is.null(init_para)) {

  prog = 1/(1 + nkmeans + nrandom)
  setTxtProgressBar(pb, prog)

  if (!check_para(p, q, g, init_para, "mcfa"))
    stop("Specified init_para is not acceptable", call. = FALSE)

  init_para <- init_para[c("g", "q", "pivec", "A", "xi",
                                          "omega", "D")]

  estd_model <- est.mcfa(init_para = init_para, Y = Y, itmax = itmax,
                         tol = tol, conv_measure = conv_measure)

  if ((class(estd_model) == "mcfa")) {
    if (estd_model$logL > maxLOGL) {
      Hmodel <- estd_model
      maxLOGL <- Hmodel$logL
    }
    # cat(sprintf("\n g = %i, q = %i, init_para,  logL %8.4f \n",
    #                 g, q,  estd_model$logL))
  }

  if (class(estd_model) == "error") {
    when <- paste("init_para")
    what <- estd_model
    warn_msg <- cbind(when, what)
    colnames(warn_msg) <- c('when', 'what')
  }
}

# Use k-means, random, and given groupings of observations
# for parameter estimates.
if ((nkmeans!=0) || (nrandom!=0) || (!is.null(init_clust))) {

  # The set of all initial partitions
  initial_partitions <- start_clust(Y, g, init_clust, nkmeans, nrandom)

  maxinit <- ncol(initial_partitions)
  if (is.null(maxinit))
    maxinit <- 1
  # for progress bar
  if (prog == 0) {
    prog <- 1/maxinit; tinit <- maxinit
  } else {
    prog  <- prog + 1/maxinit; tinit <- 1 + maxinit
  }

  # For each initial starts run EM
  for (ii in 1 : maxinit) {

    # Avoid cases where there is only one sample is assigned to a group
    if (min(table(initial_partitions[, ii])) == 1) {
      when <- paste("At start", ii)
      what <- "Initial partitioning has a single sample cluster."
      warn_msg <- rbind(warn_msg, cbind(when, what))
      next
    }

    method_for_init_para <- init_method
    if (is.null (init_method))
      method_for_init_para <- c("rand-A", "eigen-A", "gmf")

    # If there are duplicate initial partitions, don't run eigen-A
    if (ii <= maxinit - 1) {
      dup_col <- apply(initial_partitions[, (ii + 1) : maxinit, drop = FALSE],
                       2, function(x) ari(x, initial_partitions[, ii]))

      if (any(dup_col == 1))
        method_for_init_para <-
          method_for_init_para[!method_for_init_para  %in% "eigen-A"]
    }

    for (method_for_init in method_for_init_para) {

      # Initial estimates of model parameters
      init_model_para <- try(init_est_para.mcfa(Y, g, q,
                                          initial_partitions[, ii],
                                          method_for_init), silent = TRUE)

      if (class(init_model_para) == "try-error") {
        when <- paste("At start", ii)
        what <- "Failed to estimate initial parameters"
        warn_msg <- rbind(warn_msg, cbind(when, what))
        next
      }

      # EM steps
      estd_model <- est.mcfa(init_para = init_model_para, Y = Y,
  		                    	   itmax = itmax, tol = tol,
                               conv_measure = conv_measure)

      # keep the model with highest log-likelihood
      if ((class(estd_model) == "mcfa")) {

        if (estd_model$logL > maxLOGL) {
          Hmodel <- estd_model
          maxLOGL <- Hmodel$logL
        }
        # cat(sprintf("\n g = %i, q = %i, init %i logL %8.4f,
        #               maxlogL = %8.4f \n",
        #               g, q, ii, estd_model$logL, maxLOGL))
      }

      if (class(estd_model) == "error") {

        when <- paste("At start", ii)
        what <- estd_model
        warn_msg <- rbind(warn_msg, cbind(when, what))
      }
    }

    setTxtProgressBar(pb, prog)
    prog <- prog + 1/tinit
  }
}
setTxtProgressBar(pb, 1)
close(pb)

if (!exists("Hmodel")) {
  cat("Failed to Estimate a Model. See Error Messages.")
  return(warn_msg)
}

# Make A^T A = I_q
CH <- chol(t(Hmodel$A) %*% Hmodel$A)
Hmodel$A <- Hmodel$A %*% solve(CH)
Hmodel$xi <- CH %*% Hmodel$xi
for (i in 1 : g) {
  Hmodel$omega[,,i] <- CH %*% Hmodel$omega[,,i] %*% t(CH)
}
d <- (g - 1) + p + q * (p + g) + g * q * (q + 1) / 2 - q * q
Hmodel$tau <- do.call('tau.mcfa', c(list(Y = Y), Hmodel))
Hmodel$BIC <- -2 * Hmodel$logL + d * log(n)
Hmodel$clust <- apply(Hmodel$tau, 1, which.max)
Hmodel <- append(Hmodel,
		do.call('factor_scores_mcfa', c(list(Y = Y), Hmodel)))
Hmodel$call <- match.call()
if (is.null(warn_msg)) warn_msg <- "no warnings"
if (warn_messages == TRUE)
   Hmodel <- c(Hmodel, list(warn_msg = warn_msg))
class(Hmodel) <- c("emmix", "mcfa")
return(Hmodel)
}
