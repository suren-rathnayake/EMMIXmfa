mfa <- function (Y, g, q, itmax = 500,
                        nkmeans = 20, nrandom = 20,
                        tol = 1.e-5, sigma_type = 'common',
                        D_type = 'common', init_clust = NULL,
                        init_para = NULL, conv_measure = 'diff',
                        warn_messages = TRUE, ...) {


if (!is.matrix(Y))
  Y <- as.matrix(Y)

p <- ncol(Y)
n <- nrow(Y)

if (any(tolower(sigma_type) == c('c', 'common')))
  sigma_type <- 'common'
if (any(tolower(sigma_type) == c('u', 'unique')))
  sigma_type <- 'unique'
if (any(tolower(D_type) == c('c', 'common')))
  D_type <- 'common'
if (any(tolower(D_type) == c('u', 'unique')))
  D_type <- 'unique'
if (any(tolower(conv_measure) == c('d', 'diff')))
  conv_measure <- 'diff'
if (any(tolower(conv_measure) == c('r', 'ratio')))
  conv_measure <- 'ratio'

ERR <- is_valid_args.mfa(Y, g, q, itmax, nkmeans, nrandom, tol,
                         sigma_type, D_type, init_clust, init_para,
                         conv_measure, warn_messages)

if (class(ERR) == "error") {

  stop(unclass(ERR))
}

# if (class(ERR) == "error") {

#   print(paste("Error:", ERR[1]))
#   return(ERR)
# }

pb <- txtProgressBar(style = 3, char = '.')
prog <- 0
warn_msg <- NULL
max_log_like  <- -Inf
if (!is.null(init_para)) {

  prog <- 1/(1 + nkmeans + nrandom)
  setTxtProgressBar(pb, prog)

  if (all(class(init_para) != "mfa"))
    class(init_para) <- "mfa"

  if (!check_para(p, q, g, init_para, "mfa"))
    stop("incorrect specification of init_para", call. = FALSE)

  init_para <- init_para[c("g", "q", "pivec", "mu", "B", "D",
                           "sigma_type", "D_type")]

  estd_model <- est.mfa(init_para = init_para, Y = Y, g = g, q = q,
                        itmax = itmax, tol = tol,
                        conv_measure = conv_measure)

  if ((class(estd_model) == "mfa")) {
    if (estd_model$logL > max_log_like) {
      Hmodel <- estd_model
      max_log_like <- Hmodel$logL
    }
    # cat(sprintf("\n g = %i, q = %i, init_para,  logL %8.4f \n",
    #                 g, q,  estd_model$logL))
  }
  if(class(estd_model) == "warn") {
    when <- paste("init_para")
    what <- estd_model
    warn_msg <- cbind(when, what)
    colnames(warn_msg) <- c('when', 'what')
  }
}

if ( (nkmeans!=0) || (nrandom!=0) || (!is.null(init_clust))) {

  initial_partitions <- start_clust(Y, g, init_clust, nkmeans, nrandom)
  maxinit <- ncol(initial_partitions)
  if(is.null(maxinit)) maxinit <- 1

  # For the progress bar
  if (prog == 0) {
    tinit <- maxinit
    prog  <- 1/maxinit;
  } else {
    tinit <-  1 + maxinit
    prog <-  prog + 1/tinit
  }

  for (ii in 1 : maxinit) {

    if (min(table(initial_partitions[, ii])) == 1) {
      when <- paste("At start", ii)
      what <- "Initial partitioning has a single sample cluster."
      warn_msg <- rbind(warn_msg, cbind(when, what))
      next
    }

    init_model_para <- try(init_est_para.mfa(Y, g, q,
                           initial_partitions[, ii],
                           sigma_type, D_type), silent = TRUE)

    if (class(init_model_para) == "try-error") {
      when <- paste("At start", ii)
      what <- "Failed to estimate initial parameters"
      warn_msg <- rbind(warn_msg, cbind(when, what))
      next
    }

    estd_model <- est.mfa(init_para = init_model_para, Y = Y, g = g, q = q,
                          itmax = itmax, tol = tol,
                          conv_measure = conv_measure)

    if ((class(estd_model) == "mfa")) {
      if (estd_model$logL > max_log_like) {
        Hmodel <- estd_model
        max_log_like <- Hmodel$logL
      }
    # cat(sprintf("\n g = %i, q = %i, init %i logL %8.4f,
    #     max_log_like = %8.4f \n",
    #     g, q, ii, estd_model$logL, max_log_like))
    }
    if (class(estd_model) == "warn") {
      when <- paste("At start", ii)
      what <- estd_model
      warn_msg <- rbind(warn_msg, cbind(when, what))
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

if (sigma_type == 'common') {

	d <- (g - 1) + 2 * g * p + (p * q - q * (q - 1) / 2)
} else {

  if (D_type == 'common') {
	  d <- (g - 1) + 2 * g * p + g * (p * q - q * (q - 1) / 2)
  } else {
	  d <- (g-1) + g * p + p + g * (p * q - q* (q - 1) / 2)
  }
}

Hmodel$BIC <- -2 * Hmodel$logL + d * log(n)
Hmodel$clust <- apply(Hmodel$tau, 1, which.max)
Hmodel <- append(Hmodel, do.call('factor_scores_mfa',
                c(list(Y = Y), Hmodel)))
Hmodel$call <- match.call()
if (warn_messages == TRUE)
   Hmodel <- c(Hmodel, list(warn_msg = warn_msg))

class(Hmodel) <- c("emmix", "mfa")
return(Hmodel)
}
