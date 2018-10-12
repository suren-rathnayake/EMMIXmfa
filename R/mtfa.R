mtfa <- function (Y, g, q, itmax = 500,
                         nkmeans = 20, nrandom = 20,
                         tol = 1.e-5, df_init = rep(30, g),
                         df_update = TRUE, sigma_type = 'common',
                         D_type = 'common', init_clust = NULL,
                         init_para = NULL, conv_measure = 'diff',
                         warn_messages = TRUE, ...)  {

if (!is.matrix(Y))
  Y <- as.matrix(Y)

p <- ncol(Y)
n <- nrow(Y)

if (any(tolower(sigma_type) == c('c', 'common')))
  sigma_type = 'common'
if (any(tolower(sigma_type) == c('u', 'unique')))
  sigma_type = 'unique'
if (any(tolower(D_type) == c('c', 'common')))
  D_type = 'common'
if (any(tolower(D_type) == c('u', 'unique')))
  D_type = 'unique'
if (any(tolower(conv_measure) == c('d', 'diff')))
  conv_measure = 'diff'
if (any(tolower(conv_measure) == c('r', 'ratio')))
  conv_measure = 'ratio'

ERR <- is_valid_args.mtfa(Y, g, q, itmax, nkmeans, nrandom, tol, df_init,
                          df_update, sigma_type, D_type, init_clust,
                          init_para, conv_measure, warn_messages)

# if (class(ERR) == "error")
#   stop(cat(ERR))
if (class(ERR) == "error") {

  stop(unclass(ERR))
}

pb <- txtProgressBar(style=3, char='.')
prog <- 0
warn_msg <- NULL
max_log_like  <- -Inf
if (!is.null(init_para)) {

  prog <- 1 / (1 + nkmeans + nrandom)
  setTxtProgressBar(pb, prog)

  if (!check_para(p, q, g, init_para, "mtfa"))
    stop("Incorrect specification of init_para", call. = FALSE)

  v <- init_para$v
  init_para <- init_para[c("g", "q", "pivec", "mu", "B", "D",
                             "sigma_type", "D_type")]

  estd_model <- est.mtfa(init_para = init_para, Y = Y,
                           itmax = itmax, tol = tol,
                           v = v, df_update = df_update,
                           conv_measure = conv_measure)

  if (class(estd_model) == "mtfa") {

    if (estd_model$logL > max_log_like) {
        Hmodel <- estd_model
        max_log_like <- Hmodel$logL
    }
    # cat(sprintf("\n g = %i, q = %i, init_para,  logL %8.4f \n",
    #                   g, q,  estd_model$logL))
  }

  if (class(estd_model) == "error") {
    when <- paste("init_para")
    what <- estd_model
    warn_msg <- cbind(when, what)
    colnames(warn_msg) <- c('when', 'what')
  }
}

if ((nkmeans!=0) || (nrandom!=0) || (!is.null(init_clust)))  {

  initial_partitions <- start_clust(Y, g, init_clust, nkmeans, nrandom)
	maxinit <- ncol(initial_partitions)
	if (is.null(maxinit)) maxinit <- 1

    if (prog == 0) {
      prog <- 1 / maxinit;  tinit <- maxinit
    } else {
      prog <- prog + 1/maxinit; tinit <- 1 + maxinit
    }

  for (ii in 1:maxinit) {

    if (min(table(initial_partitions[, ii])) == 1) {
      when <- paste("At start", ii)
      what <- "Initial partitioning has a single sample cluster."
      warn_msg <- rbind(warn_msg, cbind(when, what))
      next
    }

    # Initial parameters estimate are same as mfa
    init_model_para <- try(init_est_para.mfa(Y, g, q,
                                                initial_partitions[, ii],
                                                sigma_type, D_type),
                                                silent = TRUE)
    if (class(init_model_para) == "try-error") {
      when <- paste("At start", ii)
      what <- "Failed to estimate initial parameters"
      warn_msg <- rbind(warn_msg, cbind(when, what))
      next
    }

    class(init_model_para) <- "mtfa"

    estd_model <- est.mtfa(init_para = init_model_para, Y = Y,
                           itmax = itmax, tol = tol,
                           v = df_init, df_update = df_update,
                           conv_measure = conv_measure)

    if((class(estd_model) == "mtfa")) {
      if (estd_model$logL > max_log_like){
        Hmodel <- estd_model
        max_log_like <- Hmodel$logL
      }
  # 	 cat(sprintf("\n g = %i, q = %i, init %i logL %8.4f,
  #                 max_log_like = %8.4f \n",
  #                 g, q, ii, estd_model$logL, max_log_like))
    }

    if(class(estd_model) == "error") {
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

if (sigma_type=='common') {
  d <- (g - 1) + 2 * g * p + (p * q - q * (q - 1) / 2)
  if (df_update)
    d <- d + g
} else {
  if (D_type == 'common') {
    d <- (g - 1) + 2 * g * p + g * (p * q - q * (q - 1) / 2)
    if (df_update)
      d <- d + g
  } else {
    d <- (g - 1) + g * p + p + g * (p * q - q * (q - 1) / 2)
    if (df_update)
      d <- d + g
  }
}
Hmodel$BIC <- -2 * Hmodel$logL + d * log(n)
Hmodel$clust <- apply(Hmodel$tau, 1, which.max)
Hmodel <- append(Hmodel, do.call('factor_scores_mtfa',
                  c(list(Y = Y), Hmodel)))
Hmodel$call <- match.call()
if (warn_messages == TRUE)
   Hmodel <- c(Hmodel, list(warn_msg = warn_msg))
Hmodel["W"] <- NULL
class(Hmodel) <- c("emmix", "mtfa")
return(Hmodel)
}
