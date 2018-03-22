is_valid_args.mfa <- function(Y, g, q, itmax, nkmeans, nrandom, tol,
														  sigma_type, D_type, init_clust, init_para,
														  conv_measure, warn_messages) {


if (itmax < 1) {
  ERR <- "Maximum number of iterations, itmax, must be greather than one."
  class(ERR) <- "error"
  return(ERR)
}

if (any(is.na(Y))) {
  ERR <- "`Y' has missing value."
  class(ERR) <- "error"
  return(ERR)
}

if (!any(is.numeric(Y))) {
  ERR <- "`Y' has a non-numeric element."
  class(ERR) <- "error"
  return(ERR)
}

p <- ncol(Y)
if ((is.null(p)) || (p == 1)) {
  ERR <- "The data must have more than one variable."
  class(ERR) <- "error"
  return(ERR)
}

if (p <= q) {
  ERR <- "The number of factors must be less than the number of variables."
  class(ERR) <- "error"
  return(ERR)
}

if (q < 0) {
  ERR <- "q must be a positive integer."
  class(ERR) <- "error"
  return(ERR)
}

if (q != round(q)) {
  ERR <- "q must be a positive integer."
  class(ERR) <- "error"
  return(ERR)
}

if (g < 0) {
  ERR <- "g must be a positive integer."
  class(ERR) <- "error"
  return(ERR)
}

if (g != round(g)) {
  ERR <- "g must be a positive integer."
  class(ERR) <- "error"
  return(ERR)
}

if (nkmeans < 0) {
  ERR <- "nkmeans must be a positive integer."
  class(ERR) <- "error"
  return(ERR)
}

if (nrandom < 0) {
  ERR <- "nrandom must be a positive integer."
  class(ERR) <- "error"
  return(ERR)
}

if (tol < 0) {
  ERR <- "tol must be a greater than zero."
  class(ERR) <- "error"
  return(ERR)
}

if ((sigma_type == 'common') && (D_type == 'unique')) {
  ERR <- "D_type = 'unique' not available with sigma_type = 'common'."
  class(ERR) <- "error"
  return(ERR)
}

if ((sigma_type != "common") && (sigma_type != "unique")) {
  ERR <- "sigma_type needs to be either 'unique' or 'common'."
  class(ERR) <- "error"
  return(ERR)
}

if ((D_type != "common") && (D_type != "unique")) {
  ERR <- "D_type needs to be either 'unique' or 'common'."
  class(ERR) <- "error"
  return(ERR)
}

if ((conv_measure != "diff") && (conv_measure != "ratio")) {
  
  ERR <- "conv_measure needs to be either 'diff' or 'ratio'."
  class(ERR) <- "error"
  return(ERR)
}

if (class(warn_messages) != "logical") {
  
  ERR <- "warn_messages must either be TRUE or FALSE."
  class(ERR) <- "error"
  return(ERR)
}

if (!is.null(init_para)) {
  if (init_para$sigma_type != sigma_type) {
    ERR <- "`init_para$sigma_type` is not same as `sigma_type`."
    class(ERR) <- "error"
    return(ERR)
  }

  if (init_para$D_type != D_type) {
    ERR <- "`init_para$D_type` is not same as `D_type`."
    class(ERR) <- "error"
    return(ERR)
  }
}
ERR <- "TRUE"
return(ERR)
}

###
is_valid_args.mtfa <- function(Y, g, q, itmax, nkmeans, nrandom,
	tol, df_init, df_update, sigma_type, D_type, init_clust, init_para,
  conv_measure, warn_messages) {


ERR <- is_valid_args.mfa(Y, g, q, itmax, nkmeans, nrandom, tol,
												 sigma_type, D_type, init_clust, init_para,
												 conv_measure, warn_messages)
if (class(ERR) == "error")
		return(ERR)

if (class(df_update) != "logical") {
  ERR <- "df_update is either TRUE or FALSE'."
  class(ERR) <- "error"
  return(ERR)
}

if (!any(is.numeric(df_init))) {
  ERR <- "df_init has a non-numeric element."
  class(ERR) <- "error"
  return(ERR)
}

ERR <- "TRUE"
return(ERR)
}
