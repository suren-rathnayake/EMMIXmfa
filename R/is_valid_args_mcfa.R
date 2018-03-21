is_valid_args.mcfa <- function (Y, g, q, itmax, nkmeans, nrandom,
																tol, initClust, init_para, init_method,
																conv_measure, warn_messages) {

for (initm in 1 : length(init_method)) {
  
  if(any(tolower(init_method[initm]) == c('e', 'eigen', 'eigena', 'eigen-a')))
    init_method[initm] <- 'eigen-A'

  if(any(tolower(init_method[initm]) == c('r', 'rand', 'randa', 'rand-a')))
    init_method[initm] <- 'rand-A'

  #if(any(tolower(init_method) == c('m', 'mf', 'mfa')))
  #  init_method <- 'mfa'

  if(any(tolower(init_method[initm]) == c('g', 'gm', 'gmf')))
    init_method[initm] <- 'gmf'
}

if(any(tolower(conv_measure) == c('d', 'diff')))
  conv_measure <- 'diff'

if(any(tolower(conv_measure) == c('r', 'ratio')))
  conv_measure <- 'ratio'

#
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
if (is.null(p)) {
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

if ( (!is.null(init_method)) && (init_method != "eigen-A") 
    && (init_method != "rand-A") && (init_method != "gmf")) {

  ERR <- "init_method is must 'eigen-A', 'rand-A', 'gmf', or NULL."
  class(ERR) <- "error"
  return(ERR)
}

if ((conv_measure != "diff") && (conv_measure != "ratio")) {
  
  ERR <- "conv_measure needs to be either 'diff' or 'ratio'."
  class(ERR) <- "error"
	return(ERR)
}

if (class(warn_messages) != "logical") {
  
  ERR <- "warn_messages must either be TRUE or FALSE'."
  class(ERR) <- "error"
  return(ERR)
}

ERR <- list(init_method = init_method, conv_measure = conv_measure)

return(ERR)
}

is_valid_args.mctfa <- function(Y, g, q, itmax, nkmeans, nrandom,
																tol, df_init, df_update, init_clust,
																init_para, init_method, conv_measure,
																warn_messages) {

#
ERR <- is_valid_args.mcfa(Y, g, q, itmax, nkmeans, nrandom, tol,
                          init_clust, init_para, init_method,
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

return(ERR)
}
