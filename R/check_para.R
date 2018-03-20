check_para <- function(p, q, g, init_para, model = "mcfa") {

if (is.null(init_para$pivec) || is.null(init_para$D))
  return (check = FALSE)

if ((model == "mcfa") || (model == "mctfa")) {

  if (!all(c("pivec", "A", "xi", "omega", "D") %in% names(init_para))) {
    message("Missing paramters in init_para")
    return (check = FALSE)
  }

  if (is.null(init_para$pivec) || is.null(init_para$A) ||
        is.null(init_para$xi) || is.null(init_para$omega) ||
        is.null(init_para$D)) {

   message("One or more paramters in init_para is NULL")
   return (check = FALSE)
  }

  if ((length(init_para$pivec) != g) || (ncol(init_para$xi) != g) ||
      (dim(init_para$omega)[3] != g)) {

    message("Paramters for one or more compoents not found in init_para")
    return(check = FALSE)
  }

  if ((ncol(init_para$A) != q) || (nrow(init_para$xi) != q) ||
      (dim(init_para$omega)[1] != q) || (dim(init_para$omega)[2] != q)) {

        message("Paramters for one or more factors not found in init_para")
        return(check = FALSE)
  }

  if ((nrow(init_para$A) != p) || (nrow(init_para$D) != p) ||
      (ncol(init_para$D) != p)) {

        message("The number of variables and init_para are not compatible")
        return(check = FALSE)
  }
}

if((model=="mfa") || (model=="mtfa")) {

  if (!all(c("pivec", "B", "mu", "D") %in% names(init_para))) {

    message("Missing paramters in init_para")
    return (check = FALSE)
  }

  if (is.null(init_para$pivec) || is.null(init_para$B) ||
      is.null(init_para$mu) || is.null(init_para$D)) {

    message("One or more paramters in init_para is NULL")
    return(check = FALSE)
  }

  if ((length(init_para$pivec) != g) || (dim(init_para$mu)[2] != g)) {

    message("Paramters for one or more compoents not found in init_para")
    return(check = FALSE)
  }

  if ((ncol(init_para$B) != q)) {

    message("Number of columns in B in init_para must be equal to q")
    return(check = FALSE)
  }

  if ((nrow(init_para$B) != p) || (nrow(init_para$D) != p) ||
     (ncol(init_para$D) != p)) {

    message("The number of variables and init_para are not compatible")
    return(check = FALSE)
  }
}

return(check = TRUE)
}

start_clust <- function(Y, g, init_clust, nkmeans, nrandom) {

n <- nrow(Y)
if (!is.null(init_clust)) {

  if ((class(init_clust) == "factor") || (class(init_clust) == "numeric")) {
    if (n != length(init_clust))
      stop("Length of init_clust must be equal to number of
                                         samples in the data")

    init_clust <- as.num.fac(init_clust)
    init_clust <- matrix(init_clust, nrow=n, ncol=1)
  }

  if (class(init_clust) == "matrix") {
    if (n != dim(init_clust)[1])
      stop("Length of init_clust must be equal to number of
                                           samples in the data")

    init_clust <- apply(init_clust, 2, as.num.fac)
    init_clust <- init_clust
  }
}

k_starts <- NULL
if (!is.null(nkmeans) && is.numeric(nkmeans) && (nkmeans > 0)) {

  k_starts <- matrix(NA, nrow = n, ncol = nkmeans)
  for (i in 1 : nkmeans)
    k_starts[, i] <- kmeans(Y, g)$cluster
}

r_starts <- NULL
if (!is.null(nrandom) && is.numeric(nrandom) &&  (nrandom > 0))
  r_starts <- matrix(sample(1 : g, n * nrandom, replace = TRUE),
                                        nrow = n, ncol = nrandom)
  
  nc <- ifelse(is.null(init_clust), 0, dim(init_clust)[2]) +
	        ifelse(is.null(k_starts), 0, dim(k_starts)[2]) +
	        ifelse(is.null(r_starts), 0, dim(r_starts)[2])

  starts <- matrix(NA, nrow = n, ncol = nc)
  ind <- 1
  if (!is.null(dim(init_clust))) {
    starts[, ind : dim(init_clust)[2]] <- init_clust
    ind <- ind + dim(init_clust)[2]
  }

  if (!is.null(dim(k_starts))) {
    starts[, ind : (ind - 1 + dim(k_starts)[2])] <- k_starts
    ind <- ind + dim(k_starts)[2]
  }

  if (!is.null(dim(r_starts))) {
    starts[, ind : (ind - 1 + dim(r_starts)[2])] <- r_starts
    ind <- ind + dim(r_starts)[2]
  }

  starts <- as.matrix(starts)
  if ((is.null(starts)) || (dim(starts)[2] == 0))
    stop("Incorrect specification of initial grouping parameters")

    if (dim(starts)[1] != n)
    stop("Length of init_clust must be equal to number of
                                             samples in the data")

  return(starts)
}
