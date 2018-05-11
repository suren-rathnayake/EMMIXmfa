is.defined <- function (sym, env) {
  sym <- deparse(substitute(sym))
  exists(sym, env)
}

chol.inv <- function(x, ...){
  C <- chol(x)
  inv_x <- chol2inv(C)
  return(inv_x)
}

perm <- function (n, r, v = 1 : n){
  if (r == 1)
    matrix(v, n, 1)
  else if (n == 1)
    matrix(v, 1, r)
  else {
    X <- NULL
    for (i in 1 : n)
      X <- rbind(X, cbind(v[i], perm(n - 1, r - 1, v[-i])))
    X
  }
}

# the minimum number of misallocations
minmis <- function(cls, hat_cls) {

  if (length(cls) != length(hat_cls))
    stop('Length of the two arguments should be equal')

  tcls <- rep(0, length(cls))
  labs <- unique(cls)
  for (j in 1 : length(labs))
    tcls[cls == labs[j]] <- j
  cls <- tcls

  tcls <- rep(0, length(hat_cls))
  labs <- unique(hat_cls)
  for(j in 1 : length(labs))
    tcls[hat_cls == labs[j]] <- j
  hat_cls <- tcls

  labs <- unique(c(hat_cls, cls))
  g <- length(labs)
  x <- perm(g, g)
  min_err <- Inf
  for (j in 1 : nrow(x)) {

    new_cls <- rep(0, length(cls))
    for (k in 1 : g)
      new_cls[cls == labs[k]] <- x[j, k]

    e <- sum(new_cls != hat_cls)
    if (e < min_err)
      min_err <- e
  }
  return(min_err)
}

ari <- function (cls, hat_cls) {

  if (length(cls) != length(hat_cls))
    stop('Length of the two arguments should be equal')

  tab <- table(cls, hat_cls)
  if (sum(diag(tab)) == length(cls))
    return(1)

  A <- sum(choose(tab, 2))
  B <- sum(choose(rowSums(tab), 2))
  C <- sum(choose(colSums(tab), 2))
  D <- choose(sum(tab), 2)

  ARI <- (A - B * C / D) / ( (B + C) / 2 - B * C / D)
  return(ARI)
}

as.num.fac <- function(x) {as.numeric(as.factor(x))}

plot.emmix <- function(x, ...) {

  if (x$q == 1) {

    plot(x$Umean, 1 : length(x$Umean),  xlim = range(x$Umean), axes = FALSE,
          xlab = expression(widehat(u)[1]), ylab = "", type = "p",
          pch = if (x$g <= 5) {20 + as.numeric(x$clust)} else {
          as.numeric(x$clust)}, col = as.numeric(x$clust),
          bg =  as.numeric(x$clust))
    axis(side = 1)
  }

  if (x$q == 2)

    plot(x$Umean[, c(1, 2)], col = as.numeric(x$clust),
          ylim = range(x$Umean[, 2]), pch = if (x$g <= 5) {
          20 + as.numeric(x$clust)} else{as.numeric(x$clust)},
          bg = as.numeric(x$clust), xlim = range(x$Umean[, 1]),
          xlab = expression(widehat(u)[1]),
          ylab=expression(widehat(u)[2]))

  if (x$q > 2)

    pairs(x$Umean, col = as.numeric(x$clust), bg = as.numeric(x$clust),
            pch = if(x$g <= 5) {
              20 + as.numeric(x$clust)} else{as.numeric(x$clust)})
}

summary.emmix <- function(object, ...) {
  
  cat("Call:\n")
  print(object$call)
  summ <- cbind(num_comp = object$g,
                num_fac = object$q,
                log_like = object$logL,
                BIC = object$BIC)
  cat("\n")
  print(summ)
}

summary.mfa <- function(object, ...) {
  class(object) <- c("emmix", class(object))
  summary(object)
}

summary.mcfa <- function(object, ...) {
  class(object) <- c("emmix", class(object))
  summary(object)
}

summary.mtfa <- function(object, ...) {
  class(object) <- c("emmix", class(object))
  summary(object)
}

summary.mctfa <- function(object, ...) {
  class(object) <- c("emmix", class(object))
  summary(object)
}
