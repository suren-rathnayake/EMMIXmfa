predict.emmix <- function(object, Y, ...) {
  
  if (missing(Y)) {
   stop("Missing newdata Y.")
  }
  
  if (any(class(object) == "mcfa")) {
     tau <- do.call("tau.mcfa", c(list(Y = Y), object))
  }
  if (any(class(object) == "mctfa")) {
    tau <- do.call("tau.mctfa", c(list(Y = Y), object))
  }
  if (any((class(object) == "mfa"))) {
    tau <- do.call("tau.mfa", c(list(Y = Y), object))
  }
  if (any((class(object) == "mtfa"))) {
    tau <- do.call("tau.mtfa", c(list(Y = Y), object))$tau
  }
  if (!( any(class(object) == "mcfa") || any(class(object) == "mctfa") ||
          any(class(object) == "mfa") || any(class(object) == "mtfa"))) {

    stop("STOP:object must be of class mcfa, mctfa, mfa or mtfa")
  } else {
    clust <- apply(tau, 1, which.max)
    return(clust)
  }
}

predict.mcfa <- function(object, Y, ...) {
  
  tau <- do.call("tau.mcfa", c(list(Y = Y), object))
  clust <- apply(tau, 1, which.max) 
  clust
}

predict.mctfa <- function(object, Y, ...) {
  
  tau <- do.call("tau.mctfa", c(list(Y = Y), object))
  clust <- apply(tau, 1, which.max) 
  clust
}

predict.mfa <- function(object, Y, ...) {

  tau <- do.call("tau.mfa", c(list(Y = Y), object))
  clust <- apply(tau, 1, which.max)
  clust
}

predict.mtfa <- function(object, Y, ...) {

  tau <- do.call("tau.mtfa", c(list(Y = Y), object))
  clust <- apply(tau, 1, which.max)
  clust
}