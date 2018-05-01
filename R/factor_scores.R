# factor_scores <- function(Y, model, tau = NULL, clust= NULL, ...) {

#   if (class(Y) == "data.frame") {
#    Y <- as.matrix(Y)
#   }
#   if (class(Y) != "matrix")
#     stop("Y needs to be a numeric matrix")

#   if (!is.null(tau)) {
#     model$tau <- tau
#   }
#   if (!is.null(clust)) {
#     model$clust <- clust
#   }

#   if (any(class(model) == "mcfa")) {
#      scores <- do.call("factor_scores_mcfa", c(list(Y = Y), model))
#   }
#   if (any(class(model) == "mctfa")) {
#     scores <- do.call("factor_scores_mctfa", c(list(Y = Y), model))
#   }
#   if (any((class(model) == "mfa"))) {
#     scores <- do.call("factor_scores_mfa", c(list(Y = Y), model))
#   }
#   if (any((class(model) == "mtfa"))) {
#     scores <- do.call("factor_scores_mtfa", c(list(Y = Y), model))
#   }

# return(scores)
# }

factor_scores <- function(model, Y, ...) UseMethod("factor_scores")

factor_scores.mcfa <- function(model, Y, tau = NULL, clust= NULL, ...) {
  
  if (class(Y) == "data.frame") {
    Y <- as.matrix(Y)
  }
  
  if (class(Y) != "matrix")
    stop("Y needs to be a numeric matrix")
  
  if (!is.null(tau)) {
    model$tau <- tau
  }
  
  if (!is.null(clust)) {
    model$clust <- clust
  }
  
  scores <- do.call("factor_scores_mcfa", c(list(Y = Y), model))
  scores
}

factor_scores.mctfa <- function(model, Y, tau = NULL, clust= NULL, ...) {
  
  if (class(Y) == "data.frame") {
    Y <- as.matrix(Y)
  }
  
  if (class(Y) != "matrix")
    stop("Y needs to be a numeric matrix")
  
  if (!is.null(tau)) {
    model$tau <- tau
  }
  
  if (!is.null(clust)) {
    model$clust <- clust
  }
  
  scores <- do.call("factor_scores_mctfa", c(list(Y = Y), model))
  scores
}

factor_scores.mfa <- function(model, Y, tau = NULL, clust= NULL, ...) {
  
  if (class(Y) == "data.frame") {
    Y <- as.matrix(Y)
  }
  
  if (class(Y) != "matrix")
    stop("Y needs to be a numeric matrix")
  
  if (!is.null(tau)) {
    model$tau <- tau
  }
  
  if (!is.null(clust)) {
    model$clust <- clust
  }
  
  scores <- do.call("factor_scores_mfa", c(list(Y = Y), model))
  scores
}

factor_scores.mtfa <- function(model, Y, tau = NULL, clust= NULL, ...) {
  
  if (class(Y) == "data.frame") {
    Y <- as.matrix(Y)
  }
  
  if (class(Y) != "matrix")
    stop("Y needs to be a numeric matrix")
  
  if (!is.null(tau)) {
    model$tau <- tau
  }
  
  if (!is.null(clust)) {
    model$clust <- clust
  }
  
  scores <- do.call("factor_scores_mtfa", c(list(Y = Y), model))
  scores
}
