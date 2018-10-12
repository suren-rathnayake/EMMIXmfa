\name{factor_scores}
\alias{factor_scores}
\alias{factor_scores.mcfa}
\alias{factor_scores.mctfa}
\alias{plot.emmix}
\title{
 Computes Factor Scores
}
\description{
This function computes factor scores for observations. 
Using factor scores,
we can represent the original data point \eqn{y_j} in a 
\emph{q}-dimensional reduced space. This is only meaningful
in the case of \code{mcfa} or \code{mctfa} models,
as the factor cores for \code{mfa} and \code{mtfa} are
white noise.

The (estimated conditional expectation of) unobservable factors
\eqn{U_{ij}} given \eqn{y_j} and the component membership
can be expressed by,
\deqn{
\hat{u}_{ij} = E_{\hat{\Psi}}\{U_{ij} \mid y_j, z_{ij} = 1\}.
}

The estimated mean \eqn{U_{ij}} (over the
component membership of \eqn{y_j})
is give as
\deqn{
\hat{u}_{j} = \sum_{i=1}^g \tau_i(y_j; \hat{\Psi}) \hat{u}_{ij},
}
where \eqn{\tau_i(y_j; \hat{\Psi})}
estimated posterior probability of \eqn{y_j}
belonging to the \eqn{i}th component.

An alternative estimate of \eqn{u_j}, the posterior expectation
of the factor corresponding to the \emph{j}th observation \eqn{y_j}, is
defined by replacing \eqn{\tau_i(y_j;\,\hat{\Psi})} by \eqn{\hat{z}_{ij}},
where
\eqn{\hat{z}_{ij} = 1}, if \eqn{\hat{\tau}_i(y_j; \hat{\Psi})}
>= \eqn{\hat{\tau_h}(y_j; \hat{\Psi})
(h=1,\,\dots,\,g; h \neq i)}, else
\eqn{\hat{z}_{ij} = 0}.
\deqn{
\hat{u}_{j}^C = \sum_{i=1}^g \hat{z}_{ij} \hat{u}_{ij}.
}

For MFA, we have
\deqn{
\hat{u}_{ij} = \hat{\beta}_i^T (y_j - \hat{\mu}_i),
}
and
\deqn{
\hat{u}_{j} = \sum_{i=1}^g \tau_i(y_j;  \hat{\Psi}) \hat{\beta}_i^T
  (y_j - \hat{\mu}_i)
}
for \eqn{j = 1, \dots, n} where
\eqn{\hat{\beta}_i = (B_iB_i^T + D_i)^{-1} B_i}.

For MCFA,
\deqn{
\hat{u}_{ij} = \hat{\xi}_i + \hat{\gamma}_i^T (y_j -\hat{A}\hat{\xi}_i),
}
\deqn{
\hat{u}_{j} = \sum_{i=1}^g\tau_i(y_j; \hat{\Psi})
\{\hat{\xi}_i + \hat{\gamma}_i^T(y_j -\hat{A}\hat{\xi}_i)\},
}
where \eqn{\gamma_i =  (A \Omega_i A + D)^{-1} A \Omega_i}.

With M\emph{t}FA and MC\emph{t}FA, the distribution of
\eqn{\hat{u}_{ij}} and of \eqn{\hat{u}_{j}}
have the same form as those of MFA and MCFA, respectively.
}

\usage{
factor_scores(model, Y, ...)
\method{factor_scores}{mcfa}(model, Y, tau = NULL, clust= NULL, ...)
\method{factor_scores}{mctfa}(model, Y, tau = NULL, clust= NULL, ...)
\method{plot}{emmix}(x, ...)
}
\arguments{
  \item{model}{
   An object of class \code{mfa, mcfa, mtfa} or \code{mctfa}.
}
\item{x}{
   An object of class \code{mfa, mcfa, mtfa} or \code{mctfa}.
}
  \item{Y}{
  Data matrix with variables in columns in the same order as used in
  model estimation.
}
  \item{tau}{
  Optional. Posterior probabilities of belonging to the components
  in the mixture model. If not provided, they will be computed based on 
  the \code{model} parameters.
}
 \item{clust}{
 Optional. Indicators of belonging to the components.
 If not provided, will be estimated using \code{tau}.
}
 \item{\dots}{
   Not used.
}
}
\details{
Factor scores can be used in visualization of the data
in the factor space.
}
\value{
  \item{Uscores}{
    Estimated conditional expected component scores of the
    unobservable factors given the data and the component membership
    (\eqn{\hat{u}_{ij}}).
    Size is \eqn{n \times q \times g}, where \code{n} is the number of sample,
    \code{q} is the number of factors and \code{g} is the number components.
}
  \item{Umean}{
   Means of the estimated conditional expected factors scores over
   estimated posterior distributions (\eqn{\hat{u}_{j}}).
   Size \eqn{n \times q}.
}
 \item{Uclust}{
  Alternative estimate of \code{Umean} where the posterior probabilities
  for each sample are replaced by component indicator vectors
  which contain one in the element corresponding to the highest posterior
  probability while others zero (\eqn{\hat{u}_{j}^C}).  Size \eqn{n \times q}.
}
}

\references{
McLachlan GJ, Baek J, and Rathnayake SI (2011). Mixtures of factor analyzers
for the analysis of high-dimensional data.
In \emph{Mixture Estimation and Applications},
KL Mengersen, CP Robert, and DM Titterington (Eds).
Hoboken, New Jersey: Wiley, pp. 171--191.

McLachlan GJ, and Peel D (2000). \emph{Finite Mixture Models}.
New York: Wiley.
}
\author{
Geoff McLachlan, Suren Rathnayake, Jungsun Baek
}
\examples{
# Fit a MCFA model to a subset
set.seed(1)
samp_size <- dim(iris)[1]
sel_subset <- sample(1 : samp_size, 50)
model <- mcfa(iris[sel_subset, -5], g = 3, q = 2, 
                          nkmeans = 1, nrandom = 0, itmax = 100)

# plot the data points in the factor space
plot(model)

# Allocating new samples to the clusters
Y <- iris[-c(sel_subset), -5]
Y <- as.matrix(Y)
clust <- predict(model, Y)

fa_scores <- factor_scores(model, Y)
# Visualizing new data in factor space
plot_factors(fa_scores, type = "Umean", clust = clust)
}
\keyword{cluster}
\keyword{multivariate}
\keyword{models}
