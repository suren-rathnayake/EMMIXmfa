\name{factor_scores}
\alias{factor_scores}
\alias{factor_scores.mcfa}
\alias{factor_scores.mctfa}
\title{
 Computes Factor Scores.
}
\description{
This function computes factor scores given a data set and
a \code{EMMIXmfa} model. Using factor scores, 
we can represent the original data point \eqn{y_j} in 
\emph{q}-dimensional reduced space by plotting the estimated
conditional expectation of the unobservable factors 
\eqn{U_{ij}} given \eqn{y_j} and the component membership
which can be expressed by,
\deqn{
\hat{u}_{ij} = E_{\hat{\Psi}}\{U_{ij} \mid y_j, z_{ij} = 1\}.
}
have obtained the estimated factor score
corresponding to \eqn{y_j} (\eqn{j = 1, \dots, n}) as 
\eqn{\hat{u}_{j}} as the mean of the estimated conditional expected
component scores \eqn{\hat{u}_{ij}} over the estimated posterior
distribution as of the component membership of \eqn{y_j},
\deqn{
\hat{u}_{j} = \sum_{i=1}^g \tau_i(y_j; \hat{\Psi}) \hat{u}_{ij}.
}
For MFA, we have 
\deqn{
\hat{u}_{ij} = \hat{\beta}_i^T (y_j - \hat{\mu}_i),
}
and
\deqn{
\hat{u}_{ij} = \sum_{i=1}^g \tau_i(y_j;  \hat{\Psi}) \hat{\beta}_i^T 
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

With M\eqn{t}FA and MC\eqn{t}FA, the distribution of 
the expression of \eqn{\hat{u}_{ij}} and \eqn{\hat{u}_{j}}
have the same form as those of MFA and MCFA, respectively.

An alternative estimate of \eqn{u_j}, the posterior expectation
of the factor corresponding to the \emph{j}th observation \eqn{y_j}, is
defined by replacing \eqn{\tau_i(y_j;\,\hat{\Psi})} by \eqn{\hat{z}_{ij}}, 
where
\eqn{\hat{z}_{ij} = 1}, if \eqn{\hat{\tau}_i(y_j; \hat{\Psi})}
>= \eqn{\hat{\tau_h}(y_j; \hat{\Psi})  
(h=1,\,\dots,\,g; h \neq i)}, else 
\eqn{\hat{z}_{ij} = 0}.


}
\usage{
factor_scores(model, Y, ...)
\method{factor_scores}{mcfa}(model, Y, tau = NULL, clust= NULL, ...)
\method{factor_scores}{mctfa}(model, Y, tau = NULL, clust= NULL, ...)
}
\arguments{
  \item{model}{
  Model of class \code{"mcfa"} or \code{"mctfa"}.
}
  \item{Y}{
  Data matrix with variables in columns in the same order as used in
  model estimation.
}
  \item{tau}{
  Optional. Posterior probabilities of belonging to the components
  in the mixture model. If not provided, they will be estimated. 
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
    unobservable factors given the data and the component membership.
    Size is \eqn{n \times q \times g}, where \code{n} is the number of sample,
    \code{q} is the number of factors and \code{g} is the number components.
}
  \item{Umean}{
   Means of the estimated conditional expected factors scores over
   estimated posterior distributions. Size \eqn{n \times q}.
}
 \item{Uassign}{
 Alternative estimate of \code{Umean} where the posterior probabilities
 for each sample are replaced by component indicator vectors 
 which contain one in the element corresponding to the highest posterior 
 probability while others zero.  Size \eqn{n \times q}.
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
%%\note{
%%}
\examples{

# Fit a MCFA model to a subset 
set.seed(1)
samp_size <- dim(iris)[1]
sel_subset <- sample(1 : samp_size, 75)
model <- mcfa(iris[sel_subset, -5], g=3, q=2, nkmeans=1, nrandom=0)

# plot the data points in the factor space
#plot(model)

# Allocating new samples to the clusters
Y <- iris[-c(sel_subset), -5]
Y <- as.matrix(Y)
clust <- predict(model, Y)

fa_scores <- factor_scores(model, Y)
# Visualizing new data in factor space
plot_factors(fa_scores, type="Umean", clust=clust)
}
\keyword{ cluster}
\keyword{ multivariate}
\keyword{ models}
