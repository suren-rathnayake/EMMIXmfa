\name{mfa}
\alias{mfa}
\alias{mfa.default}
\alias{plot.mfa}
\alias{mtfa}
\alias{mtfa.default}
\alias{plot.mtfa}
\title{
Mixtures of Factor Analyzers
}
\description{
Functions for fitting mixtures of factor analyzers (MFA) and
mixtures of \emph{t}-factor analyzers (M\emph{t}FA) to data.
Maximum Likelihood estimates of the model parameters are obtained
using the Alternating Expectation Conditional Maximization (AECM)
algorithm.

In the case of MFA, component distributions belong to the family of
multivariate normal distributions, while with M\eqn{t}FA
the component distributions correspond to multivariate
\emph{t} distributions.
}
\usage{
mfa(Y, g, q, itmax = 500, nkmeans = 20, nrandom = 20,
  tol = 1.e-5, sigma_type = 'common', D_type = 'common', init_clust = NULL,
  init_para = NULL, conv_measure = 'diff', warn_messages = TRUE, ...)
mtfa(Y, g, q, itmax = 500, nkmeans = 20, nrandom = 20,
  tol = 1.e-5, df_init = rep(30, g), df_update = TRUE,
  sigma_type = 'common', D_type = 'common', init_clust = NULL,
  init_para = NULL, conv_measure = 'diff', warn_messages = TRUE, ...)
}
\arguments{
 \item{Y}{
  A matrix or a data frame of which rows correspond to
  observations and columns to variables.
}
 \item{g}{
  Number of components.
}
 \item{q}{
  Number of factors.
}
 \item{itmax}{
  Maximum number of EM iterations.
}
 \item{nkmeans}{
  The number of times the k-means algorithm to be used in partition
  the data into \code{g} groups. These groupings are then used in
  initializing the parameters for the EM algorithm.
}
 \item{nrandom}{
  The number of random \code{g}-group partitions for the data to be used
  initializing the EM algorithm.
}
 \item{tol}{
  The EM algorithm terminates if the measure of convergence falls below
  this value.
}
 \item{sigma_type}{
  To specify whether the covariance matrices (for \code{mfa})
  or the scale matrices (for \code{mtfa}) of the components
  are constrained 
  to be the same (default, \code{sigma_type = "common"})
  or not (\code{sigma_type = "unique"}).
}
 \item{D_type}{
  To specify whether the diagonal error covariance matrix is common to all
  the components or not. If \code{sigma_type = "unique"}, then
  \code{D_type} could either be \code{"common"}
  (the default) to each component, or \code{"unique"}.
  If the \code{sigma_type = "common"}, then
  \code{D_type} must also be \code{ "common"}.
}
 \item{init_clust}{
  A vector or matrix consisting of partition of samples to be used
  in the EM algorithm. For matrix of partitions, columns must corresponds
  individual partitions of the data. Optional.
}
  \item{init_para}{
  A list containing model parameters to be used as initial
  parameter estimates for the EM algorithm. Optional.
}
 \item{conv_measure}{
  The default \code{'diff'} stops the EM iterations if
  |\eqn{l^{(k+1)}} - \eqn{l^{(k)}}| < \code{tol} where
  \eqn{l^{(k)}} is the log-likelihood at the \eqn{k}th EM iteration.
  If \code{'ratio'}, then the convergence of the EM steps is measured
  using the |(\eqn{l^{(k+1)}} - \eqn{l^{(k)}})/\eqn{l^{(k+1)}}|.
}
 \item{df_init}{
  Initial values of the degree of freedom parameters for \code{mtfa}.
}
 \item{df_update}{
  If \code{df_update = TRUE} (default), then the degree of freedom parameters
  values will be updated during the EM iterations.
  Otherwise, if \code{df_update = FALSE}, they will be fixed at the initial
  values specified in \code{df_init}.
}
 \item{warn_messages}{
  With \code{warn_messages = TRUE} (default), the output would
  include some description of the reasons where, if any, the model fitting
  function failed to provide a fit for a given set of initial
  parameter values.
}
\item{\dots}{
  Not used.
}
}
\details{
Cluster a given data set using mixtures of factor analyzers or
approach or using mixtures of \emph{t}-factor analyzers.
}
\value{
Object of class \code{c("emmix", "mfa")} or \code{c("emmix",
"mtfa")} containing the fitted model parameters is returned.
Details of the components are as fellows:
 \item{g}{
  Number of mixture components.
}
 \item{q}{
  Number of factors.
}
 \item{pivec}{
 Mixing proportions of the components.
}
 \item{mu}{
  Matrix containing estimates of component means (in columns) 
  of mixture component. Size \eqn{p \times g}.
}
 \item{B}{
  Array containing component dependent loading matrices. Size
  \eqn{p \times q \times g}.
}
  \item{D}{
  Estimates of error covariance matrices. If \code{D_type = "common"}
  was used then \code{D} is \eqn{p \times p} matrix common to
  all components,  if \code{D_type = "unique"}, then \code{D} is a
  size \eqn{p \times p \times g} array.
}
  \item{v}{
  Degrees of freedom for each component.
}
  \item{logL}{
  Log-likelihood at the convergence.
}
  \item{BIC}{
  Bayesian information criterion.
}
  \item{tau}{
  Matrix of posterior probabilities for the data
  used based on the fitted values. Matrix of size \code{n by g}.
}
  \item{clust}{
  Vector of integers 1 to g indicating cluster allocations
  of the observations.
}
  \item{Uscores}{
    Estimated conditional expected component scores of the
    unobservable factors given the data and the component membership.
    Size is Size \eqn{n \times q \times g}.
}
  \item{Umean}{
   Means of the estimated conditional expected factors scores over
   estimated posterior distributions. Size \eqn{n \times q}.
}
 \item{Uclust}{
 Alternative estimate of \code{Umean} where the posterior probabilities
 for each sample are replaced by component indicator vectors
 which contain one in the element corresponding to the highest posterior
 probability while others zero.  Size \eqn{n \times q}.
}
  \item{ERRMSG}{
  Description of messages, if any.
}
  \item{D_type}{
  Whether common or unique error covariance is used, as specified in
  model fitting.
}
  \item{df_update}{
  Whether the degree of freedom parameter 
  (\code{v}) was fixed or estimated (only for \code{mtfa}).
}
}
\references{
Ghahramani Z, and Hinton GE (1997). 
The EM algorithm for mixture of factor analyzers.
\emph{Technical Report, CRG-TR-96-1}, University of Toronto, Toronto.

McLachlan GJ, Bean RW, Ben-Tovim Jones L (2007). 
Extension of the mixture of factor analyzers model to incorporate the 
multivariate \emph{t} distribution.
\emph{Computational Statistics & Data Analysis}, \strong{51}, 5327--5338.

McLachlan GJ, Baek J, and Rathnayake SI (2011). Mixtures of factor analyzers
for the analysis of high-dimensional data.
In \emph{Mixture Estimation and Applications},
KL Mengersen, CP Robert, and DM Titterington (Eds).
Hoboken, New Jersey: Wiley, pp. 171--191.

McLachlan GJ, Peel D, and Bean RW (2003). 
Modelling high-dimensional data by mixtures of factor analyzers.
\emph{Computational Statistics & Data Analysis} \strong{41}, 379--388.
}
\author{
Suren Rathnayake, Geoffrey McLachlan
}
\seealso{
\code{\link{mcfa}}
}
\examples{
model <- mfa(iris[, -5], g=3, q=2, itmax=200, nkmeans=1, nrandom=5)
summary(model)
\donttest{
model <- mtfa(iris[, -5], g=3, q=2, itmax=200, nkmeans=1, nrandom=5)
}  
}
\keyword{cluster}
\keyword{multivariate}
\keyword{models}
