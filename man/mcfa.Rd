\name{mcfa}
\alias{mcfa}
\alias{mcfa.default}
\alias{mctfa}
\alias{mctfa.default}
\title{
Mixture of Common Factor Analyzers
}
\description{
Functions for fitting mixtures of common factor analyzers (MCFA)
models.
MCFA models are mixture of factor analyzers
(belong to the class of multivariate finite mixture models)
with a common component matrix for the
factor loadings before the transformation of the latent factors to be
white noise. It is designed specifically for the task of displaying the
observed data points in a  lower (\emph{q}-dimensional) space,
where \emph{q} is the number of factors adopted in the
factor-analytic representation of the observed vector.

The \code{mcfa} function fits mixtures common factor analyzers
where the  components distributions belong to the family of
multivariate normal distributions.
The \code{mctfa} function fits
mixtures of common \emph{t}-factor analyzers where
the component distributions corresponds to multivariate
\emph{t} distributions.
Maximum likelihood estimates of the model parameters are obtained
using the Expectation--Maximization algorithm.

}
\usage{
mcfa(Y, g, q, itmax = 500, nkmeans = 5, nrandom = 20,
  tol = 1.e-5, init_clust = NULL, init_para = NULL,
  init_method = NULL, conv_measure = 'diff',
  warn_messages = TRUE, ...)
mctfa(Y, g, q, itmax = 500, nkmeans = 5, nrandom = 20,
  tol = 1.e-5, df_init = rep(30, g), df_update = TRUE,
  init_clust = NULL, init_para = NULL, init_method = NULL,
  conv_measure = 'diff', warn_messages = TRUE, ...)
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
  \item{init_clust}{
  A vector or matrix consisting of partition of samples to be used
  in the EM algorithm. For matrix of partitions, columns must corresponds
  individual partitions of the data. Optional.
}
  \item{init_para}{
  A list containing model parameters to be used as initial
  parameter estimates for the EM algorithm. Optional.
}
  \item{init_method}{
  To determine how the initial parameter values are computed. See Details.
}
  \item{conv_measure}{
  The default \code{'diff'} stops the EM iterations if
  |\eqn{l^{(k+1)}} - \eqn{l^{(k)}}| < \code{tol} where
  \eqn{l^{(k)}} is the log-likelihood at the \eqn{k}th EM iteration.
  If \code{'ratio'}, then the convergence of the EM steps is measured
  using the |(\eqn{l^{(k+1)}} - \eqn{l^{(k)}})/\eqn{l^{(k+1)}}|.
}
  \item{df_init}{
  Initial values of the degree of freedom parameters for \code{mctfa}.
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
With \code{init_method = NULL}, the default,
model parameters are initialized using all available methods.
With the \code{init_method = "rand-A"}, the initialization of
the parameters is done using the procedure in
Baek et al. (2010) where initial values for elements of
\eqn{A} are drawn from the \eqn{N(0, 1)} distribution.
This method is appropriate when the columns of the data
are on the same scale. The
\code{init_method = "eigen-A"}
takes the first \eqn{q} eigenvectors of \eqn{Y} as the
initial value for the loading matrix \eqn{A}. 
If \code{init_method = "gmf"} then the data are factorized using
\code{gmf} with \eqn{q} factors and the resulting loading
matrix is used as the initial value for \eqn{A}.

If specified, the optional argument \code{init_para}
must be a list or an object of class \code{mcfa} or \code{mctfa}.
When fitting an \code{mcfa} model, only the
model parameters \code{q}, \code{g},
\code{pivec}, \code{A}, \code{xi},
\code{omega}, and \code{D} are extracted from
\code{init_para}, while one extra parameter
\code{nu} is extracted when fitting \code{mctfa}.
Everything else in \code{init_para} will be discarded.
}
\value{
Object of class \code{c("emmix", "mcfa")} or \code{c("emmix",
"mctfa")} containing the fitted model parameters is returned.
Details of the components are as follows:
 \item{g}{
 Number of mixture components.
}
 \item{q}{
  Number of factors.
}
 \item{pivec}{
  Mixing proportions of the components.
}
  \item{A}{
  Loading matrix. Size \eqn{p \times q}.
}
  \item{xi}{
  Matrix containing factor means for components in columns.
  Size \eqn{q \times g}.
}
  \item{omega}{
  Array containing factor covariance matrices for components.
  Size \eqn{q \times q \times g}.
}
  \item{D}{
  Error covariance matrix. Size \eqn{p \times p.}
}
  \item{Uscores}{
    Estimated conditional expected component scores of the
    unobservable factors given the data and the component membership.
    Size \eqn{n \times q \times g}.
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
  \item{clust}{Cluster labels.
}
  \item{tau}{Posterior probabilities.
}
  \item{logL}{Log-likelihood at the convergence.
}
  \item{BIC}{Bayesian information criterion.
}
  \item{warn_msg}{Description of error messages, if any.
}
}
\references{
Baek J, McLachlan GJ, and Flack LK (2010). Mixtures of factor analyzers
with common factor loadings: applications to the clustering and visualisation
of high-dimensional data. \emph{IEEE Transactions on Pattern Analysis and
Machine Intelligence} \strong{32}, 2089--2097.

Baek J, and McLachlan GJ (2011). Mixtures of common \emph{t}-factor analyzers
for clustering highdimensional microarray data.
\emph{Bioinformatics} \strong{27}, 1269--1276.

McLachlan GJ, Baek J, and Rathnayake SI (2011). Mixtures of factor analyzers
for the analysis of high-dimensional data.
In \emph{Mixture Estimation and Applications},
KL Mengersen, CP Robert, and DM Titterington (Eds).
Hoboken, New Jersey: Wiley, pp. 171--191.
}
\author{
Suren Rathnayake, Jangsun Baek, Geoff McLachlan
}
\seealso{
\code{\link{mfa}}, \code{\link{plot_factors}}
}
\examples{
mcfa_fit <- mcfa(iris[, -5], g = 3, q = 3, itmax = 25,
                  nkmeans = 5, nrandom = 5, tol = 1.e-5)

plot(mcfa_fit)
\donttest{
mctfa_fit <- mcfa(iris[, -5], g = 3, q = 3, itmax = 500,
                  nkmeans = 5, nrandom = 5, tol = 1.e-5, df_update = TRUE)
}
}
\keyword{cluster}
\keyword{multivariate}
\keyword{models}
