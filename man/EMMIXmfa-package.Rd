\name{EMMIXmfa-package}
\alias{EMMIXmfa-package}
\alias{emmixmfa-package}
\alias{EMMIXmfa}
\alias{emmixmfa}
\docType{package}
\title{
Mixture Models with Component-Wise Factor Analyzers
}
\description{
This package provides functions for fitting
mixtures of factor analyzers (MFA) and
mixtures of common factor analyzers (MCFA) models.

MFA and MCFA models belong to the class of finite mixture models,
that adopt factor models for the component-covariance matrices.
More specifically, under the factor model,
the correlation between feature variables can be estimated by their
linear dependance of these variables on a smaller small number
\emph{q} of (unobservable) latent factors.
The component distributions can be either from the family of
multivariate normals or from the family of multivariate
\emph{t}-distributions.
Maximum likelihood estimation of model parameters
is implemented using the Expectation-Maximization algorithm.

The joint distribution of the factors and errors can be taken
to be either the multivariate normal or \emph{t}-distribution.
The factor analytic representation of the component-covariance
matrices as a way of dimension reduction in that it
enable the mixture distributions
to be fitted  to data with dimension \emph{p}
relatively large compared to the sample size \emph{n}.

Unlike MFA, MCFA models can be used to display
observed data points in the \emph{q}-dimensional factor space.
MCFA also provides a greater reduction in the number of
parameters in the model.
}
\details{
\tabular{ll}{
Package: \tab EMMIXmfa\cr
Type: \tab Package\cr
Version: \tab 2.0.0\cr
Date: \tab 2018-05-01\cr
License: \tab GPL (>= 2)\cr
}
}
\author{
%Who wrote it
Suren Rathnayake, Jangsun Baek, Geoffrey McLachlan
}
\references{
McLachlan GJ, Peel D, and Bean RW (2003). Modelling high-dimensional data
by mixtures of factor analyzers.
\emph{Computational Statistics & Data Analysis} \strong{41}, 379--388.

McLachlan GJ, Bean RW, Ben-Tovim Jones L (2007). Extension of the mixture
of factor analyzers model to incorporate the multivariate \emph{t}
distribution.
\emph{Computational Statistics & Data Analysis}, \strong{51}, 5327--5338.

Baek, J. and McLachlan, G.J. (2008). Mixtures of factor analyzers with
common factor loadings for the clustering and visualisation of
high-dimensional data. \emph{Technical Report NI08018-SCH}, Preprint Series
of the Isaac Newton Institute for Mathematical Sciences, Cambridge.

Baek J, McLachlan GJ, and Flack LK (2010). Mixtures of factor analyzers
with common factor loadings: applications to the clustering and visualisation
of high-dimensional data. \emph{IEEE Transactions on Pattern Analysis and
Machine Intelligence} \strong{32}, 2089--2097.

Baek J, and McLachlan GJ (2011). Mixtures of common \emph{t}-factor analyzers
for clustering highdimensional microarray data.
\emph{Bioinformatics} \strong{27}, 1269--1276.
}
\keyword{package}
\keyword{clustering}
\keyword{multivariate}
\keyword{model}
\seealso{
\code{\link{EMMIXcskew}}, \code{\link{EMMIXuskew}}
}
\examples{
\dontrun{
set.seed(1)
Y <- scale(iris[, -5])
mfa_model <- mfa(Y, g = 3, q = 3)
mtfa_model <- mtfa(Y, g = 3, q = 3)
mcfa_model <- mcfa(Y, g = 3, q = 3)
mctfa_model <- mctfa(Y, g = 3, q = 3)
}
}
