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
This package provides functions for cluster analysis using
mixtures of factor analyzers (MFA) and mixture of common factor
analyzers (MCFA) models. 

MFA and MCFA models belong to the class of finite mixture models,
and they adopt component-wise factor analyzers to multivariate data. 
Component distributions can either be from the family of multivariate normals
or from the family of multivariate \emph{t}-distributions.
Maximum likelihood estimators of model parameters are obtained using 
the Expectation-Maximization algorithm.  

MCFA models are mixture of factor analyzers with a common component 
matrix for the  
factor loadings before the transformation of the latent factors to be 
white noise. It is designed specifically for the task of displaying the 
observed data points in a  lower (\emph{q}-dimensional) space, 
where \emph{q} is the number of factors adopted in the 
factor-analytic representation of the observed vector.
It also provides a greater reduction in the number of parameters in the model.
Component distributions can be either from the family of multivariate normals  
or from the family of multivariate \emph{t}-distributions. 
Maximum likelihood estimators of 
model parameters are obtained using the Expectation--Maximization algorithm. 
}
\details{
\tabular{ll}{
Package: \tab EMMIXmfa\cr
Type: \tab Package\cr
Version: \tab 1.3.9\cr
Date: \tab 2017-10-18\cr
License: \tab GPL\cr
}
%%~~ An overview of how to use the package, including the most important ~~
%%~~ functions ~~
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

McLachlan GJ, Baek J, and Rathnayake SI (2011). Mixtures of factor analyzers
for the analysis of high-dimensional data.
In \emph{Mixture Estimation and Applications},
KL Mengersen, CP Robert, and DM Titterington (Eds).
Hoboken, New Jersey: Wiley, pp. 171--191.
}
\keyword{ package }
\keyword{ clustering }
\keyword{ multivariate }
\keyword{ model }
\seealso{
\code{mfa}
%%~~ Optional links to other man pages, e.g. ~~
%%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
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
