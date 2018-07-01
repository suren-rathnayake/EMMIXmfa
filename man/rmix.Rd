\name{rmix}
\alias{rmix}
\title{
 Random Deviates from EMMIXmfa Models
}
\description{
Random number generator for \code{EMMIXmfa} models.
}
\usage{
rmix(n, model, ...)
}
\arguments{
  \item{model}{
  An object of class \code{'emmix'} containing a mode of
  \code{mfa, mcfa, mtfa}, or \code{mctfa}.
}
  \item{n}{
  Number of sample to generate.
}
  \item{...}{
  Not used.
  }
}
\details{
This function uses \code{rmvnorm} and
\code{rmvt} functions from the
\pkg{mvtnorm} package to generate samples
from the mixture components.

Algorithm works by first drawing a component based on
the mixture proprotion in the model, and then drawing
a sample from the component distribution.
}
\value{
  A numeric matrix with samples drawn in rows.
}
% \references{
% McLachlan GJ, and Peel D (2000).
% \emph{Finite Mixture Models}. New York: Wiley.
% }
% \author{
% Geoffrey McLachlan, Suren Rathnayake
% }
%%\note{
%%}
\examples{
\donttest{
 set.seed(1)
 model <- mcfa(iris[, -5], g=3, q=2, nkmeans=1, nrandom=1)
 dat <- rmix(n = 10, model = model)
}
}
\keyword{multivariate}
\keyword{distribution}
\keyword{datagen}
