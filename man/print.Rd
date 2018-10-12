\name{print.emmix}
\alias{print.emmix}
\alias{summary.emmix}
\alias{print.mfa}
\alias{print.mcfa}
\alias{print.mtfa}
\alias{print.mctfa}
\title{
Print Method for Class 'emmix'
}
\description{
Prints a formatted model parameters of
\code{EMMIXmfa} objects.
}
\usage{
\method{print}{emmix}(x, ...)
\method{summary}{emmix}(object, ...)
}
\arguments{
  \item{x, object}{
    An object of class \code{'emmix'}.
}
  \item{\dots}{
     Not used.
}
}
\details{
Prints the formatted model parameter values to
the screen.
}
\examples{
\donttest{
set.seed(1984)
Y <- scale(iris[, -5])
model <- mcfa(Y, g = 3, q = 3, nkmeans = 1, nrandom = 0, itmax = 100)
#
print(model)
summary(model)
}
}
\keyword{print}
