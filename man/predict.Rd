\name{predict.emmix}
\alias{predict.emmix}
\alias{predict.mfa}
\alias{predict.mcfa}
\alias{predict.mtfa}
\alias{predict.mctfa}
\title{
Extend Clustering to New Observations 
}
\description{
Given a fitted model of class 
\code{'emmix'} (or of class
\code{'mfa'}, \code{'mcfa'}, \code{'mtfa'} and 
\code{'mctfa'}) and 
a set of observations, 
the \code{predict} function provides
a clustering of new observations based on
fitted values.
}
\usage{
\method{predict}{emmix}(object, Y, ...)
}
\arguments{
  \item{object}{
    An object of class \code{'emmix'}.
}
  \item{Y}{
    A data matrix with variable in the same
    column locations as the data used in 
    fitting the model \code{object}.
}
  \item{\dots}{
     Not used.
}
}
\details{
A vector integers of length equal to number of
rows in the data. The integers range from 1 to
\eqn{g} where \eqn{g} in the number of components
in the model.
}
\examples{
\dontrun{
set.seed(42)
test <- sample(1 : nrow(iris), 100)
model <- mfa(iris[test, -5], g=3, q=3, itmax=500, nkmeans=3, nrandom=5)
pred_clust <- predict(model, iris[-test, -5])
err(pred_clust, iris[-test, 5])
}
}
\keyword{ cluster }
