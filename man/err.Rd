\name{err}
\alias{err}
\title{
 Minimum Number of Mis-Allocations.
}
\description{
Given two vectors each corresponding to a set of categories, 
this function finds the minimum number of mis-allocations
by rotating the categories.
}
\usage{
err(cls, hat_cls)
}
\arguments{
  \item{cls}{
    Vector of labels.
}
  \item{hat_cls}{
    Vector of labels same length as \code{cls}.
}
}
\details{
Rotates the categories for all possible permutations, and
returns the minimum number of mis-allocations.
The number of categories in each set of labels does not need to
be the same.
It make take several minutes to compute when the number of categories 
is large.
}
\value{
  Integer specifying the minimum number of mis-allocations. 
}
\author{
  Suren Rathnayake
}
\seealso{
\code{\link{ari}}
} 
\examples{
\dontrun{
set.seed(1984)
Y <- scale(iris[, -5])
model <- mcfa(Y, g = 3, q = 3, nkmeans = 1, nrandom = 0)
#
ari(model$clust, iris[, 5])
#
err(model$clust, iris[, 5])
}
}
\keyword{ cluster }
