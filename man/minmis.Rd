\name{minmis}
\alias{minmis}
\title{
 Minimum Number of Misallocations
}
\description{
Given two vectors each corresponding to a set of categories,
this function finds the minimum number of misallocations
by rotating the categories.
}
\usage{
minmis(cls, hat_cls)
}
\arguments{
  \item{cls}{
    A numeric or character vector of labels.
}
  \item{hat_cls}{
    A numeric or character vector of labels same length as \code{cls}.
}
}
\details{
Rotates the categories for all possible permutations, and
returns the minimum number of misallocations.
The number of categories in each set of labels does not need to
be the same.
It may take several minutes to compute when the number of categories
is large.
}
\value{
  Integer specifying the minimum number of misallocations.
}
\seealso{
\code{\link{ari}}
}
\examples{
set.seed(1984)
Y <- scale(iris[, -5])
model <- mcfa(Y, g = 3, q = 3, nkmeans = 1, nrandom = 0, itmax = 200)
ari(model$clust, iris[, 5])
minmis(model$clust, iris[, 5])
}
\keyword{ cluster }
