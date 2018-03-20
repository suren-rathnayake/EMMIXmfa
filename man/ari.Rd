\name{ari}
\alias{ari}
\title{
Computes adjusted Rand Index
}
\description{
Computes adjusted Rand Index. 
}
\usage{
ari(cls, hat_cls)
}
\arguments{
  \item{cls}{
    Vector containing labels or classes.
}
  \item{hat_cls}{
    Vector of labels same length as \code{cls}.
}
}
\details{
Measures the agreement between two set of partitions. 
The upper bound 1 implies perfect agreement. 
Expected value is zero if the partitions are random.
}
\value{
Scaler specifying how closely two partitions agree.    
}
\references{
Hubert L, and Arabie P (1985). Comparing Partitions.
\emph{Journal of the Classification} \strong{2}, 193--218.
}

\seealso{
\code{\link{err}}
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
