\name{ari}
\alias{ari}
\title{
Computes adjusted Rand Index
}
\description{
Computes adjusted Rand index. 
}
\usage{
ari(cls, hat_cls)
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
Measures the agreement between two sets of partitions. 
The upper bound of 1 implies perfect agreement. 
The expected value is zero if the partitions are random.
}
\value{
Scaler specifying how closely two partitions agree.    
}
\references{
Hubert L, and Arabie P (1985). Comparing Partitions.
\emph{Journal of the Classification} \strong{2}, 193--218.
}

\seealso{
\code{\link{minmis}}
}
\examples{
set.seed(1984)
Y <- scale(iris[, -5])
model <- mfa(Y, g = 3, q = 3, nkmeans = 1, nrandom = 0)
#
ari(model$clust, iris[, 5])
#
minmis(model$clust, iris[, 5])
}
\keyword{ cluster }
