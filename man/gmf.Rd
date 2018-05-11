\name{gmf}
\alias{gmf}
\title{
General Matrix Factorization
}
\description{
Performs a matrix factorization on the given data set.
The factorization is done using a stochastic gradient decent method.
}
\usage{
gmf(Y, q, maxit = 1000, lambda = 0.01, cor_rate = 0.9)
}
\arguments{
  \item{Y}{
   data matrix containing all numerical values.
  }
  \item{maxit}{
  maximum number of iterations.
  }
  \item{q}{
  number of factors.
  }
  \item{lambda}{
  initial learning rate.
  }
    \item{cor_rate}{
  correction rate.
  }
}
\details{
Unsupervised matrix factorization of a \eqn{n \times p} data matrix
\eqn{Y} can be expressed as,
\deqn{
 Y^{\top} \approx A B^{\top},
}
where \eqn{A} is a \eqn{p \times q} matrix and \eqn{B} is
\eqn{n \times q} matrix.
With this matrix factorization method, one replaces 
the \eqn{i}th row in matrix \eqn{Y} by the \eqn{i}th row in matrix \eqn{B}.
The matrices \eqn{A} and \eqn{B} are chosen to minimize an objective
function \eqn{f(Y, A, B)} with under constraints specific 
to the matrix factorization method.

It is imperative that columns of the data matrix be on the same scale.
Otherwise, it may not be possible to obtain a factorization
of the data using this approach.
}
\value{
A list containing,
\item{A }{A numeric matrix of size \eqn{p \times q}}
\item{B }{A numeric matrix of size \eqn{n \times q} matrix}
%% ...
}
\references{
Nikulin V, Huang T-H, Ng SK, Rathnayake SI, & McLachlan GJ (2011).
A very fast algorithm for matrix factorization.
\emph{Statistics & Probability Letters} \strong{81}, 773--782.
}

\examples{
lst <- gmf(iris[, -5], q = 2, maxit = 100)
}
\keyword{multivariate}
\keyword{methods}
