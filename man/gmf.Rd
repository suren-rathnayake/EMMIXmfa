\name{gmf}
\alias{gmf}
%- Also NEED an '\alias' for EACH other topic documented here.
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
%- maybe also 'usage' for other objects documented here.
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
Unsupervised matrix factorization of the data matrix 
\eqn{X} can be expressed as,
\deqn{
 X \approx A B,
}
where \eqn{A} is a \eqn{p \times q} matrix and \eqn{B} is 
\eqn{q \times n} matrix.
With matrix factorization one expects to replace the 
sample \eqn{x_{i}}, the \eqn{i}th row
in matrix \eqn{X}, by the \eqn{i}th column, \eqn{b_{i}}, in matrix \eqn{B}.
The matrices \eqn{A} and \eqn{B} are chosen to minimize objective
function \eqn{f(X, A, B)} with under specific constraints
to the matrix factorization method. 

It is imperative that columns of the data matrix be on the same scale.
Otherwise, it may not be possible to obtain a factorization
of the data using this approach.
}
\value{
A list containing,
\item{A }{\eqn{p \times q} matrix}
\item{B }{\eqn{q \times n} matrix}
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
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~methods }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~algebra }% __ONLY ONE__ keyword per line
