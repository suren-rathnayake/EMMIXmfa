\name{plot_factors}
\alias{plot_factors}
\title{
 Plot Function for Factor Scores
}
\description{
 Plot functions for factor scores.
}
\usage{
plot_factors(scores, type = "Umean",
    clust=if (exists('clust', where = scores)) scores$clust else NULL,
    limx = NULL, limy = NULL)
}
\arguments{
  \item{scores}{
  A list containing factor scores specified by
  \code{Umean}, \code{Uclust} or \code{Uscores}, or a
  model of class \code{mcfa, mctfa, mfa}, or \code{mtfa}.
}
  \item{type}{
  What type of factor scores are to be plotted. See Details.
}
  \item{clust}{
  Indicators of belonging to components. If available, they will be
  portrayed in plots.
  If not provided, looks for \code{clust} in \code{scores},
  and sets to \code{NULL} if still not available.
}
  \item{limx}{
  Numeric vector. Values in \code{limx} will only be used in setting
  the x-axis range for 1-D and 2-D plots.
}
  \item{limy}{
  Numeric vector. Values in \code{limy} will only be used in setting
  the y-axis range for 1-D and 2-D plots.
}
}
\details{
When the factor scores were obtained using \code{mcfa}
or \code{mctfa}, then a visualization of the group structure 
can be obtained by plotting the factor scores.
In the case of \code{mfa} and \code{mtfa}, the factor scores
simply corresponds to white noise.

The \code{type} should either be \code{"Uscores"}, \code{"Uclust"} or
the default \code{"Umean"}. See \code{factor_scores} for a detailed
description of the factor scores.
}
\references{
McLachlan GJ, Baek J, and Rathnayake SI (2011). 
Mixtures of factor analyzers
for the analysis of high-dimensional data.
In \emph{Mixture Estimation and Applications},
KL Mengersen, CP Robert, and DM Titterington (Eds).
Hoboken, New Jersey: Wiley, pp. 171--191.

McLachlan GJ, and Peel D (2000).
\emph{Finite Mixture Models}. New York: Wiley.
}
\author{
Geoffrey McLachlan, Suren Rathnayake, Jungsun Baek
}
\examples{
# Visualizing data used in model estimation
set.seed(1)
inds <- dim(iris)[1]
indSample <- sample(1 : inds, 50)
model <- mcfa (iris[indSample, -5], g = 3, q = 2, 
                nkmeans = 1, nrandom = 0, itmax = 150)
minmis(model$clust, iris[indSample, 5])

#same as plot_factors(model, tyep = "Umean", clust = model$clust)
plot(model)

#can provide alternative groupings of samples via plot_factors
plot_factors(model, clust = iris[indSample, 5])

#same as plot_factors(model, tyep = "Uclust")
plot(model, type = "Uclust")

Y <- iris[-c(indSample), -5]
Y <- as.matrix(Y)
clust <- predict(model, Y)
minmis(clust, iris[-c(indSample), 5])

fac_scores <- factor_scores(model, Y)
plot_factors(fac_scores, type = "Umean", clust = clust)
plot_factors(fac_scores, type = "Umean", clust = iris[-c(indSample), 5])
}
\keyword{cluster}
\keyword{multivariate}
\keyword{models}
