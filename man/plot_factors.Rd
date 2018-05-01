\name{plot_factors}
\alias{plot_factors}
\title{
 Plot Function for Factor Scores.
}
\description{
 Plot function for factor scores given factor
 score matrix of fitted model.
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
The \code{type} should either be \code{"Uscores"}, \code{"Uclust"} or
the default \code{"Umean"}.
If \code{type = "Uscores"}, then the estimated conditional expected component
scores of the unobservable factors given the data and the component membership
are plotted.
If \code{type = "Umean"}, then the means of the estimated conditional
expected factors scores over estimated posterior distributions
are plotted.
If \code{type = "Uclust"}, then an alternative estimate of \code{"Umean"},
where the posterior probabilities are replaced by component
indicator vector, is plotted.
}

\references{
McLachlan GJ, and Peel D (2000).
\emph{Finite Mixture Models}. New York: Wiley.

McLachlan GJ, Baek J, and Rathnayake SI (2011). Mixtures of factor analyzers
for the analysis of high-dimensional data.
In \emph{Mixture Estimation and Applications},
KL Mengersen, CP Robert, and DM Titterington (Eds).
Hoboken, New Jersey: Wiley, pp. 171--191.
}

\author{
Geoffrey McLachlan, Suren Rathnayake, Jungsun Baek
}

\examples{

# Visualizing data used in model estimation
set.seed(1)
inds <- dim(iris)[1]
indSample <- sample(1 : inds, 75)
model <- mcfa (iris[indSample, -5], g = 3, q = 2, nkmeans = 1, nrandom = 0)
err (model$clust, iris[indSample, 5])

#same as plot_factors(model, tyep = "Umean", clust = model$clust)
plot (model)

#can provide alternative groupings of samples via plot_factors
plot_factors (model, clust = iris[indSample, 5])

#same as plot_factors(model, tyep = "Uclust")
plot (model, type = "Uclust")

Y <- iris[-c(indSample), -5]
Y <- as.matrix(Y)
clust <- predict(model, Y)
err(clust, iris[-c(indSample), 5])

fac_scores <- factor_scores(model, Y)
plot_factors (fac_scores, type="Umean", clust = clust)
plot_factors (fac_scores, type="Umean", clust = iris[-c(indSample), 5])

}
\keyword{ cluster}
\keyword{ multivariate}
\keyword{ models}
