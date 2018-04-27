[![Build Status](https://travis-ci.org/suren-rathnayake/EMMIXmfa.svg?branch=master)](https://travis-ci.org/suren-rathnayake/EMMIXmfa)

# EMMIXmfa

This package provides functions for cluster analysis using
mixtures of factor analyzers (MFA) and mixture of common factor
analyzers (MCFA) models. 

MFA and MCFA models belong to the class of finite mixture models,
and they adopt component-wise factor analyzers to multivariate data. 
Component distributions can either be from the family of multivariate normals
or from the family of multivariate _t_-distributions.
Maximum likelihood estimators of model parameters are obtained using 
the Expectation-Maximization algorithm.  

## Installation

Installingdirectly from suing the `devtools` package.

```
library(devtools)
install_github("suren-rathnayake/EMMIXmfa")
```

## Usage Example

Fitting a MFA model with there components using two factors for the Iris data available in
R can be done using,  
```
# Fit a Gaussian mixture model using MFA
mfa_fit <- mfa(Y = iris[, -5], g = 3, q = 2, sigma_type = "common", D_type = "common")

# Fit a Gaussian mixture model using MCFA
set.seed(1984)
mcfa_fit <- mcfa(Y = iris[, -5], g = 3, q = 2)
```

The groupings can be visualized in the _q_-dimensional factor space.
```
plot_factors(mcfa_fit)
```
![Plot of the factor scores](https://raw.githubusercontent.com/suren-rathnayake/misc/master/iris_mcfa_q2.png)

Adjust Rand Index
```
ari(mcfa_fit$clust, iris[, 5])
```
Functions `mfa` and `mcfa` fits multivariate normals to the data, fitting _t_-distributions can be achieved
using `mtfa` and `mctfa` function. Further, there are functions to generate data from a `emmix`
models (`rmix`), estimate factor scores (`factor_scores`), estimate adjusted Rand Index (`ari`),
find the number of misallocations (`err`), among others.
