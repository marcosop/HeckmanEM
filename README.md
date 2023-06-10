
# HeckmanEM: Fit Normal or Student-t Heckman Selection Models

## Overview

It performs maximum likelihood estimation for the Heckman selection
model (Normal or Student-t) using an EM-algorithm [(Lachos et al.,
2021)](https://doi.org/10.1016/j.jmva.2021.104737). It also performs
influence diagnostic through global and local influence for four
possible perturbation schema.

## Installation

You can install `HeckmanEM` from CRAN when available, or you can install
the development version from GitHub with the following commands:

``` r
# Install from CRAN (when available)
install.packages("HeckmanEM")

# Or the development version from GitHub
# install.packages("devtools")
devtools::install_github("marcosop/HeckmanEM")
```

## Available Functions

The HeckmanEM package provides tools for fitting and analyzing the
Heckman Selection model. Here is a list of the functions available in
the package:

- `HeckmanEM()`: Fits a Heckman Selection model, adjusting for potential
  selection bias.

- `HeckmanEM.criteria()`: Calculates the AIC, AICc, and BIC model
  selection criteria for a fitted Heckman Selection model.

- `HeckmanEM.envelope()`: Generates an envelope plot for a fitted
  Heckman Selection model, providing a visual model fit assessment.

- `HeckmanEM.infomat()`: Estimates standard errors for the parameters of
  a fitted model, quantifying parameter uncertainty.

- `Influence()`: Performs an influence analysis for a `HeckmanEM` object
  to assess the effect of excluding individual observations.

- `CaseDeletion()`: Carries out case deletion analysis for a `HeckmanEM`
  object, studying the impact of observation removal on parameter
  estimates.

- `rHeckman()`: Generates a random sample from a Heckman Selection model
  (Normal or Student-t), facilitating data simulation.

Each of these functions comes with complete documentation, including
examples of usage. You can access this documentation using the `help()`
function in R.

## Tutorial

``` r
# First, let's generate some data from the Heckman Selection model.

n    <- 100
nu   <- 3
cens <- 0.25

# We set a seed to make the results reproducible
set.seed(13)

# We generate the covariates matrices
w <- cbind(1, runif(n, -1, 1), rnorm(n))
x <- cbind(w[,1:2])

# We calculate the quantile of the Student's T distribution
c <- qt(cens, df = nu)

# We set the values of the parameters
sigma2   <- 1
beta     <- c(1, 0.5)
gamma    <- c(1, 0.3, -.5)
gamma[1] <- -c * sqrt(sigma2)

# We generate the data
datas <- rHeckman(x, w, beta, gamma, sigma2, rho = 0.6, nu, family = "T")
y     <- datas$y
cc    <- datas$cc

# Now, let's fit a Heckman Selection model to the data.
heckmodel <- HeckmanEM(y, x, w, cc, family = "Normal", iter.max = 50)

# We perform a case deletion analysis on the fitted model.
global <- CaseDeletion(heckmodel)

# We create a plot of the result.
plot(global)

# We perform an influence analysis with different types of perturbations.

# Case when the weight of an observation is altered
local_case <- Influence(heckmodel, type = "case-weight")
# Influential observations
local_case$influent 
# Visualization
plot(local_case)

# Case when an observation is scaled
local_scale <- Influence(heckmodel, type = "scale")
# Influential observations
local_scale$influent
# Visualization
plot(local_scale)

# Case when the response of an observation is altered
local_response <- Influence(heckmodel, type = "response")
# Influential observations
local_response$influent
# Visualization
plot(local_response)

# Exploratory case in which a covariate is altered (column 2 in this case)
local_explore <- Influence(heckmodel, type = "exploratory", colx = 2)
# Influential observations
local_explore$influent
# Visualization
plot(local_explore)
```

## Authors

This package has been developed by:

- Marcos Prates (Author, Creator, Translator), email:
  <marcosop@est.ufmg.br>
- Victor Lachos (Author)
- Dipak Dey (Author)
- Marcos Oliveira (Author, Contributor)
- Christian Galarza (Contributor)
- Katherine Loor (Contributor)

## Contributions

Contributions to the software are welcome. If you find a bug or have an
idea for an enhancement, please open an issue or submit a pull request.

## Contact

If you have any questions, feel free to open an issue or send an email
to the repository owner.

## License

This software is open source under the MIT License. Please see the
LICENSE file for more information.

<!-- ## Usage -->
<!-- `library(HeckmanEM)` will load **HeckmanEM**, a function to fit Heckman Selection models. The user can choose between two families of distribution: **Normal** or **Student-t**. -->
<!-- See `?HeckmanEM` for examples. -->
