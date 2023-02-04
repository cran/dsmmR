README
================

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- Use `build_readme()` to edit the .md file from the .rmd file -->
<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/dsmmR)](https://CRAN.R-project.org/package=dsmmR)
[![R-CMD-check](https://github.com/Mavrogiannis-Ioannis/dsmmR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Mavrogiannis-Ioannis/dsmmR/actions/workflows/R-CMD-check.yaml)
[![CRAN
logs](https://cranlogs.r-pkg.org/badges/dsmmR)](https://CRAN.R-project.org/package=dsmmR)

<!-- badges: end -->

# dsmmR

The **dsmmR** R package allows the user to estimate, simulate and define
different Drifting semi-Markov model (DSMM) specifications.

## Installation

``` r
# Install the released version from CRAN
install.packages('dsmmR')
# Or the development version from GitHub
# install.packages("devtools")
devtools::install_github("Mavrogiannis-Ioannis/dsmmR")
```

## High-level documentation

The main functions of **dsmmR** are the following:

- `fit_dsmm()` : estimate a DSMM (parametric or non-parametric
  estimation is possible).
- `parametric_dsmm()` : define a parametric DSMM.
- `nonparametric_dsmm()` : define a non-parametric DSMM.
- `simulate()` : simulate a sequence from a DSMM.
- `get_kernel()` : obtain the Drifting semi-Markov kernel.

### Theory overview

Drifting semi-Markov models are best suited to capture non-homogeneities
which evolve in a linear (or polynomial) way. For example, through this
approach we account for non-homogeneities that occur from the intrinsic
evolution of the system or from the interactions between the system and
the environment.

For a detailed introduction in Drifting semi-Markov models consider the
documentation through `?dsmmR`.

For an extensive description of this approach, consider visiting the
complete documentation of the package on the [official CRAN
page](https://CRAN.R-project.org/package=dsmmR).

### Estimation

The easiest way to use **dsmmR** is through the main function
`dsmm_fit()` in the non-parametric case. This function can estimate a
Drifting semi-Markov model from a sequence of states (i.e. a character
vector in R). Example data is included in the package, defined in the
DNA sequence `lambda`. Also some parameters need to be specified before
using `dsmm_fit()`, most notably the polynomial *degree* and the model
of our choice. The model is chosen by defining whether the sojourn times
*f* and the transition matrices *p* are drifting or not.

``` r
# Loading the package
library(dsmmR)

# Obtaining the sequence
data("lambda", package = "dsmmR")
sequence <- c(lambda)

# Obtaining the states
states <- sort(unique(sequence))

# Defining the polynomial degree
degree <- 1 # we define a linear evolution in time (state jumps of the embedded Markov chain)

# Defining the model 
f_is_drifting <- TRUE # sojourn time distributions are drifting in time (state jumps of the EMC)
p_is_drifting <- FALSE # transition matrices are not drifting in time (state jumps of the EMC)
# When both f and p are drifting, we have Model 1.

# Fitting the Drifting semi-Markov model
fitted_model <- fit_dsmm(sequence = sequence,
                         states = states,
                         degree = degree,
                         f_is_drifting = f_is_drifting,
                         p_is_drifting = p_is_drifting)
```

For more details about the estimation, consider viewing the extended
documentation through `?fit_dsmm`.

### Defining drifting semi-Markov models

When defining a DSMM object we need to input parameters like the
polynomial degree, the state space, the DSMM size (length of the
embedded Markov chain), the sojourn times *f*, the transition matrices
*p* and more.

For more information, consider the documentation through
`?parametric_dsmm` and `?nonparametric_dsmm`.

### Simulation

After fitting a DSMM (or defining it through `nonparametric_dsmm()` or
`parametric_dsmm()`), we can simulate a sequence from that DSMM. This is
pretty straightforward:

``` r
sim_seq <- simulate(fitted_model)
```

Since we follow an object oriented approach, providing the previous
object `fitted_model` is the only necessary attribute.

For more information, consider the documentation through
`?simulate.dsmm`.

### Drifting semi-Markov kernel

In order to account for the large dimension of the DSM kernel, a
separate function was necessary. You can obtain the DSM kernel through
the command:

``` r
kernel <- get_kernel(fitted_model)
```

The dimensionality of the DSM kernel can be reduced further through the
attributes of the function.

For more information, consider the documentation through `?get_kernel`.

## Further reading

Regarding semi-Markov models, the book [Semi-Markov Chains and Hidden
Semi-Markov Models toward
Applications](https://doi.org/10.1007/978-0-387-73173-5) gives a good
overview of the topic and also combines the flexibility of the
semi-Markov chain with the known advantages of hidden semi-markov
models.

If you are not familiar with Drifting Markov models, they were first
introduced in [Drifting Markov models with Polynomial Drift and
Applications to DNA Sequences](https://doi.org/10.2202/1544-6115.1326),
while a comprehensive overview is provided in [Reliability and Survival
Analysis for Drifting Markov Models: Modeling and
Estimation](https://doi.org/10.1007/s11009-018-9682-8).

### Community Guidelines

For third parties wishing to contribute to the software, or to report
issues or problems about the software, they can do so directly through
the [development github page of the
package](https://github.com/Mavrogiannis-Ioannis/dsmmR).

### Notes

Automated tests are in place in order to aid the user with any false
input made and, furthermore, to ensure that the functions used return
the expected output. Moreover, through strict automated tests, it is
made possible for the user to properly define their own `dsmm` objects
and make use of them with the generic functions of the package.

If you are in need of support, please contact the maintainer at
<mavrogiannis.ioa@gmail.com>.

## References

Barbu, V. S., Limnios, N. (2008). Semi-Markov Chains and Hidden
Semi-Markov Models Toward Applications - Their Use in Reliability and
DNA Analysis. New York: Lecture Notes in Statistics, vol. 191, Springer.

Vergne, N. (2008). Drifting Markov models with Polynomial Drift and
Applications to DNA Sequences. Statistical Applications in Genetics
Molecular Biology 7 (1).

Barbu V. S., Vergne, N. (2019). Reliability and survival analysis for
drifting Markov models: modelling and estimation. Methodology and
Computing in Applied Probability, 21(4), 1407-1429.
