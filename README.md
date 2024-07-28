README
================

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- Use `build_readme()` to edit the .md file from the .rmd file -->
<!-- badges: start -->

[![CRAN status](https://www.r-pkg.org/badges/version/dsmmR)](https://CRAN.R-project.org/package=dsmmR)
[![R-CMD-check](https://github.com/Mavrogiannis-Ioannis/dsmmR/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/Mavrogiannis-Ioannis/dsmmR/actions/workflows/R-CMD-check.yaml)
[![CRAN logs](https://cranlogs.r-pkg.org/badges/dsmmR)](https://CRAN.R-project.org/package=dsmmR)
[![Codecov test coverage](https://codecov.io/gh/Mavrogiannis-Ioannis/dsmmR/branch/master/graph/badge.svg)](https://app.codecov.io/gh/Mavrogiannis-Ioannis/dsmmR?branch=master)
[![Contributor Covenant](https://img.shields.io/badge/Contributor%20Covenant-2.1-4baaaa.svg)](https://usethis.r-lib.org/CODE_OF_CONDUCT.html)
<!-- badges: end -->

# Developer Version

## 1.0.5

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
complete documentation of the package on the [official CRAN page](https://CRAN.R-project.org/package=dsmmR).

### Estimation

The easiest way to use **dsmmR** is through the main function
`dsmm_fit()` in the non-parametric case. This function can estimate a
Drifting semi-Markov model from a sequence of states (i.e. a character
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
# When f is drifting and p is not drifting, we have Model 3.

# Fitting the drifting semi-Markov model on the sequence.
fitted_model <- fit_dsmm(sequence = sequence,
                         states = states,
                         degree = degree,
                         f_is_drifting = f_is_drifting,
                         p_is_drifting = p_is_drifting)
```

For more details about the estimation, consider viewing the extended
documentation through `?fit_dsmm`.

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

In order to account for the dimension of the DSM kernel, a separate
function was necessary. You can obtain the DSM kernel through the
command:

``` r
kernel <- get_kernel(fitted_model)
```

The dimensionality of the DSM kernel can be reduced further through the
attributes of the function.

For more information, consider the documentation through `?get_kernel`.

### Defining drifting semi-Markov models

We can put together all the previous concepts in the showcase of
parametric estimation. First, we will define the drifting transition
matrices and the drifting sojourn time distributions. Then, we will
create a `dsmm_parametric` object, we will simulate a sequence from it
and then finally we will estimate a drifting semi-Markov model from that
simulated sequence.

For more information, consider the documentation through
`?parametric_dsmm` and `?nonparametric_dsmm`.

First of all we load the package,

``` r
library(dsmmR)
```

and then we define the states and we set the degree equal to 1.

``` r
states <- c("a", "b", "c")
s <- length(states)
degree <- 1
```

Since degree is equal to 1, we then define the 2 drifting transition
matrices:

``` r
p_dist_1 <- matrix(c(0,   0.4,  0.6,
                     0.5, 0,    0.5,
                     0.3, 0.7,  0   ), ncol = s, byrow = TRUE)
p_dist_2 <- matrix(c(0,   0.55, 0.45,
                     0.25, 0,   0.75,
                     0.5, 0.5,  0   ), ncol = s, byrow = TRUE)
p_dist <- array(c(p_dist_1, p_dist_2), dim = c(s, s, degree + 1))
```

Let us also consider the case where only the parameters of the
distributions modeling the sojourn times are drifting across the
sequence. Note that distributions like the Negative Binomial and the
Discrete Weibull require two parameters, which we define in two matrices
for each distribution.

``` r
f_dist_1 <- matrix(c(NA,   "nbinom",   "unif",
                   "geom",  NA,        "pois",
                   "pois", "dweibull",  NA   ), nrow = s, ncol = s, byrow = TRUE)
f_dist_1_pars_1 <- matrix(c(NA,  4,   3,
                            0.7, NA,  5,
                            3,   0.6, NA), nrow = s, ncol = s, byrow = TRUE)
f_dist_1_pars_2 <- matrix(c(NA,  0.5, NA,
                            NA,  NA,  NA,
                            NA,  0.8, NA), nrow = s, ncol = s, byrow = TRUE)
f_dist_2 <- f_dist_1 
f_dist_2_pars_1 <- matrix(c(NA,  3,   5,
                            0.3, NA,  2,
                            5,   0.3, NA), nrow = s, ncol = s, byrow = TRUE)
f_dist_2_pars_2 <- matrix(c(NA,  0.4, NA,
                            NA,  NA,  NA,
                            NA,  0.5, NA), nrow = s, ncol = s, byrow = TRUE)

f_dist <- array(c(f_dist_1, f_dist_2), dim = c(s, s, degree + 1))
f_dist_pars <- array(c(f_dist_1_pars_1, f_dist_1_pars_2,
                       f_dist_2_pars_1, f_dist_2_pars_2), 
                     dim = c(s, s, 2, degree + 1))
```

Then, defining a `dsmm_parametric` object is done simply through the
function `parametric_dsmm()`:

``` r
dsmm_model <- parametric_dsmm(
    model_size = 10000,
    states = states,
    initial_dist = c(0.6, 0.3, 0.1),
    degree = degree,
    p_dist = p_dist,
    f_dist = f_dist,
    f_dist_pars = f_dist_pars,
    p_is_drifting = TRUE,
    f_is_drifting = TRUE
)
```

We can then simulate a sequence from this parametric object like-so:

``` r
sim_seq <- simulate(dsmm_model, klim = 30, seed = 1)
```

To fit this sequence with a drifting semi-Markov model, one can use:

``` r
fitted_model <- fit_dsmm(sequence = sim_seq,
                         states = states,
                         degree = degree,
                         f_is_drifting = TRUE,
                         p_is_drifting = TRUE,
                         estimation = 'parametric',
                         f_dist = f_dist)
```

Finally, the drifting transition matrix is estimated as:

``` r
print(fitted_model$dist$p_drift, digits = 2)
```

with output:

``` r
, , p_0

     a    b    c
a 0.00 0.40 0.60
b 0.51 0.00 0.49
c 0.27 0.73 0.00

, , p_1

     a    b    c
a 0.00 0.54 0.46
b 0.23 0.00 0.77
c 0.51 0.49 0.00
```

and the parameters for the drifting sojourn time distributions are:

``` r
print(fitted_model$dist$f_drift_parameters, digits = 2)
```

with output:

``` r
, , 1, fpars_0

     a    b   c
a   NA 3.66 3.0
b 0.65   NA 4.8
c 3.09 0.62  NA

, , 2, fpars_0

   a    b  c
a NA 0.46 NA
b NA   NA NA
c NA 0.84 NA

, , 1, fpars_1

     a    b   c
a   NA 2.74 5.0
b 0.31   NA 2.1
c 5.02 0.29  NA

, , 2, fpars_1

   a    b  c
a NA 0.38 NA
b NA   NA NA
c NA 0.50 NA
```

## Further reading

Regarding semi-Markov models, the book [Semi-Markov Chains and Hidden Semi-Markov Models toward Applications](https://doi.org/10.1007/978-0-387-73173-5) gives a good
overview of the topic and also combines the flexibility of the
semi-Markov chain with the known advantages of hidden semi-markov
models.

If you are not familiar with Drifting Markov models, they were first
introduced in [Drifting Markov models with Polynomial Drift and Applications to DNA Sequences](https://doi.org/10.2202/1544-6115.1326),
while a comprehensive overview is provided in [Reliability and Survival Analysis for Drifting Markov Models: Modeling and Estimation](https://doi.org/10.1007/s11009-018-9682-8).

### Community Guidelines

For third parties wishing to contribute to the software, or to report
issues or problems about the software, they can do so directly through
the [development github page of the package](https://github.com/Mavrogiannis-Ioannis/dsmmR).

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
DNA Analysis. New York: Lecture Notes in Statistics, vol. 191, Springer.

Vergne, N. (2008). Drifting Markov models with Polynomial Drift and
Applications to DNA Sequences. Statistical Applications in Genetics
Molecular Biology 7 (1).

Barbu V. S., Vergne, N. (2019). Reliability and survival analysis for
drifting Markov models: modelling and estimation. Methodology and
Computing in Applied Probability, 21(4), 1407-1429.

## Acknowledgements

We acknowledge the DATALAB Project
<https://lmrs-num.math.cnrs.fr/projet-datalab.html> (financed by the
European Union with the European Regional Development fund (ERDF) and by
the Normandy Region) and the HSMM-INCA Project (financed by the French
Agence Nationale de la Recherche (ANR) under grant ANR-21-CE40-0005).
