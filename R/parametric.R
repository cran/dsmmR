# '''
#    This file concerns itself with the creation and definition of the
#    parametric drifting semi-Markov model.
# '''
#
#' @title Parametric Drifting semi-Markov model specification
#' @aliases dsmm_parametric parametric
#' @description Creates a parametric model specification for a Drifting
#' semi-Markov model. Returns an object of class
#' \code{(dsmm_parametric, dsmm)}.
#'
#' @param model_size Positive integer that represents the size of
#' the Drifting semi-Markov model \eqn{n}. It is equal to the length of a
#' theoretical embedded Markov chain
#' \eqn{(J_{t})_{t\in \{0,\dots,n\}}}, without the last state.
#' @param states Character vector that represents the state space \eqn{E}.
#'     It has length equal to \eqn{s = |E|}.
#' @param initial_dist Numerical vector of \eqn{s} probabilities, that
#'     represents the initial distribution for each state in the state space
#'     \eqn{E}.
#' @param degree Positive integer that represents the polynomial degree \eqn{d}
#'     for the Drifting semi-Markov model.
#' @param f_is_drifting Logical. Specifies if \eqn{f} is drifting or not.
#' @param p_is_drifting Logical. Specifies if \eqn{p} is drifting or not.
#' @param p_dist Numerical array, that represents the probabilities of the
#'    transition matrix \eqn{p} of the embedded Markov chain
#'    \eqn{(J_{t})_{t\in \{0,\dots,n\}}} (it is defined
#'    the same way in the \link{nonparametric_dsmm} function).
#'    It can be defined in two ways:
#' \itemize{
#'  \item
#'      If \eqn{p} \strong{is not} drifting, it has dimensions of
#'      \eqn{s \times s}.
#'  \item
#'      If \eqn{p} \strong{is} drifting, it has dimensions of
#'      \eqn{s \times s \times (d+1)}
#'  (see more in \emph{Details, Defined Arguments}.)}
#' @param f_dist Character array, that represents the discrete sojourn time
#'     distribution \eqn{f} of our choice.
#'     \code{NA} is allowed for state transitions
#'     that we do not wish to have a sojourn time distribution
#'     (e.g. all state transition to the same state should have \code{NA}
#'     as their value). The list of possible values is:
#'     \code{["unif", "geom", "pois", "dweibull", "nbinom", NA]}.
#'     It can be defined in two ways:
#'     \itemize{
#'     \item
#'          If \eqn{f} \strong{is not} drifting, it has dimensions of
#'          \eqn{s \times s}.
#'     \item
#'         If \eqn{f} \strong{is} drifting, it has dimensions of
#'         \eqn{s \times s \times (d+1)}
#'     (see more in \emph{Details, Defined Arguments}.)}
#' @param f_dist_pars Numerical array, that represents the parameters of the
#'     sojourn time distributions given in \code{f_dist}. \code{NA} is allowed,
#'     in the case that the distribution of our choice does not require
#'     a parameter. It can be defined in two ways:
#' \itemize{
#' \item If \eqn{f} \strong{is not} drifting, it has dimensions of
#'     \eqn{s \times s \times 2}, specifying \strong{two} possible
#'     parameters required for the discrete distributions.
#' \item If \eqn{f} \strong{is} drifting, it has dimensions of
#'     \eqn{s \times s \times 2 \times (d+1)},
#'     specifying \strong{two} possible parameters required for the discrete
#'     distributions, but for every single one of the \eqn{i = 0, \dots, d}
#'     sojourn time distributions \eqn{f_{\frac{i}{d}}} that are required.
#'     (see more in \emph{Details, Defined Arguments}.)}
#'
#' @details
#' \strong{Defined Arguments}
#'
#' For the parametric case, we explicitly define:
#' \enumerate{
#' \item The \emph{transition matrix} of the embedded Markov chain
#'     \eqn{(J_{t})_{t\in \{0,\dots,n\}}}, given in
#'     the attribute \code{p_dist}:
#' \itemize{
#' \item If \eqn{p} \strong{is not drifting}, it contains the values:
#'     \deqn{p(u,v), \forall u, v \in E,}
#'     given in an array with dimensions of \eqn{s \times s},
#'     where the first dimension corresponds to the previous state \eqn{u} and
#'     the second dimension corresponds to the current state \eqn{v}.
#' \item If \eqn{p} \strong{is drifting}, for \eqn{i \in \{ 0,\dots,d \}},
#'     it contains the values:
#'     \deqn{p_{\frac{i}{d}}(u,v), \forall u, v \in E,}
#'     given in an array with dimensions of \eqn{s \times s \times (d + 1)},
#'     where the first and second dimensions are defined as in the non-drifting
#'     case, and the third dimension corresponds to the \eqn{d+1} different
#'     matrices \eqn{p_{\frac{i}{d}}.}
#' }
#' \item The \emph{conditional sojourn time distribution}, given in the
#'     attribute \code{f_dist}:
#' \itemize{
#' \item If \eqn{f} \strong{is not drifting}, it contains the discrete
#'     distribution \emph{names} (as characters or \code{NA}), given in an
#'     array with dimensions of \eqn{s \times s},
#'     where the first dimension corresponds to the previous state \eqn{u},
#'     the second dimension corresponds to the current state \eqn{v}.
#' \item If \eqn{f} \strong{is drifting}, it contains the discrete
#'     distribution \emph{names} (as characters or \code{NA}) given in an
#'     array with dimensions of \eqn{s \times s \times (d + 1)},
#'     where the first and second dimensions are defined as in the
#'     non-drifting case, and the third dimension corresponds to the
#'     \eqn{d+1} different arrays \eqn{f_{\frac{i}{d}}.}
#' }
#' \item The \emph{conditional sojourn time distribution parameters},
#'     given in the attribute \code{f_dist_pars}:
#' \itemize{
#' \item If \eqn{f} \strong{is not drifting}, it contains the
#'     \emph{numerical values} (or \code{NA}) of the corresponding
#'     distributions defined in \code{f_dist}, given in
#'     an array with dimensions of \eqn{s \times s},
#'     where the first dimension corresponds to the previous state \eqn{u},
#'     the second dimension corresponds to the current state \eqn{v}.
#' \item If \eqn{f} \strong{is drifting}, it contains the
#'     \emph{numerical values} (or \code{NA}) of the corresponding
#'     distributions defined in \code{f_dist}, given in an array
#'     with dimensions of \eqn{s \times s \times (d + 1)},
#'     where the first and second dimensions are defined as in the
#'     non-drifting case, and the third dimension corresponds to the
#'     \eqn{d+1} different arrays \eqn{f_{\frac{i}{d}}.}
#' }
#' }
#'
#' \strong{Sojourn time distributions}
#'
#' In this package, the available distributions for the modeling of the
#' conditional sojourn times, of the Drifting semi-Markov model, used through
#' the argument \code{f_dist}, are the following:
#' \itemize{
#' \item Uniform \eqn{(n)}:
#'
#'     \eqn{f(x) = 1/n}, for \eqn{x = 1, 2, \dots, n}, where \eqn{n} is a
#'     positive integer.
#'     This can be specified through the following:
#'     \itemize{
#'     \item \code{f_dist = "unif"}
#'     \item \code{f_dist_pars} = (\eqn{n}, \code{NA})
#'     (\eqn{n} as defined here).
#'     }
#'
#' \item Geometric \eqn{(p)}:
#'
#'     \eqn{f(x) = p (1-p)^{x-1}}, for
#'     \eqn{x = 1, 2, \dots,} where \eqn{p \in (0, 1)}
#'     is the probability of success.
#'     This can be specified through the following:
#'     \itemize{
#'     \item \code{f_dist} = \code{"geom"}
#'     \item \code{f_dist_pars} = (\eqn{p}, \code{NA})
#'     (\eqn{p} as defined here).
#'     }
#' \item Poisson \eqn{(\lambda)}:
#'
#'     \eqn{f(x) = \frac{\lambda^{x-1} exp(-\lambda)}{(x-1)!}}, for
#'     \eqn{x = 1, 2, \dots,} where \eqn{\lambda > 0}.
#'     This can be specified through the following:
#'     \itemize{
#'     \item \code{f_dist} = \code{"pois"}
#'     \item \code{f_dist_pars} = (\eqn{\lambda}, \code{NA})
#'     }
#' \item Negative binomial \eqn{(\alpha, p)}:
#'
#'     \eqn{f(x)=\frac{\Gamma(x+\alpha-1)}{\Gamma(\alpha)(x-1)!}
#'     p^{\alpha}(1-p)^{x-1}}, for \eqn{x = 1, 2,\dots,} where
#'     \eqn{\Gamma} is the Gamma function,
#'     \eqn{\alpha \in (0, +\infty) } is the parameter describing the
#'     target for number of successful trials, or the dispersion parameter
#'     (the shape parameter of the gamma mixing distribution).
#'     \eqn{p} is the probability of success, \eqn{0 < p < 1}.
#'     \itemize{
#'     \item \code{f_dist} = \code{"nbinom"}
#'     \item \code{f_dist_pars} = (\eqn{\alpha, p})
#'     (\eqn{p} as defined here)
#'     }
#' \item Discrete Weibull  of type 1 \eqn{(q, \beta)}:
#'
#'     \eqn{f(x)=q^{(x-1)^{\beta}}-q^{x^{\beta}}}, for \eqn{x=1,2,\dots,} with
#'     \eqn{q \in (0, 1)} is the first parameter (probability)
#'     and \eqn{\beta \in (0, +\infty)} is the second parameter.
#'     This can be specified through the following:
#'     \itemize{
#'     \item \code{f_dist} = \code{"dweibull"}
#'     \item \code{f_dist_pars} = (\eqn{q, \beta})
#'     (\eqn{q} as defined here)
#'     }
#' }
#'
#' From these discrete distributions, by using \code{"dweibull", "nbinom"}
#' we require two parameters. It's for this reason that the attribute
#' \code{f_dist_pars} is an array of dimensions
#' \eqn{s \times s \times 2} if \eqn{f}
#' \strong{is not drifting} or \eqn{s \times s \times 2 \times (d+1)}
#' if \eqn{f} \strong{is drifting}.
#'
#'
#' @return Returns an object of the S3 class \code{dsmm_parametric, dsmm}.
#' It has the following attributes:
#' \itemize{
#' \item \code{dist} : List. Contains 3 arrays, passing down from the arguments:
#'   \itemize{
#'     \item \code{p_drift} or \code{p_notdrift}, corresponding to whether the
#'     defined \eqn{p} transition matrix is drifting or not.
#'     \item \code{f_drift_parametric} or \code{f_notdrift_parametric},
#'     corresponding to whether the
#'     defined \eqn{f} sojourn time distribution is drifting or not.
#'     \item \code{f_drift_parameters} or \code{f_notdrift_parameters},
#'     which are the defined \eqn{f} sojourn time distribution parameters,
#'     depending on whether \eqn{f} is drifting or not.
#'   }
#' \item \code{initial_dist} : Numerical vector. Passing down from the arguments.
#' It contains the initial distribution of the Drifting semi-Markov model.
#' \item \code{states} : Character vector. Passing down from the arguments.
#' It contains the state space \eqn{E}.
#' \item \code{s} : Positive integer. It contains the number of states in the
#' state space, \eqn{s = |E|}, which is given in the attribute \code{states}.
#' \item \code{degree} : Positive integer. Passing down from the arguments.
#' It contains the polynomial degree \eqn{d} considered for the drifting of
#' the model.
#' \item \code{model_size} : Positive integer. Passing down from the arguments.
#' It contains the size of the Drifting semi-Markov model \eqn{n}, which
#' represents the length of the embedded Markov chain
#' \eqn{(J_{t})_{t\in \{0,\dots,n\}}}, without the last state.
#' \item \code{f_is_drifting} : Logical. Passing down from the arguments.
#' Specifies if \eqn{f} is drifting or not.
#' \item \code{p_is_drifting} : Logical. Passing down from the arguments.
#' Specifies if \eqn{p} is drifting or not.
#' \item \code{Model} : Character. Possible values:
#'     \itemize{
#'         \item \code{"Model_1"} : Both \eqn{p} and \eqn{f} are drifting.
#'         \item \code{"Model_2"} : \eqn{p} is drifting and \eqn{f}
#'               is not drifting.
#'         \item \code{"Model_3"} : \eqn{f} is drifting and \eqn{p}
#'               is not drifting.
#'     }
#' \item \code{A_i} : Numerical Matrix. Represents the polynomials
#'     \eqn{A_i(t)} with degree \eqn{d} that are used for solving
#'     the system \eqn{MJ = P}. Used for the methods defined for the
#'     object. Not printed when viewing the object.
#' }
#'
#' @seealso
#' Methods applied to this object: \link{simulate.dsmm}, \link{get_kernel}.
#'
#' For the non-parametric Drifting semi-Markov model specification:
#' \link{nonparametric_dsmm}.
#'
#' For the theoretical background of Drifting semi-Markov models: \link{dsmmR}.
#'
#' @references
#' V. S. Barbu, N. Limnios. (2008). semi-Markov Chains and Hidden semi-Markov
#' Models Toward Applications - Their Use in Reliability and DNA Analysis.
#' New York: Lecture Notes in Statistics, vol. 191, Springer.
#'
#' Vergne, N. (2008). Drifting Markov models with Polynomial Drift and
#' Applications to DNA Sequences. Statistical Applications in Genetics
#' Molecular Biology 7 (1).
#'
#' Barbu V. S., Vergne, N. (2019). Reliability and survival analysis for
#' drifting Markov models: modeling and estimation.
#' Methodology and Computing in Applied Probability, 21(4), 1407-1429.
#'
#' T. Nakagawa and S. Osaki. (1975). The discrete Weibull distribution.
#' IEEE Transactions on Reliability, R-24, 300-301.
#'
#' @export
#'
#' @examples
#' # Setup.
#' # We can define states in a flexible way, including spaces.
#' states <- c("Dollar $", " /1'2'3/ ", " Z E T A ", "O_M_E_G_A")
#' s <- length(states)
#' d <- 1
#'
#' # ===========================================================================
#' # Defining parametric drifting semi-Markov models.
#' # ===========================================================================
#'
#' # ---------------------------------------------------------------------------
#' # Defining the drifting distributions for Model 1.
#' # ---------------------------------------------------------------------------
#'
#' # `p_dist` has dimensions of: (s, s, d + 1).
#' # Sums over v must be 1 for all u and i = 0, ..., d.
#'
#' # First matrix.
#' p_dist_1 <- matrix(c(0,   0.1, 0.4, 0.5,
#'                      0.5, 0,   0.3, 0.2,
#'                      0.3, 0.4, 0,   0.3,
#'                      0.8, 0.1, 0.1, 0),
#'                    ncol = s, byrow = TRUE)
#'
#' # Second matrix.
#' p_dist_2 <- matrix(c(0,   0.3, 0.6, 0.1,
#'                      0.3, 0,   0.4, 0.3,
#'                      0.5, 0.3, 0,   0.2,
#'                      0.2, 0.3, 0.5, 0),
#'                    ncol = s, byrow = TRUE)
#'
#' # get `p_dist` as an array of p_dist_1 and p_dist_2.
#' p_dist_model_1 <- array(c(p_dist_1, p_dist_2), dim = c(s, s, d + 1))
#'
#' # `f_dist` has dimensions of: (s, s, d + 1).
#'
#' # First matrix.
#' f_dist_1 <- matrix(c(NA,         "unif", "dweibull", "nbinom",
#'                      "geom",      NA,    "pois",     "dweibull",
#'                      "dweibull", "pois",  NA,        "geom",
#'                      "pois",      NA,    "geom",      NA),
#'                    nrow = s, ncol = s, byrow = TRUE)
#'
#'
#' # Second matrix.
#' f_dist_2 <- matrix(c(NA,     "pois", "geom", "nbinom",
#'                      "geom",  NA,    "pois", "dweibull",
#'                      "unif", "geom",  NA,    "geom",
#'                      "pois", "pois", "geom",  NA),
#'                    nrow = s, ncol = s, byrow = TRUE)
#'
#'
#' # get `f_dist` as an array of `f_dist_1` and `f_dist_2`
#' f_dist_model_1 <- array(c(f_dist_1, f_dist_2), dim = c(s, s, d + 1))
#'
#'
#' # `f_dist_pars` has dimensions of: (s, s, 2, d + 1).
#' # First array of coefficients, corresponding to `f_dist_1`.
#' # First matrix.
#' f_dist_1_pars_1 <- matrix(c(NA,  5,  0.4, 4,
#'                             0.7, NA, 5,   0.6,
#'                             0.2, 3,  NA,  0.6,
#'                             4,   NA, 0.4, NA),
#'                           nrow = s, ncol = s, byrow = TRUE)
#' # Second matrix.
#' f_dist_1_pars_2 <- matrix(c(NA,  NA, 0.2, 0.6,
#'                             NA,  NA, NA,  0.8,
#'                             0.6, NA, NA,  NA,
#'                             NA,  NA, NA,  NA),
#'                           nrow = s, ncol = s, byrow = TRUE)
#'
#' # Second array of coefficients, corresponding to `f_dist_2`.
#' # First matrix.
#' f_dist_2_pars_1 <- matrix(c(NA,  6,   0.4, 3,
#'                             0.7, NA,  2,   0.5,
#'                             3,   0.6, NA,  0.7,
#'                             6,   0.2, 0.7, NA),
#'                           nrow = s, ncol = s, byrow = TRUE)
#' # Second matrix.
#' f_dist_2_pars_2 <- matrix(c(NA, NA, NA, 0.6,
#'                             NA, NA, NA, 0.8,
#'                             NA, NA, NA, NA,
#'                             NA, NA, NA, NA),
#'                           nrow = s, ncol = s, byrow = TRUE)
#'
#' # Get `f_dist_pars`.
#' f_dist_pars_model_1 <- array(c(f_dist_1_pars_1, f_dist_1_pars_2,
#'                                f_dist_2_pars_1, f_dist_2_pars_2),
#'                              dim = c(s, s, 2, d + 1))
#'
#' # ---------------------------------------------------------------------------
#' # Parametric object for Model 1.
#' # ---------------------------------------------------------------------------
#'
#' obj_par_model_1 <- parametric_dsmm(
#'     model_size = 10000,
#'     states = states,
#'     initial_dist = c(0.8, 0.1, 0.1, 0),
#'     degree = d,
#'     p_dist = p_dist_model_1,
#'     f_dist = f_dist_model_1,
#'     f_dist_pars = f_dist_pars_model_1,
#'     p_is_drifting = TRUE,
#'     f_is_drifting = TRUE
#' )
#'
#' # p drifting array.
#' p_drift <- obj_par_model_1$dist$p_drift
#' p_drift
#'
#' # f distribution.
#' f_dist_drift <- obj_par_model_1$dist$f_drift_parametric
#' f_dist_drift
#'
#' # parameters for the f distribution.
#' f_dist_pars_drift <- obj_par_model_1$dist$f_drift_parameters
#' f_dist_pars_drift
#'
#' # ---------------------------------------------------------------------------
#' # Defining Model 2 - p is drifting, f is not drifting.
#' # ---------------------------------------------------------------------------
#'
#' # `p_dist` has the same dimensions as in Model 1: (s, s, d + 1).
#' p_dist_model_2 <- array(c(p_dist_1, p_dist_2), dim = c(s, s, d + 1))
#'
#' # `f_dist` has dimensions of: (s, s).
#' f_dist_model_2 <- matrix(c( NA,       "pois",  NA,       "nbinom",
#'                             "geom",    NA,    "geom",    "dweibull",
#'                             "unif",   "geom",  NA,       "geom",
#'                             "nbinom", "unif", "dweibull", NA),
#'                          nrow = s, ncol = s, byrow = TRUE)
#'
#' # `f_dist_pars` has dimensions of: (s, s, 2),
#' #  corresponding to `f_dist_model_2`.
#'
#' # First matrix.
#' f_dist_pars_1_model_2 <- matrix(c(NA,  0.2, NA,  3,
#'                                   0.2, NA,  0.2, 0.5,
#'                                   3,   0.4, NA,  0.7,
#'                                   2,   3,   0.7, NA),
#'                                 nrow = s, ncol = s, byrow = TRUE)
#'
#' # Second matrix.
#' f_dist_pars_2_model_2 <- matrix(c(NA,  NA, NA,  0.6,
#'                                   NA,  NA, NA,  0.8,
#'                                   NA,  NA, NA,  NA,
#'                                   0.2, NA, 0.3, NA),
#'                                 nrow = s, ncol = s, byrow = TRUE)
#'
#'
#' # Get `f_dist_pars`.
#' f_dist_pars_model_2 <- array(c(f_dist_pars_1_model_2,
#'                                f_dist_pars_2_model_2),
#'                              dim = c(s, s, 2))
#'
#' # ---------------------------------------------------------------------------
#' # Parametric object for Model 2.
#' # ---------------------------------------------------------------------------
#'
#' obj_par_model_2 <- parametric_dsmm(
#'     model_size = 10000,
#'     states = states,
#'     initial_dist = c(0.8, 0.1, 0.1, 0),
#'     degree = d,
#'     p_dist = p_dist_model_2,
#'     f_dist = f_dist_model_2,
#'     f_dist_pars = f_dist_pars_model_2,
#'     p_is_drifting = TRUE,
#'     f_is_drifting = FALSE
#' )
#'
#' # p drifting array.
#' p_drift <- obj_par_model_2$dist$p_drift
#' p_drift
#'
#' # f distribution.
#' f_dist_notdrift <- obj_par_model_2$dist$f_notdrift_parametric
#' f_dist_notdrift
#'
#' # parameters for the f distribution.
#' f_dist_pars_notdrift <- obj_par_model_2$dist$f_notdrift_parameters
#' f_dist_pars_notdrift
#'
#' # ---------------------------------------------------------------------------
#' # Defining Model 3 - f is drifting, p is not drifting.
#' # ---------------------------------------------------------------------------
#'
#' # `p_dist` has dimensions of: (s, s).
#' p_dist_model_3 <- matrix(c(0,   0.1,  0.3,  0.6,
#'                            0.4, 0,    0.1,  0.5,
#'                            0.4, 0.3,  0,    0.3,
#'                            0.9, 0.01, 0.09, 0),
#'                          ncol = s, byrow = TRUE)
#'
#' # `f_dist` has the same dimensions as in Model 1: (s, s, d + 1).
#' f_dist_model_3 <- array(c(f_dist_1, f_dist_2), dim = c(s, s, d + 1))
#'
#'
#' # `f_dist_pars` has the same dimensions as in Model 1: (s, s, 2, d + 1).
#' f_dist_pars_model_3 <- array(c(f_dist_1_pars_1, f_dist_1_pars_2,
#'                                f_dist_2_pars_1, f_dist_2_pars_2),
#'                              dim = c(s, s, 2, d + 1))
#'
#' # ---------------------------------------------------------------------------
#' # Parametric object for Model 3.
#' # ---------------------------------------------------------------------------
#'
#' obj_par_model_3 <- parametric_dsmm(
#'     model_size = 10000,
#'     states = states,
#'     initial_dist = c(0.3, 0.2, 0.2, 0.3),
#'     degree = d,
#'     p_dist = p_dist_model_3,
#'     f_dist = f_dist_model_3,
#'     f_dist_pars = f_dist_pars_model_3,
#'     p_is_drifting = FALSE,
#'     f_is_drifting = TRUE
#' )
#'
#' # p drifting array.
#' p_notdrift <- obj_par_model_3$dist$p_notdrift
#' p_notdrift
#'
#' # f distribution.
#' f_dist_drift <- obj_par_model_3$dist$f_drift_parametric
#' f_dist_drift
#'
#' # parameters for the f distribution.
#' f_dist_pars_drift <- obj_par_model_3$dist$f_drift_parameters
#' f_dist_pars_drift
#'
#' # ===========================================================================
#' # Parametric estimation using methods corresponding to an object
#' #     which inherits from the class `dsmm_parametric`.
#' # ===========================================================================
#' ### Comments
#' ### 1.  Using a larger `klim` and a larger `model_size` will increase the
#' ###     accuracy of the model, with the need of larger memory requirements
#' ###     and computational cost.
#' ### 2.  For the parametric estimation it is recommended to use a common set
#' ###     of distributions while only the parameters are drifting. This results
#' ###     in higher accuracy.
#'
#'
#' # ---------------------------------------------------------------------------
#' # Defining the distributions for Model 1 - both p and f are drifting.
#' # ---------------------------------------------------------------------------
#'
#' # `p_dist` has dimensions of: (s, s, d + 1).
#' # First matrix.
#' p_dist_1 <- matrix(c(0,   0.2, 0.4, 0.4,
#'                      0.5, 0,   0.3, 0.2,
#'                      0.3, 0.4, 0,   0.3,
#'                      0.5, 0.3, 0.2, 0),
#'                    ncol = s, byrow = TRUE)
#'
#' # Second matrix.
#' p_dist_2 <- matrix(c(0,   0.3, 0.5, 0.2,
#'                      0.3, 0,   0.4, 0.3,
#'                      0.5, 0.3, 0,   0.2,
#'                      0.2, 0.4, 0.4, 0),
#'                    ncol = s, byrow = TRUE)
#'
#' # get `p_dist` as an array of p_dist_1 and p_dist_2.
#' p_dist_model_1 <- array(c(p_dist_1, p_dist_2), dim = c(s, s, d + 1))
#'
#' # `f_dist` has dimensions of: (s, s, d + 1).
#' # We will use the same sojourn time distributions.
#' f_dist_1 <- matrix(c( NA,        "unif",   "dweibull", "nbinom",
#'                      "geom",      NA,      "pois",     "dweibull",
#'                      "dweibull", "pois",    NA,        "geom",
#'                      "pois",     'nbinom', "geom",      NA),
#'                    nrow = s, ncol = s, byrow = TRUE)
#'
#' # get `f_dist`
#' f_dist_model_1 <- array(f_dist_1, dim = c(s, s, d + 1))
#'
#' # `f_dist_pars` has dimensions of: (s, s, 2, d + 1).
#' # First array of coefficients, corresponding to `f_dist_1`.
#' # First matrix.
#' f_dist_1_pars_1 <- matrix(c(NA,  7, 0.4, 4,
#'                             0.7, NA, 5,   0.6,
#'                             0.2, 3,  NA,  0.6,
#'                             4,   4,  0.4, NA),
#'                           nrow = s, ncol = s, byrow = TRUE)
#' # Second matrix.
#' f_dist_1_pars_2 <- matrix(c(NA,  NA,  0.2, 0.6,
#'                             NA,  NA,  NA,  0.8,
#'                             0.6, NA,  NA,  NA,
#'                             NA,  0.3, NA,  NA),
#'                           nrow = s, ncol = s, byrow = TRUE)
#'
#' # Second array of coefficients, corresponding to `f_dist_2`.
#' # First matrix.
#' f_dist_2_pars_1 <- matrix(c(NA,  6,  0.5, 3,
#'                             0.5, NA, 4,   0.5,
#'                             0.4, 5,  NA,  0.7,
#'                             6,   5,  0.7, NA),
#'                           nrow = s, ncol = s, byrow = TRUE)
#' # Second matrix.
#' f_dist_2_pars_2 <- matrix(c(NA,  NA,  0.4, 0.5,
#'                             NA,  NA,  NA,  0.6,
#'                             0.5, NA,  NA,  NA,
#'                             NA,  0.4, NA,  NA),
#'                           nrow = s, ncol = s, byrow = TRUE)
#' # Get `f_dist_pars`.
#' f_dist_pars_model_1 <- array(c(f_dist_1_pars_1, f_dist_1_pars_2,
#'                                f_dist_2_pars_1, f_dist_2_pars_2),
#'                              dim = c(s, s, 2, d + 1))
#'
#' # ---------------------------------------------------------------------------
#' # Defining the parametric object for Model 1.
#' # ---------------------------------------------------------------------------
#'
#' obj_par_model_1 <- parametric_dsmm(
#'     model_size = 4000,
#'     states = states,
#'     initial_dist = c(0.8, 0.1, 0.1, 0),
#'     degree = d,
#'     p_dist = p_dist_model_1,
#'     f_dist = f_dist_model_1,
#'     f_dist_pars = f_dist_pars_model_1,
#'     p_is_drifting = TRUE,
#'     f_is_drifting = TRUE
#' )
#'
#' cat("The object has class of (",
#'     paste0(class(obj_par_model_1),
#'            collapse = ', '), ").")
#'
#'
#' # ---------------------------------------------------------------------------
#' # Generating a sequence from the parametric object.
#' # ---------------------------------------------------------------------------
#'
#' # A larger klim will increase accuracy and also the computation cost.
#' klim <- 20
#' sim_seq <- simulate(obj_par_model_1, klim = klim, seed = 1)
#'
#'
#' # ---------------------------------------------------------------------------
#' # Fitting the generated sequence under the same distributions.
#' # ---------------------------------------------------------------------------
#'
#' fit_par_model1 <- fit_dsmm(sequence = sim_seq,
#'                            states = states,
#'                            degree = d,
#'                            f_is_drifting = TRUE,
#'                            p_is_drifting = TRUE,
#'                            estimation = 'parametric',
#'                            f_dist = f_dist_model_1)
#'
#' cat("The object has class of (",
#'     paste0(class(fit_par_model1),
#'            collapse = ', '), ").")
#'
#' cat("\nThe estimated parameters are:\n")
#' fit_par_model1$dist$f_drift_parameters
#'
#'
parametric_dsmm <- function(model_size,
                            states,
                            initial_dist,
                            degree,
                            f_is_drifting,
                            p_is_drifting,
                            p_dist,
                            f_dist,
                            f_dist_pars) {
    # Check for the validity of parameters given.
    if (missing(model_size)) {
        stop("\nPlease provide a positive integer for the",
             " `model_size` parameter.")
    }
    stopifnot(is_integer(model_size))
    model_size <- as.integer(model_size)
    # State Space
    if (missing(states)) {
        stop("\nPlease provide a character vector for the `states`",
             " of the State Space.")
    }
    stopifnot(valid_states(states))
    s <- length(states)
    # Initial Distribution
    if (missing(initial_dist)) {
        stop("\nPlease provide a numeric vector for the ",
             "initial distribution defined in `initial_dist`.")
    }
    stopifnot(valid_initial_dist(initial_dist = initial_dist, s = s))
    names(initial_dist) <- states
    # degree
    if (missing(degree)) {
        stop("\nPlease provide the polynomial `degree` as an integer.")
    }
    stopifnot(valid_degree(degree))
    degree <- as.integer(degree)
    # f_is_drifting
    if (missing(f_is_drifting)) {
        stop("\nPlease provide whether the sojourn time distribution f ",
             "is drifting through the logical parameter `f_is_drifting`.")
    } else if (!is_logical(f_is_drifting)) {
        stop("\nThe logical parameter `f_is_drifting` should be ",
             "either TRUE or FALSE.")
    }
    # p_is_drifting
    if (missing(p_is_drifting)) {
        stop("\nPlease provide whether the transition matrix p is drifting",
             " through the logical parameter `p_is_drifting`.")
    } else if (!is_logical(p_is_drifting)) {
        stop("\nThe logical parameter `p_is_drifting` should be ",
             "either TRUE or FALSE.")
    }
    # p_dist
    if (missing(p_dist)) {
        stop("\nPlease provide the `p_dist` array.")
    }
    # f_dist
    if (missing(f_dist)) {
        stop("\nPlease provide the `f_dist` array.")
    }
    # f_dist_pars
    if (missing(f_dist_pars)) {
        stop("\nPlease provide the `f_dist_pars` array.")
    }
    stopifnot(valid_model(p_is_drifting = p_is_drifting,
                          f_is_drifting = f_is_drifting),
              valid_p_dist(p_dist = p_dist,
                           s = s, degree = degree,
                           p_is_drifting = p_is_drifting,
                           states = states),
              valid_fdist_parametric(fdist = f_dist,
                                     states = states, s = s,
                                     degree = degree,
                                     params = f_dist_pars,
                                     f_is_drifting = f_is_drifting)
    )
    # Assign names for p_dist, f_dist, f_dist_pars.
    # Empty the names first so the process works as intended.
    dimnames(p_dist) <- NULL
    dimnames(f_dist) <- NULL
    dimnames(f_dist_pars) <- NULL
    dimnames(p_dist) <- list(as.list(states), as.list(states))
    dimnames(f_dist) <- list(as.list(states), as.list(states))
    # 2 matrices for the parameters.
    dimnames(f_dist_pars) <- list(as.list(states), as.list(states),
                                                as.list(1:2))
    if (p_is_drifting) {
        dimnames(p_dist)[[3]] <- as.list(names_i_d(degree, 'p'))
    }
    if (f_is_drifting) {
        dimnames(f_dist)[[3]] <- as.list(names_i_d(degree, 'f'))
        dimnames(f_dist_pars)[[4]] <- as.list(names_i_d(degree, 'fpars'))
    }
    # Get `dist` with regards to `model`.
    model <- get_model(p_is_drifting = p_is_drifting,
                       f_is_drifting = f_is_drifting)
    if (model == "Model_1") {
        dist <- list('p_drift' = p_dist, 'f_drift_parametric' = f_dist,
                     'f_drift_parameters' = f_dist_pars)
    } else if (model == "Model_2") {
        dist <- list('p_drift' = p_dist, 'f_notdrift_parametric' = f_dist,
                     'f_notdrift_parameters' = f_dist_pars)
    } else if (model == "Model_3") {
        dist <- list('p_notdrift' = p_dist, 'f_drift_parametric' = f_dist,
                     'f_drift_parameters' = f_dist_pars)
    }
    # Assign the values to the object.
    obj <- list(
        "dist" = dist,
        "initial_dist" = initial_dist,
        "states" = states,
        "s" =  s,
        "degree" = degree,
        "model_size" = model_size,
        "f_is_drifting" = f_is_drifting,
        "p_is_drifting" = p_is_drifting,
        'Model' = model,
        "A_i" = get_A_i(degree, model_size) # Not shown through `print`.
    )
    class(obj) <- c("dsmm_parametric", "dsmm")
    return(obj)
}

