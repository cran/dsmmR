# '''
#    1. This file contains the generic functions that are used for the objects
#    defined in this package.
#    These objects are of the class `dsmm`, which acts like the parent class,
#    and then the 4 child classes:
#       `dsmm_fit_nonparametric`, `dsmm_fit_parametric` and
#       `dsmm_nonparametric` and `dsmm_parametric`.
#    The `dsmm` class is only used if there is no need to classify a function
#    for the three child classes.
#    2. It is worth noting that, from the following functions, ONLY the
#    generics have documentation available to the user,
#    apart from the functions :
#       `get_kernel()`, `simulate.dsmm()`, `is.dsmm()`,
#       `is.dsmm_fit()`, `is.dsmm_nonparametric()` and `is.dsmm_parametric()`,
#    which are all available to the user.
# '''


# ______________________________________________________________________________
# Checking the validity of the attributes passed into the functions.
# ______________________________________________________________________________
check_attributes <- function(obj) UseMethod('check_attributes', obj)

check_attributes.dsmm_fit_nonparametric <- function(obj) {
    # Check whether `f_is_drifting`, `p_is_drifting` are correctly used.
    if (!is_logical(f_is_drifting <- obj$f_is_drifting)) {
        stop("\nThe logical parameter `f_is_drifting` ",
             "should be either TRUE or FALSE.")
    }
    if (!is_logical(p_is_drifting <- obj$p_is_drifting)) {
        stop("\nThe logical parameter `p_is_drifting` ",
             "should be either TRUE or FALSE.")
    }
    # # numerical_est
    # if (!is_logical(numerical_est <- obj$numerical_est)) {
    #     stop("\nThe logical parameter `numerical_est` should ",
    #          "be either TRUE or FALSE.")
    # }
    # Check names of `dist`, in order for `obj$dist[[ i ]]`
    # to work in the next section of `stopifnot`.
    pname <- if (p_is_drifting) 'p_drift' else 'p_notdrift'
    fname <- if (f_is_drifting) 'f_drift' else 'f_notdrift'
    if (!identical((names_tmp <- names(obj$dist)),
                   (names_real <- c(pname, fname)))) {
        stop('\n`', paste0(substitute(obj), '$dist`'),
             ' should have the named distributions in order, as in: \n',
             paste0(1:2, ". ", names_real, collapse = ', '),
             ', when it has :\n',
             paste0(1:length(names_tmp), ". ",
                    names_tmp, collapse = ', '))
    }
    stopifnot(
        valid_model(p_is_drifting = p_is_drifting,
                    f_is_drifting = f_is_drifting
                    # , numerical_est = numerical_est
        ),
        valid_seq(seq = (seq <- obj$seq)),
        valid_model_size(model_size = (model_size <- obj$model_size),
                         length_seq = (l <- length(seq))),
        valid_soj_times(soj_times = (soj_times <- obj$soj_times),
                        length_seq = l),
        valid_k_max(k_max = (k_max <- obj$k_max), soj_times = soj_times),
        valid_states(states = (states <- obj$states)),
        valid_length_states(s = (s <- obj$s), states = states),
        valid_initial_dist(obj$initial_dist, s),
        valid_degree(degree = (degree <- obj$degree), model_size = model_size),
        valid_estimation_fit(estimation = obj$estimation,
                             degree = degree,
                             states = states,
                             s = s,
                             k_max = k_max,
                             f_is_drifting = f_is_drifting,
                             p_is_drifting = p_is_drifting,
                             pdist = obj$dist[[1]],
                             fdist = obj$dist[[2]]
        )
    )
    # Check names of the list.
    names_fit <- c("dist", "seq", "soj_times", "initial_dist", "states",
                   "s", "degree", "k_max", "model_size", "f_is_drifting",
                   "p_is_drifting", "Model", "estimation", "A_i", "J_i")
    if (!identical((nobj <- names(obj)), names_fit)) {
        msg <- c("\nThe attributes of the object defined are not",
                 " the same as for a parametric object defined ",
                 "through `nonparametric_dsmm()`.\n")
        extranames <- nobj[which(!nobj %in% names_fit)]
        if (length(extranames) > 0) {
            msg <- c(msg, paste0("\nThey have the extra names of:\n", '"',
                                 paste(extranames, collapse = '", "'), '"'))
        }
        lackingnames <- names_fit[which(!names_fit %in% nobj)]
        if (length(lackingnames) > 0) {
            msg <- c(msg, paste0("\n\nThey are lacking the names of:\n", '"',
                                 paste(lackingnames, collapse = '", "'), '"'))
        }
        message(msg)
        return(FALSE)
    }
    TRUE
}

check_attributes.dsmm_fit_parametric <- function(obj) {
    # Check whether `f_is_drifting`, `p_is_drifting` are correctly used.
    if (!is_logical(f_is_drifting <- obj$f_is_drifting)) {
        stop("\nThe logical parameter `f_is_drifting` ",
             "should be either TRUE or FALSE.")
    }
    if (!is_logical(p_is_drifting <- obj$p_is_drifting)) {
        stop("\nThe logical parameter `p_is_drifting` ",
             "should be either TRUE or FALSE.")
    }
    # # numerical_est
    # if (!is_logical(numerical_est <- obj$numerical_est)) {
    #     stop("\nThe logical parameter `numerical_est` should ",
    #          "be either TRUE or FALSE.")
    # }
    # Check names of `dist`, in order for `obj$dist[[ i ]]`
    # to work in the next section of `stopifnot`.
    pname <- if (p_is_drifting) 'p_drift' else 'p_notdrift'
    fname <-
        if (f_is_drifting) 'f_drift_parametric' else 'f_notdrift_parametric'
    fparname <-
        if (f_is_drifting) 'f_drift_parameters' else 'f_notdrift_parameters'
    if (!identical((names_tmp <- names(obj$dist)),
                   (names_real <- c(pname, fname, fparname)))) {
        stop('\n`', paste0(substitute(obj), '$dist`'),
             ' should have the named distributions in order, as in: \n',
             paste0(1:2, ". ", names_real, collapse = ', '),
             ', when it has :\n',
             paste0(1:length(names_real), ". ",
                    names_tmp, collapse = ', '))
    }
    stopifnot(
        valid_model(p_is_drifting = p_is_drifting,
                    f_is_drifting = f_is_drifting
                    # , numerical_est = numerical_est
        ),
        valid_seq(seq = (seq <- obj$seq)),
        valid_model_size(model_size = (model_size <- obj$model_size),
                         length_seq = (l <- length(seq))),
        valid_soj_times(soj_times = (soj_times <- obj$soj_times),
                        length_seq = l),
        valid_k_max(k_max = (k_max <- obj$k_max), soj_times = soj_times),
        valid_states(states = (states <- obj$states)),
        valid_length_states(s = (s <- obj$s), states = states),
        valid_initial_dist(obj$initial_dist, s),
        valid_degree(degree = (degree <- obj$degree), model_size = model_size),
        valid_estimation_fit(estimation = obj$estimation,
                             degree = degree,
                             states = states,
                             s = s,
                             k_max = k_max,
                             f_is_drifting = f_is_drifting,
                             p_is_drifting = p_is_drifting,
                             pdist = obj$dist[[1]],
                             fdist = obj$dist[[2]]
        )
    )
    # Check names of the list.
    names_fit <- c("dist", "seq", "soj_times", "initial_dist", "states",
                   "s", "degree", "k_max", "model_size", "f_is_drifting",
                   "p_is_drifting", "Model", "estimation", "A_i", "J_i")
    if (!identical((nobj <- names(obj)), names_fit)) {
        msg <- c("\nThe attributes of the object defined are not",
                 " the same as for a parametric object defined ",
                 "through `nonparametric_dsmm()`.\n")
        extranames <- nobj[which(!nobj %in% names_fit)]
        if (length(extranames) > 0) {
            msg <- c(msg, paste0("\nThey have the extra names of:\n", '"',
                                 paste(extranames, collapse = '", "'), '"'))
        }
        lackingnames <- names_fit[which(!names_fit %in% nobj)]
        if (length(lackingnames) > 0) {
            msg <- c(msg, paste0("\n\nThey are lacking the names of:\n", '"',
                                 paste(lackingnames, collapse = '", "'), '"'))
        }
        message(msg)
        return(FALSE)
    }
    TRUE
}

check_attributes.dsmm_nonparametric <- function(obj) {
    # Check whether `f_is_drifting`, `p_is_drifting` are correctly used.
    if (!is_logical(f_is_drifting <- obj$f_is_drifting)) {
        stop("\nThe logical parameter `f_is_drifting` should be ",
             "either TRUE or FALSE.")
    }
    if (!is_logical(p_is_drifting <- obj$p_is_drifting)) {
        stop("\nThe logical parameter `p_is_drifting` should be ",
             "either TRUE or FALSE.")
    }
    # Check names of `dist`, in order for `obj$dist[[ i ]]`
    # to work in the next section of `stopifnot`.
    pname <- if (p_is_drifting) 'p_drift' else 'p_notdrift'
    fname <- if (f_is_drifting) 'f_drift' else 'f_notdrift'
    if (any((names_tmp <- names(obj$dist)) !=
            (names_real <- c(pname, fname)))) {
        stop('\n`', paste0(substitute(obj), '$dist`'),
             ' should have the named distributions in order, as in: \n',
             paste0(1:2, ". ", names_real, collapse = ', '),
             ', when it has :\n',
             paste0(1:2, ". ", names_tmp, collapse = ', '))
    }
    # Check whether the object attributes are correctly used.()
    stopifnot(
        is_integer(model_size <- (model_size <- obj$model_size)),
        is_integer(k_max <- obj$k_max),
        valid_states(states <- obj$states),
        valid_length_states(s <- obj$s, states),
        valid_degree(degree <- obj$degree, model_size = model_size),
        valid_model(p_is_drifting, f_is_drifting),
        valid_p_dist(p_dist = obj$dist[[1]],
                     states = states,
                     s = s, degree = degree,
                     p_is_drifting = p_is_drifting),
        valid_fdist_nonparametric(f_dist = obj$dist[[2]],
                                  states = states, s = s,
                                  degree = degree,
                                  f_is_drifting = f_is_drifting,
                                  k_max = k_max)
    )
    # Check names of the list.
    names_nonpar <- c(
        "dist", "initial_dist", "states", "s", "degree", "k_max",
        "model_size", "f_is_drifting", "p_is_drifting", 'Model', "A_i")
    if (!identical((nobj <- names(obj)), names_nonpar)) {
        msg <- c("\nThe attributes of the object defined are not",
                 " the same as for a parametric object defined ",
                 "through `nonparametric_dsmm()`.\n")
        extranames <- nobj[which(!nobj %in% names_nonpar)]
        if (length(extranames) > 0) {
            msg <- c(msg, paste0("\nThey have the extra names of:\n", '"',
                                 paste(extranames, collapse = '", "'), '"'))
        }
        lackingnames <- names_nonpar[which(!names_nonpar %in% nobj)]
        if (length(lackingnames) > 0) {
            msg <- c(msg, paste0("\n\nThey are lacking the names of:\n", '"',
                                 paste(lackingnames, collapse = '", "'), '"'))
        }
        message(msg)
        return(FALSE)
    }
    TRUE
}

check_attributes.dsmm_parametric <- function(obj) {
    if (!is_logical(f_is_drifting <- obj$f_is_drifting)) {
        stop("\nThe logical parameter `f_is_drifting` should be either ",
             "TRUE or FALSE.")
    }
    if (!is_logical(p_is_drifting <- obj$p_is_drifting)) {
        stop("\nThe logical parameter `p_is_drifting` should be either ",
             "TRUE or FALSE.")
    }
    # Check names of dist, in order for `obj$dist[[ i ]]` to work in
    # the next section of `stopifnot`.
    pname <- if (p_is_drifting) 'p_drift' else 'p_notdrift'
    fname <- if (f_is_drifting) 'f_drift_parametric' else 'f_notdrift'
    fparname <-
        if (f_is_drifting) 'f_drift_parameters' else 'f_notdrift_parameters'
    if (any((names_tmp <- names(obj$dist)) !=
            (names_real <- c(pname, fname, fparname)))) {
        stop('\n`', paste0(substitute(obj), '$dist`'),
             ' should have the named distributions in order, as in: \n',
             paste0(1:2, ". ", names_real, collapse = ', '),
             ', when it has :\n',
             paste0(1:2, ". ", names_tmp, collapse = ', '))
    }
    # Check whether `numerical_est`, `f_is_drifting` and `p_is_drifting`
    # are correctly used.
    stopifnot(is_integer(obj = (model_size <- obj$model_size)),
              valid_model(p_is_drifting, f_is_drifting),
              valid_states(states = (states <- obj$states)),
              valid_length_states(s = (s <- obj$s), states),
              valid_initial_dist(initial_dist = obj$initial_dist,
                                 s = s),
              valid_degree(degree = (degree <- obj$degree),
                           model_size = model_size),
              valid_p_dist(p_dist = obj$dist[[1]],
                           states = states,
                           s = s, degree = degree,
                           p_is_drifting = p_is_drifting),
              valid_fdist_parametric(fdist = obj$dist[[2]],
                                     params = obj$dist[[3]],
                                     s = s, degree = degree,
                                     states = states,
                                     f_is_drifting = f_is_drifting)
    )
    # Check names of the list.
    names_par <- c(
        "dist", "initial_dist", "states", "s", "degree", "model_size",
        "f_is_drifting", "p_is_drifting", 'Model', "A_i")
    if (!identical((nobj <- names(obj)), names_par)) {
        msg <- c("\nThe attributes of the object defined are not",
                 " the same as for a parametric object defined ",
                 "through `parametric_dsmm()`.\n")
        extranames <- nobj[which(!nobj %in% names_par)]
        if (length(extranames) > 0) {
            msg <- c(msg, paste0("\nThey have the extra names of:\n", '"',
                                 paste(extranames, collapse = '", "'), '"'))
        }
        lackingnames <- names_par[which(!names_par %in% nobj)]
        if (length(lackingnames) > 0) {
            msg <- c(msg, paste0("\nThey are lacking the names of:\n",
                                 paste(lackingnames, collapse = '", "'), '"'))
        }
        message(msg)
        return(FALSE)
    }
    TRUE
}


# ______________________________________________________________________________
# Checking the class of an object and if it satisfies the necessary conditions.
# ______________________________________________________________________________

#' @title Check if an object has a valid \code{dsmm} class
#' @description Checks for the validity of the specified attributes  and the
#'    inheritance of the S3 class \code{dsmm}. This class acts like a parent
#'    class for the classes \code{dsmm_fit_nonparametric,}
#'    \code{dsmm_fit_parametric, dsmm_parametric} and \code{dsmm_nonparametric}.
#' @param obj Arbitrary \code{R} object.
#' @seealso \link{is.dsmm_fit_nonparametric}, \link{is.dsmm_fit_nonparametric},
#'      \link{is.dsmm_parametric}, \link{is.dsmm_nonparametric}
#' @return TRUE or FALSE.
#' @export
is.dsmm <- function(obj) {
    # Check for missing object, `obj`.
    if (missing(obj)) {
        stop("\nPlease input the `obj` parameter of class `dsmm`. ",
             "This can be done through the functions\n `fit_dsmm()`, ",
             "`dsmm_parametric()` and `dsmm_nonparametric()`.")
    }
    # Check for class of object.
    if (!inherits(obj, 'dsmm')) {
        message("\n`obj` needs to be of class `dsmm`.",
                "\nThis can be done automatically through the functions:\n",
                "`fit_dsmm()`, `dsmm_parametric()` ",
                "and `dsmm_nonparametric()`.")
        return(FALSE)
    }
    check_attributes(obj)
}

#' @title Check if an object has a valid \code{dsmm_fit_nonparametric} class
#' @description Checks for the validity of the specified attributes and the
#'     inheritance of the S3 class \code{dsmm_fit_nonparametric}.
#'     This class inherits methods from the parent class \code{dsmm}.
#' @param obj Arbitrary \code{R} object.
#' @seealso \link{is.dsmm}, \link{is.dsmm_fit_parametric},
#'      \link{is.dsmm_nonparametric}, \link{is.dsmm_parametric}
#' @return TRUE or FALSE.
#' @export
is.dsmm_fit_nonparametric <- function(obj) {
    # Check for missing object, `obj`.
    if (missing(obj)) {
        stop("\nPlease input the `obj` parameter of class ",
             "`dsmm_fit_nonparametric`.",
             "\nThis can be done automatically through",
             " the function\n`fit_dsmm()`, when ",
             "`estimation` = 'nonparametric'.")
    }
    # Check for class of object.
    if (!inherits(obj, c('dsmm_fit_nonparametric'))) {
        message("\n`obj` does not have the class ",
                "('dsmm_fit_nonparametric', 'dsmm').",
                "\nThis can be done automatically through",
                " the function\n`fit_dsmm()`, when ",
                "`estimation` = 'nonparametric'.")
        return(FALSE)
    }
    check_attributes(obj)
}

#' @title Check if an object has a valid \code{dsmm_fit_parametric} class
#' @description Checks for the validity of the specified attributes and the
#'     inheritance of the S3 class \code{dsmm_fit_parametric}.
#'     This class inherits methods from the parent class \code{dsmm}.
#' @param obj Arbitrary \code{R} object.
#' @seealso \link{is.dsmm}, \link{is.dsmm_fit_nonparametric},
#'    \link{is.dsmm_parametric}, \link{is.dsmm_nonparametric}
#' @return TRUE or FALSE.
#' @export
is.dsmm_fit_parametric <- function(obj) {
    # Check for missing object, `obj`.
    if (missing(obj)) {
        stop("\nPlease input the `obj` parameter of class ",
             "`dsmm_fit_parametric`.\n",
             "This can be done automatically through",
             " the function\n`fit_dsmm()`, when ",
             "`estimation` = 'parametric'.")
    }
    # Check for class of object.
    if (!inherits(obj, c('dsmm_fit_parametric'))) {
        message("\n`obj` does not have the class ",
                "('dsmm_fit_parametric', 'dsmm').",
                "\nThis can be done automatically through",
                " the function\n`fit_dsmm()`, when ",
                "`estimation` = 'parametric'.")
        return(FALSE)
    }
    check_attributes(obj)
}

#' @title Check if an object has a valid \code{dsmm_nonparametric} class
#' @description Checks for the validity of the specified attributes and the
#'      inheritance of the S3 class \code{dsmm_nonparametric}.
#'      This class inherits methods from the parent class \code{dsmm}.
#' @seealso \link{is.dsmm}, \link{is.dsmm_fit_nonparametric},
#'     \link{is.dsmm_fit_parametric}, \link{is.dsmm_parametric}
#' @param obj Arbitrary \code{R} object.
#' @return TRUE or FALSE.
#' @export
is.dsmm_nonparametric <- function(obj) {
    # Check for missing object, `obj`.
    if (missing(obj)) {
        stop("\nPlease input the `obj` parameter of class `dsmm_nonparametric`.",
             "\nThis can be done automatically ",
             "through the function `nonparametric_dsmm()`.")
    }
    # Check for class of object.
    if (!inherits(obj, c('dsmm_nonparametric'))) {
        message("\n`obj` does not have the class ('dsmm_nonparametric', 'dsmm').",
                "\nThis can be done automatically ",
                "through the function `nonparametric_dsmm()`.")
        return(FALSE)
    }
    check_attributes(obj)
}

#' @title Check if an object has a valid \code{dsmm_parametric} class
#' @description Checks for the validity of the specified attributes and the
#'      inheritance of the S3 class \code{dsmm_parametric}.
#'      This class inherits methods from the parent class \code{dsmm}.
#' @param obj Arbitrary \code{R} object.
#' @seealso \link{is.dsmm}, \link{is.dsmm_fit_parametric},
#'      \link{is.dsmm_fit_nonparametric}, \link{is.dsmm_nonparametric}
#' @return TRUE or FALSE.
#' @export
is.dsmm_parametric <- function(obj) {
    if (missing(obj)) {
        stop("\nPlease input the `obj` parameter of class `dsmm_parametric`.",
             "\nThis can be done automatically ",
             "through the function `parametric_dsmm()`.")
    }
    # Check for class of object.
    if (!inherits(obj, c('dsmm_parametric'))) {
        message("\n`obj` does not have the class ('dsmm_parametric', 'dsmm').",
                "\nThis can be done automatically ",
                "through the function `parametric_dsmm()`.")
        return(FALSE)
    }
    check_attributes(obj)
}


# ______________________________________________________________________________
# Get the kernel q_(t/n) (u,v,l) that is necessary for the `simulate` function.
# ______________________________________________________________________________

#' @title Obtain the Drifting semi-Markov kernel
#'
#' @description
#' This is a generic method that computes and returns the drifting
#' semi-Markov kernel as a numerical array of dimensions
#' \eqn{s \times s \times k_{max} \times (n + 1)}.
#'
#' @param obj An object that inherits from the S3
#' classes \code{dsmm},
#' \code{dsmm_fit_parametric}, or
#' \code{dsmm_fit_nonparametric},
#' \code{dsmm_nonparametric} or \code{dsmm_parametric}.
#' @param t Optional, but recommended. Positive integer specifying
#' the instance \eqn{t} of the visited states.
#' @param u Optional. Can be either of the two options below:
#' \itemize{
#'   \item Character specifying the previous state \eqn{u}, e.g. \code{u = "a"}.
#'   \item Positive integer, specifying a state in the state space \eqn{E}.
#'   For example, if \eqn{E = \{a, c, g, t\}} and \code{u = 1}, it corresponds
#'   to the state \eqn{a}, if \code{u = 2}, it corresponds to the state \eqn{c}.
#' }
#' @param v Optional. Can be either of the two options below:
#' \itemize{
#'   \item Character specifying the next state \eqn{v}, e.g. \code{v = "c"}.
#'   \item Positive integer, specifying a state in the state space \eqn{E}.
#'   For example, if \eqn{E = \{a, c, g, t\}} and \code{v = 3}, it corresponds
#'   to the state \eqn{c}, if \code{v = 4}, it corresponds to the state \eqn{t}.
#' }
#' @param l Optional. Positive integer specifying the sojourn time \eqn{l}
#' that is spent in the previous state \eqn{u}.
#' @param klim Optional. Positive integer. Used only when \code{obj} inherits
#' from the S3 classes \code{dsmm_parametric} or \code{dsmm_fit_parametric}.
#' Specifies the time horizon used to approximate the \eqn{d + 1} sojourn time
#' distributions if \eqn{f} is drifting, or just \eqn{1} sojourn time
#' distribution if \eqn{f} is \emph{not drifting}.
#' Default value is 100.
#'
#' A larger value will result in a considerably larger
#' kernel, which has dimensions of \eqn{s \times s \times klim \times (n + 1)},
#' which will increase the memory requirements and will slow down considerably
#' the \code{simulate.dsmm()} method.
#' However, this will lead to better estimations through \code{fit_dsmm()}.
#' (\link{dsmm_parametric}, \link{fit_dsmm}, \link{simulate.dsmm})
#'
#'
#' @details
#' The drifting semi-Markov kernel is given as the probability that,
#' given at the instance \eqn{t} the previous state
#' is \eqn{u}, the next state state \eqn{v} will be reached
#' with a sojourn time of \eqn{l}:
#' \deqn{q_{\frac{t}{n}}(u,v,l) = P(J_{t}=v,X_{t}=l|J_{t-1}=u),}
#' where \eqn{n} is the model size, defined as the length of the embedded
#' Markov chain \eqn{(J_{t})_{t\in \{0,\dots,n\}}} minus the last state,
#' \eqn{J_t} is the visited state at the instant \eqn{t} and
#' \eqn{X_{t} = S_{t}-S_{t-1}} is the sojourn time of the state \eqn{J_{t-1}}.
#' Specifically, it is given as the sum of a linear combination:
#' \deqn{q_{\frac{t}{n}}(u,v,l)=
#'      \sum_{i = 0}^{d}A_{i}(t)\ q_{\frac{i}{d}}(u,v,l),}
#' where \eqn{A_i, i = 0, \dots, d} are \eqn{d + 1} polynomials with degree
#' \eqn{d} that satisfy certain conditions (see \link{dsmmR}) and
#' \eqn{q_{\frac{i}{d}}(u,v,l), i = 0, \dots, d}
#' are \eqn{d + 1} semi-Markov kernels.
#' Three possible model specifications are described below.
#' We will use the exponentials \eqn{(1), (2), (3)} to distinguish between
#' the drifting semi-Markov kernel \eqn{q_\frac{t}{n}} and the
#' semi-Markov kernels \eqn{q_\frac{i}{d}} used in
#' Model 1, Model 2 and Model 3.
#'
#' \strong{\emph{Model 1}}
#'
#' In this case, both \eqn{p} and \eqn{f} are "drifting" between \eqn{d + 1}
#' fixed points of the model, hence the "drifting" in drifting semi-Markov
#' models. Therefore, the semi-Markov kernels \eqn{q_{\frac{i}{d}}^{\ (1)}} are
#' equal to:
#'
#' \deqn{q_{\frac{i}{d}}^{\ (1)}(u,v,l) =
#'      {p_{\frac{i}{d}}(u,v)}{f_{\frac{i}{d}}(u,v,l)},}
#'
#' where for \eqn{i = 0, \dots, d} we have \eqn{d + 1} Markov Transition
#' matrices \eqn{p_{\frac{i}{d}}(u,v)}, and \eqn{d + 1} sojourn time
#' distributions \eqn{f_{\frac{i}{d}}(u,v,l)}, where \eqn{d} is the
#' polynomial degree.
#'
#' Thus, the drifting semi-Markov kernel will be equal to:
#'
#' \deqn{q_{\frac{t}{n}}^{\ (1)}(u,v,l) =
#' \sum_{i = 0}^{d} A_i(t)\ q_{\frac{i}{d}}^{\ (1)}(u,v,l) =
#' \sum_{i = 0}^{d} A_i(t)\ p_{\frac{i}{d}}(u,v)f_{\frac{i}{d}}(u,v,l)
#' }
#'
#'
#' \strong{\emph{Model 2}}
#'
#' In this case, \eqn{p} is drifting and \eqn{f} \strong{is not drifting}.
#' Therefore, the semi-Markov kernels \eqn{q_{\frac{i}{d}}^{\ (2)}} are
#' equal to:
#' \deqn{q_{\frac{i}{d}}^{\ (2)}(u,v,l)={p_{\frac{i}{d}}(u,v)}{f(u,v,l)}.}
#'
#' Thus, the drifting semi-Markov kernel will be equal to:
#' \deqn{q_{\frac{t}{n}}^{\ (2)}(u,v,l) =
#' \sum_{i = 0}^{d} A_i(t)\ q_{\frac{i}{d}}^{\ (2)}(u,v,l) =
#' \sum_{i = 0}^{d} A_i(t)\ p_{\frac{i}{d}}(u,v)f(u,v,l)
#' }
#'
#'
#' \strong{\emph{Model 3}}
#'
#' In this case, \eqn{f} is drifting and \eqn{p} \strong{is not drifting}.
#'
#' Therefore, the semi-Markov kernels \eqn{q_{\frac{i}{d}}^{\ (3)}}
#' are now described as:
#' \deqn{q_{\frac{i}{d}}^{\ (3)}(u,v,l)={p(u,v)}{f_{\frac{i}{d}}(u,v,l)}.}
#'
#' Thus, the drifting semi-Markov kernel will be equal to:
#' \deqn{q_{\frac{t}{n}}^{\ (3)}(u,v,l) =
#' \sum_{i = 0}^{d} A_i(t)\ q_{\frac{i}{d}}^{\ (3)}(u,v,l) =
#' \sum_{i = 0}^{d} A_i(t)\ p(u,v)f_{\frac{i}{d}}(u,v,l)
#' }
#'
#' @return An array with dimensions of
#' \eqn{s \times s \times k_{max} \times (n + 1)}, giving the
#' value of the drifting semi-Markov kernel \eqn{q_{\frac{t}{n}}(u,v,l)} for
#' the corresponding \eqn{(u,v,l,t)}. If any of \eqn{u,v,l} or \eqn{t} were
#' specified, their dimension in the array becomes 1.
#'
#' @export
#'
#' @seealso
#' For the objects required to calculate this kernel:
#' \link{fit_dsmm}, \link{parametric_dsmm}, \link{nonparametric_dsmm}.
#'
#' For sequence simulation through this kernel: \link{simulate.dsmm}.
#'
#' For the theoretical background of drifting semi-Markov models: \link{dsmmR}.
#'
#' @examples
#' # Setup.
#' states <- c("Rouen", "Bucharest", "Samos", "Aigio", "Marseille")
#' seq <- create_sequence(states, probs = c(0.3, 0.1, 0.1, 0.3, 0.2))
#' obj_model_2 <- fit_dsmm(
#'     sequence = seq,
#'     states = states,
#'     degree = 3,
#'     f_is_drifting = FALSE,
#'     p_is_drifting = TRUE
#' )
#'
#' # Get the kernel.
#' kernel_model_2 <- get_kernel(obj_model_2)
#' cat(paste0("If no further arguments are made, kernel has dimensions ",
#'            "for all u, v, l, t:\n",
#'            "(s, s, k_max, n + 1) = (",
#'            paste(dim(kernel_model_2), collapse = ", "), ")"))
#'
#' # Specifying `t`.
#' kernel_model_2_t <- get_kernel(obj_model_2, t = 100)
#' # kernel_model_2_t[ , , , t = 100]
#' cat(paste0("If we specify t, the kernel has dimensions for ",
#'            "all the remaining u, v, l:\n(s, s, k_max) = (",
#'            paste(dim(kernel_model_2_t), collapse = ", "), ")"))
#'
#' # Specifying `t` and `u`.
#' kernel_model_2_tu <- get_kernel(obj_model_2, t = 2, u = "Aigio")
#' # kernel_model_2_tu["Aigio", , , t = 2]
#' cat(paste0("If we specify t and u, the kernel has dimensions for ",
#'            "all the remaining v, l:\n(s, k_max) = (",
#'            paste(dim(kernel_model_2_tu), collapse = ", "), ")"))
#'
#' # Specifying `t`, `u` and `v`.
#' kernel_model_2_tuv <- get_kernel(obj_model_2, t = 3,
#'                                  u = "Rouen", v = "Bucharest")
#' # kernel_model_2_tuv["Rouen", "Bucharest", , t = 3]
#' cat(paste0("If we specify t, u and v, the kernel has dimensions ",
#'            "for all l:\n(k_max) = (",
#'            paste(length(kernel_model_2_tuv), collapse = ", "), ")"))
#'
#' # It is possible to ask for any valid combination of `u`, `v`, `l` and `t`.
get_kernel <- function(obj, t, u, v, l, klim = 100) {
    # Check for missing values & the validity of the attributes given.
    if (missing(obj)) {
        stop("\nPlease provide the object `obj`.")
    }
    stopifnot(check_attributes(obj))
    n <- obj$model_size
    N <- n + 1
    k_max <- obj$k_max
    states <- obj$states
    s <- length(states)
    if (!missing(t)) {
        if (!is_integer(t)) {
            stop("\nAttribute `t` should be a positive integer, specifying",
                 " the instance of the visited state of your choice.\n",
                 "Currently, it is equal to: ", t)
        }
        if (t > N) {
            stop("\n`t` cannot be larger than n + 1, where n = ", n)
        }
    }
    if (!missing(u)) {
        if (is.character(u)) {
            stopifnot(valid_state(u, states))
        } else if (is_integer(u) && u > s) {
            stop("\nThe previous state `u` is specified as the numbered ",
                 "state `", u, "` in the state space of total length s = ",
                 s, ".")
        }
        u <- states[which(states == u)]
    }
    if (!missing(v)) {
        if (is.character(v)) {
            stopifnot(valid_state(v, states))
        } else if (is_integer(v) && v > s) {
            stop("\nThe previous state `v` is specified as the numbered",
                 " state `", v,
                 "` in the state space of total length s = ", s, ".")
        }
        v <- states[which(states == v)]
    }
    if (!missing(l)){
        if (!is_integer(l)) {
            stop("\nAttribute `l` should be a positive integer specifying the",
                 " sojourn time.\n",
                 "Currently, it is equal to: ", l)
        }
        if (l > k_max) {
            stop("\nThe maximum sojourn time specified, `", l,
                 "` cannot be more than ", k_max, ".")
        }
    }
    if (!is_integer(klim)) {
        stop("\nAttribute `klim` should be a positive integer.\n",
             "Currently, it is equal to: ", klim)
    } else if (klim > 200) {
        warning("\nAttribute `klim` = ", klim, " will result in",
                " considerably expensive memory requirements for",
                " a larger model size.")
    }
    obj$k_max <- klim
    # Possible defined parameters are passed down to `UseMethod`.
    UseMethod(generic = 'get_kernel', object = obj)
}

#' @export
get_kernel.dsmm <- function(obj, t, u, v, l, ...) {
    # =========================================================================.
    # '''
    #    This function is valid for both `dsmm_fit_nonparametric`
    #    AND `dsmm_nonparametric`.
    # '''
    # =========================================================================.
    # Get the correct distributions `dist`.
    D <- obj$degree + 1L
    dist <- obj$dist
    n <- obj$model_size
    k_max <- obj$k_max
    states <- obj$states
    s <- length(states)
    Ai <- obj$A_i # this should be a matrix with dim(A_i) == c(d+1, n+1)
    Model <- obj$Model
    pdist <- dist[[1]]
    fdist <- dist[[2]]
    # Get J_i with regards to Model.
    if (Model == 'Model_1') {
        Ji <- fdist * c(apply(pdist, c(3),
                              function(M_uv) rep(M_uv, times = k_max)))
    } else if (Model == "Model_2") {
        # dim(fdist) == (s, s, k_max) --> (s, s, k_max, d + 1).
        f_vector <- rep(fdist, D)
        # dim(pdist) == (s, s, d+1) --> (s, s, k_max, d + 1)
        p_vector <- apply(pdist, MARGIN = c(3),
                          FUN = function(x) rep(x, k_max))
        Ji <- f_vector * p_vector
    } else if (Model == "Model_3") {
        # dim(pdist) == (s, s) --> (s, s, k_max, d+1)
        p_vector <- rep(pdist, k_max * D)
        Ji <- fdist * p_vector
    }
    # Get kernel q_(t/n) (u,v,l).
    dim(Ji) <- c(s*s*k_max, D)
    # dimnames(Ji) <-list(
    #   as.list(paste0(rep(E, s), sort(rep(E, s)),
    #    sort(rep(1:k_max, s * s)))),
    #   as.list(names_i_d(degree, "q"))
    # )
    kernel <- get_valid_kernel(Ji, Ai, s, n, k_max, states)
    if (!is_prob(kernel)) stop("\nThe final kernel is not stochastic.")
    kernel[u, v, l, t]
}

#' @export
get_kernel.dsmm_parametric <- function(obj, t, u, v, l, klim = 100) {
    # '''
    #    This function returns the kernel q_(t/n) with variables (t, u, v, l).
    # '''
    # Does not need missing() case, since there is only one argument
    # and it needs to be passed down from the generic function.
    stopifnot(check_attributes(obj))
    # Get the correct distributions `dist`.
    pdist <- obj$dist[[1]] # [u,v]
    fdist <- obj$dist[[2]] # [u,v,]
    fpar <- obj$dist[[3]]  # [u,v,,]
    #  A large klim will create a large kernel
    #  which will slow down the computation
    #  for the method simulate.dsmm()` considerably.
    f_vector <- get_fdist_parametric(fdist = fdist,
                                     params = fpar,
                                     klim = klim)
    # f_vector[which(f_vector < 1e-10)] <- 0 # this is not necessary...
    p_vector <- rep(pdist, each = klim)
    Ji <- f_vector * p_vector
    degree <- obj$degree
    D <- degree + 1
    dim(Ji) <- c(klim, s <- obj$s, s, D)
    Ji <- aperm(Ji, c(2, 3, 1, 4))
    dim(Ji) <- c(s * s * klim, D)
    Ai <- get_A_i(degree, n <- obj$model_size)
    kernel <- Ji %*% Ai
    dim(kernel) <- c(s, s, klim, (N <- n + 1))
    # dimnames(kernel) <- list(
    #     as.list(states <- obj$states),
    #     as.list(states),
    #     as.list(paste0('l = ', 1:klim)),
    #     as.list(paste0('t = ', 0:n))
    # )
    kernel <- apply(kernel, c(1,4), function(vl_values) {
        if (any(vl_values < 0)) {
            return(get_f(vl_values))
        } else {
            return(vl_values)
        }
    })
    dim(kernel) <- c(s, klim, s, N)
    kernel <- aperm(kernel, c(3,1,2,4))
    dimnames(kernel) <- list(
        as.list(states <- obj$states),
        as.list(states),
        as.list(paste0('l = ', 1:klim)),
        as.list(paste0('t = ', 0:n))
    )
    kernel[u, v, l, t]
}

#' @export
get_kernel.dsmm_fit_parametric <- function(obj, t, u, v, l, klim = 100) {
    # '''
    #    This function returns the kernel q_(t/n) with variables (t, u, v, l).
    #    The same technique is used as in `get_kernel.dsmm_parametric()`
    # '''
    # Does not need missing() case, since there is only one argument
    # and it needs to be passed down from the generic function.
    get_kernel.dsmm_parametric(obj, t, u, v, l, klim)
}


# ______________________________________________________________________________
# Printing methods for each of the child classes.
# ______________________________________________________________________________
#' @export
print.dsmm <- function(x, ...) {
    # '''
    #    This function acts on both classes of `dsmm_fit_parametric` and
    #    `dsmm_fit_nonparametric`.
    #    It was made for the purpose of NOT printing certain
    #    parameters. For example, the polynomials `A_i` are too long
    #    and would obscure the printing of the object.
    # '''
    check_attributes(x)
    nm <- names(x)
    for (i in seq_along(nm)) {
        if (!nm[i] %in% c('A_i', 'dist', 'J_i', 'soj_times', 'seq')) {
            cat(paste0("\n$", nm[i], "\n"))
            if (nm[i] %in% c('k_max', 'model_size', 's',
                             'degree', 'states',
                             'f_is_drifting', 'p_is_drifting',
                             'numerical_est', 'Model', 'estimation')) {
                cat(x[[i]], "\n")
            } else {
                print(x[[i]])
            }
        } else if (nm[i] %in% c('soj_times', 'seq')) {
            cat(paste0("\n$", nm[i], "\n"))
            print(head((xx <- x[[i]]), n = 100L))
            cat(" ... [ output truncated at 100 values -- ommited ",
                length(xx) - 100, " entries ]\n")
        } else if (nm[i] == 'dist') {
            nd <- names(x$dist)
            for (j in seq_along(nd)) {
                cat(paste0("\n\n\n$dist$", nd[j], "\n\n"))
                print(x$dist[[j]])
            }
        }
    }
    cat("\n\nClass:", paste0(class(x), collapse = ", "))
}

#' @export
print.dsmm_nonparametric <- function(x, ...) {
    # '''
    #    This function was made for the purpose of NOT printing
    #    certain parameters. For example, the polynomials `A_i`
    #    are too long and would obscure the printing of the object.
    # '''
    check_attributes(x)
    nm <- names(x)
    for (i in seq_along(nm)) {
        if (!nm[i] %in% c('A_i', 'dist')) {
            cat(paste0("\n$", nm[i], "\n"))
            if (nm[i] %in% c('k_max', 'model_size', 's', 'degree', 'states',
                             'f_is_drifting', 'p_is_drifting', 'Model')) {
                cat(x[[i]], "\n")
            } else {
                print(x[[i]])
            }
        } else if (nm[i] == 'dist') {
            nd <- names(x$dist)
            for (j in seq_along(nd)) {
                cat(paste0("\n\n\n$dist$", nd[j], "\n\n"))
                print(x$dist[[j]])
            }
        }
    }
    cat("\n\nClass:", paste0(class(x), collapse = ", "))
}

#' @export
print.dsmm_parametric <- function(x, ...) {
    # '''
    #    This function was made for the purpose of NOT printing certain
    #    parameters. For example, the polynomials `A_i` are too long
    #    and would obscure the printing of the object.
    # '''
    check_attributes(x)
    nm <- names(x)
    for (i in seq_along(nm)) {
        if (!nm[i] %in% c('A_i', 'dist')) {
            cat(paste0("\n$", nm[i], "\n"))
            if (nm[i] %in% c('model_size', 's', 'degree', 'states',
                             'f_is_drifting', 'p_is_drifting', 'Model')) {
                cat(x[[i]], "\n")
            } else {
                print(x[[i]])
            }
        } else if (nm[i] == 'dist') {
            nd <- names(x$dist)
            for (j in seq_along(nd)) {
                cat(paste0("\n\n\n$dist$", nd[j], "\n\n"))
                print(x$dist[[j]])
            }
        }
        # else if (nm[i] == 'f_distribution_parameters') {
        #     print(x$f_distribution_parameters)
        # }
    }
    cat("\n\nClass:", paste0(class(x), collapse = ", "))
}

# ______________________________________________________________________________
# Simulate a sequence from any `dsmm` object.
# ______________________________________________________________________________
#' @title Simulate a sequence under a drifting semi-Markov kernel.
#'
#' @description Generic function that simulates a number of states \code{nsim}
#' under the rule of a drifting semi-Markov kernel, which is retrieved from the
#' object \code{obj}, which in turn inherits from the S3 class \code{dsmm}.
#'
#' @param object An object of S3 class \code{dsmm},
#' \code{dsmm_fit_nonparametric}, \code{dsmm_nonparametric},
#' \code{dsmm_fit_parametric} or \code{dsmm_parametric}.
#'
#' @param nsim Optional. An integer specifying the number of simulations to be made
#' from the drifting semi-Markov kernel. The maximum value of \code{nsim} is the
#' model size which is specified in \code{obj}, which is also the default value.
#' We define a special case for \code{nsim = 0}, where only the initial distribution
#' is considered and only the simulation of its sojourn time will be made, without
#' the next state.
#'
#' @param seq_length Optional. A positive integer that will ensure the simulated
#' sequence will not have a \emph{total length} greater than \code{seq_length}
#' (however, it is possible for the total length to be \emph{less} than
#' \code{seq_length}).
#'
#' @param seed Optional. An integer specifying the initialization of the random
#' number generator.
#'
#' @param klim Optional. Positive integer. Passed down to \code{get_kernel}
#' for the parametric object, with class \code{dsmm_parametric}.
#' Default value is \eqn{100}.
#'
#' @param ... Optional. Attributes passed down from the \code{simulate} method.
#'
#' @seealso
#' About random number generation in R: \code{\link[base:RNG]{RNG}}.
#'
#' Fitting a model through a sequence from this function: \link{fit_dsmm}.
#'
#' For the theoretical background of drifting semi-Markov models: \link{dsmmR}.
#'
#' @return A character vector based on \code{nsim} simulations, with a
#' maximum length of \code{seq_length}.
#' @export
#' @examples
#' # Setup.
#' seq <- create_sequence("DNA", len = 1000)
#' states <- sort(unique(seq))
#' d <- 1
#' obj_model_3 <- fit_dsmm(sequence = seq,
#'                         states = states,
#'                         degree = d,
#'                         f_is_drifting = TRUE,
#'                         p_is_drifting = FALSE)
#'
#' # Using the method `simulate.dsmm()`.
#' simulated_seq <- simulate(obj_model_3, seed = 1)
#' short_sim <- simulate(obj = obj_model_3, nsim = 10, seed = 1)
#' cut_sim <- simulate(obj = obj_model_3, seq_length = 10, seed = 1)
#' str(simulated_seq)
#' str(short_sim)
#' str(cut_sim)
simulate.dsmm <- function(object, nsim = NULL, seed = NULL,
                          seq_length = NULL, klim = 100, ...) {
    # Parameters Setup.
    if (missing(object)) {
        stop("\nPlease provide an objectect of class `dsmm`.")
    } else if (!inherits(object, c('dsmm', 'dsmm_fit', 'dsmm_nonparemetric',
                                   'dsmm_parametric'))) {
        stop("\nPlease provide an object of class `dsmm` to use for the",
             " function `simulate()`.",
             "\nThe object can be created through the functions ",
             "`parametric_dsmm()`, `nonparametric_dsmm()` and `fit_dsmm()`.")
    }
    # Check if nsim and sequence length are given at the same time.
    if (!is.null(nsim) && !is.null(seq_length)) {
        stop("\nPlease specify only one of `nsim` or `seq_length` for ",
             "the simulation.")
    }
    # Check `nsim`.
    if (is.null(nsim)) {
        nsim <- object$model_size
    } else if (!(is_integer(nsim) || nsim == 0)) {
        stop("\nThe number of simulations `nsim` ",
             "needs to be a positive integer or 0.")
    }
    # Check `seq_length`.
    if (!is.null(seq_length) && !is_integer(seq_length)) {
        stop("\nThe final length of the sequence `seq_length` ",
             "needs to be a positive integer.")
    } else if (is_integer(seq_length)) {
        nsim <- seq_length
    }
    # Set the RNG seed.
    if (!is.null(seed)) {
        set.seed(seed)
    } # Otherwise, set.seed is automatic through the user's `Sys.time()`.
    if (!is_integer(klim)) {
        stop("\nAttribute `klim` should be a positive integer.")
    }
    stopifnot(check_attributes(object))
    # In order to get the first letter of the sequence.
    alpha <- object$initial_dist
    states <- object$states
    initial_state <- sample(states, prob = alpha, size = 1)
    s <- length(states)
    n <- object$model_size
    if (nsim > n) {
        stop("\nThe number of simulations `nsim` = ", nsim,
             " cannot be larger than the model size, n = ", n)
    }
    kernel <- get_kernel(object, klim = klim)
    k_max <- dim(kernel)[3] # We get `k_max` even for the parametric case.
    # Get the rest of the sequence.
    vl_names <- paste(states, sort(rep(1:k_max, s)))
    vl_vector <- c()
    u <- initial_state
    if (nsim == 0) {
        vl <- get_vl(t = 1, u = u, vl_names = vl_names,
                     kernel = kernel, states = states)
        sim_seq <- as.vector(rep(vl[1], vl[2]))
        return(sim_seq)
    }
    for (time in 1:nsim) { # max(nsim) == n!
        u <- states[which(states == u)]
        vl_vector <- c(vl_vector,
                       vl <- get_vl(t = time, u = u, vl_names = vl_names,
                                    kernel = kernel, states = states))
        u <- vl[1]
    }
    l_vl <- length(vl_vector)
    seq <- c(initial_state, vl_vector[seq(1, l_vl, 2)])
    X <- as.numeric(c(vl_vector[seq(2, l_vl, by = 2)], 1))
    sim_seq <- as.vector(unlist(sapply(seq_along(seq),
                                       function(i) rep(seq[i], X[i]))))
    if (is.null(seq_length)) {
        return(sim_seq)
    } else {
        return(sim_seq[1:seq_length])
    }
}

