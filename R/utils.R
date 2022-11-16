# '''
#    1. This file contains functions that are used for supporting the other
#       functions in the package.
#    2. It is worth noting that none of the following functions are available
#       to the user directly, apart from the last one, `create_sequence()`.
#    3. For the user to gain access to these functions, he has to use
#       dsmmR:::foo()
# '''

# ==============================================================================
# Check for equality.
# ==============================================================================
all_equal <- function(obj1, obj2) isTRUE(all.equal(obj1, obj2))

all_equal_numeric <- function(obj1, obj2, tol = sqrt(.Machine$double.eps))
    all(abs(obj1 - obj2) < tol)

all_equal_numeric_vector <- function(obj1, obj2,
                                     tol = sqrt(.Machine$double.eps))
    abs(obj1 - obj2) < tol

# ==============================================================================
# Check for types and classes of objects.
# ==============================================================================
is_vector <- function(obj)
    is.vector(obj) && !is.list(obj) && !is.expression(obj)

is_logical <- function(obj) {
    # '''
    #    Returns TRUE or FALSE if `obj` is a TRUE or FALSE logical vector
    #    with ___length 1____.
    # '''
    is_vector(obj) && (length(obj) == 1) && (isTRUE(obj) || isFALSE(obj))
}


is_number <- function(obj) {
    # '''
    #   Checks if `obj` is a vector of numbers with length = 1.
    # '''
    is_vector(obj) && (length(obj) == 1) && is.numeric(obj) &&
        !is.logical(obj) && is.finite(obj) && !is.nan(obj)
}


is_double <- function(obj) {
    # '''
    #   Checks if `obj` is a double vector with length = 1.
    # '''
    is_number(obj) && is.double(obj)
}


is_double_vector <- function(obj) # Used in `valid_initial_dist()`.
    is_vector(obj) && is.double(obj) && !is_number(obj)


is_double_matrix <- function(obj) {
    # '''
    #     This is used in `is_fdist_parametric` for checking `params` array.
    #     The `params` array should contain NA values.
    # '''
    is.matrix(obj) && !any(is.nan(obj)) &&
        !any(is.infinite(obj)) && all(is.double(obj))
}

is_double_array <- function(obj) {
    # '''
    #     This is used in `is_fdist_parametric` for checking `params` array.
    #     The `params` array should contain NA values.
    # '''
    is.array(obj) && !is.matrix(obj) && !any(is.nan(obj)) &&
        !any(is.infinite(obj)) && all(is.double(obj))
}

is_integer <- function(obj) {
    # '''
    #   Checks if `obj` is a number that can be represented as an integer.
    #   E.g. for the double 10 not equal to the integer 10L, it returns TRUE.
    # '''
    is_number(obj) && obj > 0 && obj - as.integer(obj) == 0
}

is_integer_vector <- function(obj) # used in `valid_soj_times()`.
    is_vector(obj) && all(sapply(obj, is_integer))

is_prob <- function(prob) {
    # '''
    #    This function checks whether 'obj' is `sufficiently close`
    #     to the domain [0, 1]
    #    (`sufficiently` is with accordance to `sqrt(.Machine$double.eps)`.
    #    E.g. 0.00001, 1.00000001, -0.00000001, 0 and 1 all return TRUE, but
    #         -0.0000001 returns FALSE (not `sufficiently` far from 0).
    # '''
    !( # The opposite of the below is a valid way to define a probability.
        ( (!all_equal_numeric(min <- min(prob), 0)) && min < 0) ||
            # `prob` not close to 0 & less than 0 OR
            ( (!all_equal_numeric(max <- max(prob), 1)) && max > 1)
        # `prob` not close to 1 & more than 1
    )
}

is_prob_vector <- function(prob) {
    # '''
    #    This function checks whether 'obj' is `sufficiently close`
    #     to the domain [0, 1]
    #    (`sufficiently` is with accordance to `sqrt(.Machine$double.eps)`.
    #    E.g. 0.00001, 1.00000001, -0.00000001, 0 and 1 all return TRUE, but
    #         -0.0000001 returns FALSE (not `sufficiently` far from 0).
    # '''
    !( # The opposite of the below is a valid way to define a probability.
        ( (!all_equal_numeric_vector(prob, 0)) & prob < 0) |
            # `prob` not close to 0 & less than 0 OR
            ( (!all_equal_numeric_vector(prob, 1)) & prob > 1)
        # `prob` not close to 1 & more than 1
    )
}


# ==============================================================================
# Check for the valid use of the different arguments that the user can input
# ==============================================================================
valid_sequence <- function(sequence, s) {
    # '''
    #    Used in `fit_dsmm`.
    # '''
    if (!is.character(sequence)){
        stop("\nThe `sequence` argument should be a character vector.")
    } else if (length(unique(sequence)) < s) {
        # check for (very) short sequence.
        stop("\nThe `sequence` argument should have more than s = ",
             s, " values.")
    }
    TRUE
}

valid_seq <- function(seq) {
    # '''
    #    This function is for checking the correct use of the 'trimmed'
    #    sequence that has all sojourn times equal to 1. This is not
    #    about the correct use of the original sequence given e.g. for
    #    `fit_dsmm`.
    #     Intended for use in the function `is.dsmm`.
    # '''
    if (!is.character(seq)){
        stop("\n`seq` argument should be a character vector.")
    } else if (!all_equal(seq, rle(seq)$values)) {
        # check for seq not having any recurrent states.
        stop("\n`seq` argument should have all sojourn times equal to 1.")
    }
    TRUE
}

valid_soj_times <- function(soj_times, length_seq) {
    # '''
    #    Intended for use in `is.dsmm`.
    # '''
    if (!is_integer_vector(soj_times)) {
        stop("\nAttribute `soj_times` should be an integer vector.")
    } else if (!all_equal_numeric(length(soj_times), length_seq)) {
        stop("\nAttribute `soj_times` should have length equal to that of ",
             "attribute `seq`.")
    }
    TRUE
}

valid_k_max <- function(k_max, soj_times) {
    # '''
    #   This functions checks for correct behavior of the `k_max` parameter.
    #   Intended for use in `fit.dsmm`, `is.dsmm` & `nonparametric_dsmm`.
    # '''
    if (is_integer(k_max)) {
        if (k_max == max(soj_times)) {
            return(TRUE)
        } else {
            stop("\nThe maximum sojourn time attribute `k_max` should be",
                 " an integer equal to the maximum of `soj_times`.")
        }
    }
    stop("\nThe maximum sojourn time `k_max` should be an integer.")
}

valid_model_size <- function(model_size, length_seq) {
    # '''
    #    Intended for use in `is.dsmm`.
    #    `length_seq` is with regards to the sequence without the
    #    sojourn times, e.g. from the `rle()`function.
    # '''
    if (is_integer(model_size)) {
        if (!all_equal_numeric(model_size, length_seq - 1)) {
            stop("\n`model_size` should be equal to the length of `seq`.")
        }
    }
    TRUE
}

valid_states <- function(states) {
    # '''
    #   This functions checks for correct behavior of the `states` parameter.
    # '''
    if (!is.character(states) ||
        !is_vector(states) ||
        !all_equal(states, unique(states))) {
        stop("\nThe `states` of the State Space should be a character vector",
             " of unique values.")
    } else if (length(states) == 1L) {
        # Check for a single state.
        stop("\nState Space `states` should have more than one states.")
    } else if (!all_equal(states, unique(states))) {
        # Check for uniqueness of states.
        stop("\nThe state space given in `states` contains duplicate values.")
    } else if(min(prob <- (sapply(strsplit(
        x = gsub(" ", "", states), split = ""), length))) == 0) {
        # Check for a state equal to an empty string "" or `character(0)`.
        stop("\nState Space `states` includes a state without a given name,",
             "at position ", which(prob == 0))
    }
    TRUE
}

valid_state <- function(state, states) {
    # '''
    #    Intended for use in `get_kernel`.
    # '''
    if (!is.character(state)) {
        stop("\nAttribute `state` should be a character.")
    }
    stopifnot(valid_states(states))
    if (!state %in% states) {
        stop("\nState `", state, "` is not included in the State",
             " Space that is given,\nE = (", paste(states, collapse = ', '),
             ").")
    }
    TRUE
}

valid_length_states <- function(s, states) {
    # '''
    #    Intended for use in `is.dsmm`.
    # '''
    if (!(is_integer(s) && s == (lstates <- length(states)))) {
        stop("\nAttribute `s` = ", s,
             " should be an integer equal to the length of `states` = ",
             lstates)
    }
    TRUE
}

valid_degree <- function(degree, model_size) {
    # '''
    #   This functions checks for correct behavior of the `degree` parameter.
    # '''
    if (!is_integer(degree)) {
        stop("\nThe the polynomial `degree` should be a positive integer.")
    }
    if (!missing(model_size)) {
        if (degree > model_size) {
            stop("\nThe polynomial `degree` = ", degree,
                 " should not be larger than the model size = ",
                 model_size, ".")
        }
    }
    TRUE
}

valid_initial_dist <- function(initial_dist, s) {
    # '''
    #   This functions checks for correct behavior of the `initial_dist`
    #   parameter.
    # '''
    if (!is_double_vector(initial_dist)) {
        stop("\nThe `initial_dist` should be a vector of 'double' ",
             "numeric values.")
    } else if ((n_init <- length(initial_dist)) != s) {
        stop("\nThe `initial_dist` has ", n_init, " values and it ",
             "should be equal",
             " to the number of `states` given, ", s, ".")
    } else if (!is_prob(initial_dist)) {
        stop("\nThe `initial_dist` should contain numeric values ",
             "between 0 and 1.")
    } else if (!sum(initial_dist) == 1) {
        stop("\nThe sum of the `initial_dist` should be exactly ",
             "equal to 1.")
    }
    TRUE
}

valid_initial_dist_char <- function(obj) {
    # '''
    #     Returns TRUE or FALSE if `obj` is a single character in the list
    #     ['unif', 'freq'].
    #     Used in `fit_dsmm()`.
    # '''
    if (!(is.character(obj) && length(obj) == 1 &&
          (obj == "unif" || obj == "freq"))) {
        stop('\n`initial_dist` should be either "unif" or "freq".')
    }
    TRUE
}

valid_model <- function(p_is_drifting, f_is_drifting) {
    # '''
    #     Intended for use in `is.dsmm_fit` & `fit_dsmm`, as well as
    #    `dsmm_parametric`, `dsmm_nonparametric`,
    #    `is.dsmm_parametric` and `is.dsmm_nonparametric`
    # '''
    if (!p_is_drifting && !f_is_drifting) {
        # Neither p or f are drifting.
        stop("\nAt least p OR f should be drifting to use a Drifting",
             "semi-Markov Model",
             " for the estimation. Otherwise, a semi-Markov Model",
             "should be used.")
    }
    TRUE
}

valid_estimation <- function(estimation, fpar, s, f_is_drifting, degree,
                             states) {
    # '''
    #    This function is used in `fit_dsmm` with the purpose of
    #    checking whether the attributes `estimation` and `fpar`
    #    correspond to either "nonparametric" and NULL
    #    or to "nonparametric" and `f_parametric`.
    # '''
    if (!is.character(estimation) ||
        (estimation != "nonparametric" && estimation != "parametric")) {
        stop('`estimation` attribute should be either "nonparametric" for the',
             ' non-parametric estimation or "parametric" for the parametric',
             ' estimation.')
    }
    if (estimation == "nonparametric") {
        if (!is.null(fpar)) {
            stop('When `estimation = "nonparametric"`, attribute `fpar` ',
                 'should be NULL')
        }
    } else if (estimation == "parametric") {
        stopifnot(valid_fdist_parametric(fdist = fpar, params = NULL, s = s,
                                         f_is_drifting = f_is_drifting,
                                         degree = degree, states = states))
    }
    TRUE
}

valid_p_dist <- function(p_dist, s, degree, p_is_drifting, states) {
    # '''
    #    Intended for use in `dsmm_nonparametric` & `dsmm_parametric`.
    # '''
    D <- degree + 1L
    if (p_is_drifting) {
        # Case when p is drifting, a.k.a. `p_is_drifting` = TRUE.
        if (!is_double_array(p_dist)) {
            stop("\nSince `p_is_drifting` is TRUE, `p_dist` should be a ",
                 "numeric array with dimensions:\n(s, s, degree + 1).")
        } else if (!is_prob(p_dist)) {
            non_prob <- which(!is_prob_vector(p_dist))
            stop("\nArray `p_drift` contains values that are not ",
                 "probabilities:\n",
                 paste0(p_dist[non_prob], collapse = ", "))
        } else if (!all_equal((dimension <- dim(p_dist)),
                              c(s, s, D))) {
            stop("\nSince `p_is_drifting` = TRUE, `p_dist` should be an",
                 "array with dimensions:\n(s, s, degree + 1) = (",
                 paste(c(s, s, D), collapse = ', '),
                 ") when it has dimensions of: (",
                 paste(dimension, collapse = ', '), ").")
        } else if (!all_equal_numeric(
            sum(p_dist * array(diag(s), dim = c(s, s, D))), 0)) {
            stop("\nThe diagonal values of `p_dist` should all be",
                 " equal to 0.")
        } else if (!all_equal_numeric(
            (p_drift_rowsum <-
             c(sapply(1:D, function(u) rowSums(p_dist[ , , u])))),
            1)) {
            # Case where, for every d + 1 matrix, the sums over v for every u
            #  are not equal to 1.
            logical_vector <- sapply(p_drift_rowsum,
                                     function(u) !all_equal(u, 1))
            possible_cases <- paste0("u = ", rep(states, s - 1), ", ",
                                     rep(names_i_d(degree, kernel_name = "p"),
                                         each = s))
            stop("\nFor the transition probability matrices",
                 " contained in `p_dist`, the probabilities of transistioning",
                 " from previous state `u` over all the possible next states",
                 " `v` should be equal to 1, for all the (d + 1) matrices",
                 " p_(i/d).\nWhat follows are the cases that violate",
                 " this principle:\n",
                 paste("\n", possible_cases[logical_vector]))
        } else if (any(same <- sapply(1:D, function(d1) {
            sapply(1:D, function(d2) {
                if (d1 != d2) {
                    return(all_equal(p_dist[, , d1], p_dist[, , d2]))
                }
                FALSE
            })}))) {
            # case where p is NOT drifting, but given as drifting.
            which_same <- which(same)
            arrays <- array(names_i_d(degree, 'p'), dim = c(D, D))
            same_arrays <- arrays[which_same] # always even in number...
            odd <- same_arrays[x <- seq(1, length(same_arrays), by = 2)]
            even <- same_arrays[x + 1]
            same_pairs <- paste0("\n\n", odd, " and ", even, collapse = ',\n ')
            stop("\nThe values given in `p_dist`",
                 " are the same for the arrays:",
                 same_pairs, ".")
        }
    } else {
        # Case when p is not drifting, a.k.a. `p_is_drifting` = FALSE.
        if (!is_double_matrix(p_dist)) {
            stop("\nSince `p_is_drifting` = FALSE, `p_dist` should be",
                 " a square matrix with dimensions: \n(s, s)")
        } else if (!is_prob(p_dist)) {
            non_prob <- which(!is_prob_vector(p_dist))
            stop("\nArray `p_drift` contains values that are not ",
                 "probabilities:\n",
                 paste0(p_dist[non_prob], collapse = ", "))
        } else if (!all((dimension <- dim(p_dist)) == c(s, s))) {
            stop("\nSince `p_is_drifting` = FALSE, `p_dist` should be",
                 " a square matrix with dimensions: \n(s, s) = (",
                 paste(c(s,s), collapse = ', '),
                 ") when it has dimensions: (",
                 paste(dimension, collapse = ', '), ").")
        } else if (!all_equal(sum(p_dist *  diag(s)), 0)) {
            stop("\nThe diagonal values of `p_dist` should all be equal to 0.")
        } else if (!all_equal(as.vector(p_notdrift_rowsum <- rowSums(p_dist)),
                              rep(1, s))) {
            # Case where the sums over v for every u are not equal to 1.
            logical_vector <- sapply(p_notdrift_rowsum,
                                     function(u) !all_equal(u, 1))
            stop("\nFor the transition probability matrice contained in ",
                 "`p_dist`the probabilities of transistioning from ",
                 "previous state `u` over all the possible next states",
                 "`v` should be equal to 1. \nWhat follows are ",
                 "the states  that violate this principle.\n",
                 paste("\n", states[logical_vector]))
        }
    }
    TRUE
}

# ==============================================================================
# Functions specific to the nonparametric object.
# ==============================================================================
valid_fdist_nonparametric <- function(f_dist, states, s, degree,
                                      f_is_drifting, k_max) {
    # '''
    #    Intended for use in `dsmm_nonparametric`.
    # '''
    D <- degree + 1L
    if (f_is_drifting) {
        # Case when f is drifting, a.k.a. `f_is_drifting` == TRUE.
        if (!is_double_array(f_dist)) {
            stop("\nSince `f_is_drifting` is TRUE, `f_dist` should be a ",
                 "numeric array with dimensions:\n(s, s, k_max, degree + 1).")
        } else if (!is_prob(f_dist)) {
            non_prob <- which(!is_prob_vector(f_dist))
            stop("\nArray `f_drift` contains values that are not ",
                 "probabilities:\n",
                 paste0(f_dist[non_prob], collapse = ", "))
        } else if (!all((dimension <- dim(f_dist)) == c(s, s, k_max, D))) {
            stop("\nSince `f_is_drifting` is TRUE, the `f_dist` should be an ",
                 "array with dimensions:\n(s, s, k_max, degree + 1) = (",
                 paste(c(s, s, k_max, D), collapse = ', '),
                 '), when it has',
                 " dimensions (", paste(dimension, collapse = ', '),
                 ")\nThese dimensions",
                 " should correspond to previous states `u`,",
                 " current states `v`,\n",
                 "maximum sojourn time of `k_max` and degree + 1.")
        } else if (!sum(f_dist * array(diag(s),
                                       dim = c(s, s, k_max, D))) ==  0) {
            # Diagonal values are not equal to 0.
            stop("\nThe diagonal values of `f_dist` should all be equal to 0.")
        } else if (!all_equal_numeric(
            f_drift_l_sum <- apply(f_dist, c(1,2,4), sum),
            no_diag <- array(array(1, dimension[1:2]) -
                             base::diag(dimension[1]),
                             dim = dimension[-3]))
        ) { # Sums over l are not equal to 1.
            array_names <- sapply(
                names_i_d(as.integer(degree), "f"), function(d)
                    sapply(states, function(v)
                        sapply(states, function(u)
                            paste(c(paste("u =", u),
                                    paste("v =", v),
                                    paste("Array =", d)),
                                  collapse = ', '))))
            diffs <- sapply(f_drift_l_sum - no_diag, function(diff)
                !all_equal_numeric(diff, 0))
            stop('\nThe sums over l of `f_dist` are not equal to 1 for the ',
                'following cases of u, v and array:\n',
                paste0(array_names[diffs], collapse = '\n'))
        } else if (any(same <- sapply(1:D, function(d1) {
            sapply(1:D, function(d2) {
                if (d1 != d2) {
                    return(all_equal(f_dist[, , , d1], f_dist[, , , d2]))
                }
                FALSE
            })}))) {
            # case where f is NOT drifting, but given as drifting.
            which_same <- which(same)
            arrays <- array(names_i_d(degree, 'f'), dim = c(D, D))
            same_arrays <- arrays[which_same] # always even in number...
            odd <- same_arrays[x <- seq(1, length(same_arrays), by = 2)]
            even <- same_arrays[x + 1]
            same_pairs <- paste0("\n\n", odd, " and ", even, collapse = ',\n ')
            stop("\nThe values given in `f_dist` are the same for the arrays:",
                 same_pairs, ".")
        }
    } else {
        # Case when f is not drifting, a.k.a. `f_is_drifting` = FALSE.
        if (!is_double_array(f_dist)) {
            stop("\nSince `f_is_drifting` is FALSE, `f_dist` should be a ",
                 "numeric array with dimensions:\n(s, s, k_max).")
        } else if (!is_prob(f_dist)) {
            non_prob <- which(!is_prob_vector(f_dist))
            stop("\nArray `f_drift` contains values that are not ",
                 "probabilities:\n",
                 paste0(f_dist[non_prob], collapse = ", "))
        } else if (!all((dimension <- dim(f_dist)) == c(s, s, k_max))) {
            stop("\nSince `f_is_drifting` is FALSE, the `f_dist` should be",
                 " an array with dimensions:\n",
                 "(s, s, k_max) = (", paste(c(s, s, k_max), collapse = ', '),
                 ") when it has dimensions of (",
                 paste(dimension, collapse = ", "), ").\n",
                 "These dimensions should correspond to previous states `u`, ",
                 "current states `v`, and maximum sojourn time of `k_max`.")
        } else if (!all_equal_numeric(sum(f_dist *
                                          array(diag(s), dimension)), 0)) {
            stop("\nThe diagonal values of `f_dist` should all be equal to 0.")
        } else if (!all_equal_numeric(
            f_notdrift_l_sum <- apply(f_dist, c(1,2), sum),
            no_diag <- array(1, dimension[1:2]) - base::diag(dimension[1]))
        ) {
            array_names <- sapply(states, function(v)
                sapply(states, function(u)
                    paste(c(paste("u =", u),
                            paste("v =", v)),
                          collapse = ', ')))
            diffs <- sapply(f_notdrift_l_sum - no_diag, function(diff)
                !all_equal_numeric(diff, 0))
            stop('\nThe sums over l of `f_dist` are not equal to 1 for the ',
                 'following cases of u and v:\n',
                 paste0(
                     array_names[diffs],
                     collapse = '\n')
            )
        }
    }
    TRUE
}

# Specific to the parametric object. -------------------------------------------
# Check the validity of the f distribution given.
valid_fdist_parametric <- function(fdist, params, degree, s, f_is_drifting,
                                   states) {
    # '''
    #    * To be used in `dsmm_parametric` as a check for valid parameters given.
    #    * This is used before `get_fdist_parametric`.
    #    * Also used in `fit_dsmm()` for the parametric estimation, for which
    #       we define a special case when `params` attribute is equal to NULL.
    #       In this case, we check specifically for when `fdist` is a character
    #       array.
    # '''
    # fdist check.
    # both drift & nondrift cases.
    D <- degree + 1L
    if (!is.array(fdist) && !is.character(fdist)) {
        stop("\n`f_dist` should be a character array.")
    } else if (all(is.na(fdist))) {
        stop("\n`f_dist` should not have all if its values equal to NA.")
    }  else if (!all((unqf <- unique(fdist[!is.na(fdist)])) %in%
                     (dist_vector <- c('dweibull', 'geom',
                                       'nbinom', 'pois', 'unif')))) {
        stop('\n`f_dist` should only contain the following values: \nNA, "',
             paste0(dist_vector, collapse = '", "'), '", when it contains:\n',
             unqf)
    }
    # Checking f.
    if (f_is_drifting) {
        if (!all((tmp_d <- dim(fdist)) == (real_dim <- c(s, s, D)))) {
            stop("\nSince `f_is_drifting` = TRUE, `f_dist` should have ",
                 "dimensions: \n(s, s, degree + 1) = (",
                 paste(real_dim, collapse = ', '),
                 ') when it has dimensions: \n(',
                 paste(tmp_d, collapse = ', '), ").")
        } else if (!all(apply(fdist, c(3), function(matrix_uv)
            all(is.na(diag(matrix_uv)))))) {
            stop("\n`f_dist` should have all the diagonal values equal to NA.")
        }
    } else {
        if (!all((tmp_d <- dim(fdist)) == (real_dim <- c(s, s)))) {
            stop("\nSince `f_is_drifting` = FALSE,",
                 " `f_dist` should be a square matrix",
                 " with dimensions: \n(s, s) = (",
                 paste(real_dim, collapse = ', '),
                 ') when it has dimensions: \n(',
                 paste(tmp_d, collapse = ', '), ").")
        } else if (!all(is.na(diag(fdist)))) {
            stop("\n`f_dist` should have all the diagonal values equal to NA.")
        }
    }
    # Parametric ESTIMATION case   -   `fit_dsmm()`.
    if (is.null(params)) {
        return(TRUE) # checks for parameters do not occur.
    }
    # Parametric OBJECT case   -   `parametric_dsmm()`.
    if (!is_double_array(params)) {
        stop("\n`f_dist_parameters` should be a array with double ",
             "and NA values.")
    }  else if (all(is.na(params))) {
        stop("\n`f_dist_parameters` should not have all if its values",
             "equal to NA.")
    }
    # Checking parameters
    if (f_is_drifting) {
        if (!all((par_tmp_d <- dim(params)) == (par_real_d <- c(s, s, 2L, D)))) {
            stop("\nSince `f_is_drifting` = TRUE, `f_dist_parameters`",
                 " should have dimensions:",
                 "\n(s, s, 2, degree + 1) = (",
                 paste(par_real_d, collapse = ', '),
                 ') when it has dimensions: \n(',
                 paste(par_tmp_d, collapse = ', '), ").")
        } else if (!all(apply(params, c(3, 4), function(matrix_uvd)
            all(is.na(diag(matrix_uvd)))))) {
            stop("\n`f_dist_parameters` should have all the diagonal",
                 " values equal to NA.")
        } else if (any(same <- sapply(1:D, function(d1) {
            sapply(1:D, function(d2) {
                if (d1 != d2) {
                    return(all_equal(params[, , ,d1], params[, , ,d2]))
                }
                FALSE
            })}))) {
            # case where f is NOT drifting, but given as drifting.
            which_same <- which(same)
            arrays <- array(names_i_d(degree, 'f_pars'), dim = c(D, D))
            same_arrays <- arrays[which_same] # always even in number...
            odd <- same_arrays[x <- seq(1, length(same_arrays), by = 2)]
            even <- same_arrays[x + 1]
            same_pairs <- paste0("\n\n", odd, " and ", even, collapse = ',\n ')
            stop("\nThe values given in `f_dist_pars`",
                 " are the same for the arrays:",
                 same_pairs, ".")
        }
        # Check for the validity of the parameters, for every i, j, d.
        is_valid_f_param <- sapply(1:D, function(d) {
            dist_ij <- fdist[,,d]
            params_ij <- params[,,,d]
            sapply(1:s, function(j) {
                dist_i <- dist_ij[,j]
                params_i <- params_ij[,j,]
                sapply(1:s, function(i) {
                    row <- c(dist_i[i], params_i[i,])
                    valid_parameters(row = row, i = i, j = j, d = d,
                                     degree = degree, states = states)
                })
            })
        })
    } else {
        if (!all((par_tmp_d <- dim(params)) == (par_real_d <- c(s, s, 2L)))) {
            stop("\nSince `f_is_drifting` = FALSE, ",
                 "`f_dist_parameters` should be an ",
                 "array with dimensions: \n(s, s, 2) = (",
                 paste(par_real_d, collapse = ', '),
                 ') when it has dimensions: \n(',
                 paste(par_tmp_d, collapse = ', '), ").")
        } else if (!all(
            apply(params, c(3),
                  function(matrix_uvd) all(is.na(diag(matrix_uvd)))))) {
            stop("\n`f_dist_parameters` should have all the diagonal",
                 " values equal to NA.")
        }
        is_valid_f_param <- sapply(1:s, function(j) {
            dist_i <- fdist[,j]
            params_i <- params[,j,]
            sapply(1:s, function(i) {
                row <- c(dist_i[i], params_i[i,])
                valid_parameters(row = row, i = i, j = j, d = 0,
                                 degree = degree, states = states)
            })
        })
    }
    all(is_valid_f_param)
}

# These functions are used in conjunction, therefore they are packed together.
# Check the ''row'' values.
valid_parameters <- function(row, i, j, d, degree, states) {
    # '''
    #    This is to be used in function `valid_fdist_parametric`.
    #    `row` argument is supposed to have 3 values:
    #    The first one : a string in ('dweibull', 'geom', 'pois', 'nbinom',
    #    'unif').
    #    The second & third one : a numeric value to be given to the
    #    corrresponding distributions.
    # '''
    if (all(is.na(row))) return(TRUE)
    distr <- row[1]
    p1 <- as.numeric(row[2])
    p2 <- as.numeric(row[3])
    msg <- paste0("\nCurrently for the sojourn time distribution ",
                  names_i_d(d = degree, kernel_name = "f")[d],
                  ", for the previous state u = ", states[i],
                  " and the next state v = ", states[j],
                  ", we have that the first parameter is equal to ", p1,
                  " and the second parameter is equal to ", p2, ".")
    if (distr == "unif") {
        if (is.na(p1) || !(is.na(p2))) {
            stop("\nFor Uniform distributions, only the first parameter",
                 " must be specified.", msg)
        } else if (!((p1 > 0) && ((p1 %% 1) == 0))) {
            stop("\nFor Uniform distributions, the value of the parameter",
                 " must be a positive integer.", msg)
        }
    } else if (distr == "geom") {
        if (is.na(p1) || !(is.na(p2))) {
            stop("\nFor Geometric distributions, the first parameter must",
                 " be specified ",
                 "and the second parameter must be NA.", msg)
        } else if (p1 <= 0 || p1 >= 1) {
            stop("\nFor Geometric distributions, the value of the parameter",
                 " must be between [0, 1] (probability of success).", msg)
        }
    } else if (distr == "pois") {
        if (is.na(p1) || !(is.na(p2))) {
            stop("\nFor Poisson distributions, the first parameter must",
                 " be specified and the second parameter must be NA.", msg)
        } else if (p1 <= 0) {
            stop("\nFor Poisson distributions, the first parameter ",
                 "specified must be a positive number.", msg)
        }
    } else if (distr == "dweibull") {
        if (anyNA(c(p1, p2))) {
            stop("\nFor Discrete Weibull distributions, both ",
                 "parameters must be specified.", msg)
        } else if (p1 <= 0 || p1 >= 1 || p2 <= 0) {
            stop("\nFor Discrete Weibull distributions, the value ",
                 "of the first parameter must be between [0, 1],",
                 " and the second parameter must be a positive number.",
                 msg)
        }
    } else if (distr == "nbinom") {
        if (anyNA(c(p1, p2))) {
            stop("\nFor Negative Binomial distributions, both parameters ",
                 "must be specified.", msg)
        } else if (p1 <= 0 || p2 <= 0 || p2 >= 1) {
            stop("\nFor Negative Binomial distributions, the first parameter ",
                 "must be a positive number (parameter of overdispersion), ",
                 "and the value of the second parameter must ",
                 "be between [0, 1] ",
                 "(the parameter is the probability of success).", msg)
        }
    }
    TRUE
}

# Get the numerical simulated values for a big limit of `klim` sojourn times.
get_fdist_parametric <- function(fdist, params, klim) {
    # '''
    #    This is intended to be used in S3 methods `get_kernel` and `simulate`,
    #    for the class `dsmm_parametric`.
    #    Therefore, no check is required for `fdist` and `params`.
    # '''
    m <- matrix(c(fdist, params[,,1,], params[,,2,]), ncol = 3)
    # No `nrow` given, because by not defining the number of rows,
    # we simultaneously define the drifting & non-drifting cases.
    ans <- apply(m, c(1), function(row, klim) get_dist_param(row, klim),
                 k = klim)
    return(ans)
}

# Get the ''row'' values.
get_dist_param <- function(row, k) {
    # '''
    #    This is a helper function intended for use in `get_fdist_parametric`.
    # '''
    if (all(is.na(row))) return(rep(0, k))
    dist <- row[1]
    par1 <- as.numeric(row[2])
    par2 <- as.numeric(row[3]) # if it is NA, it is not given as a parameter.
    if (dist == 'unif') {
        unif_value <- 1/par1
        distklim <- sapply(1:k, function(x) if (x <= par1) unif_value else 0)
    } else if (dist == 'geom') {
        distklim <- dgeom(0:(k - 1), prob = par1)
    } else if (dist == 'pois') {
        distklim <- dpois(0:(k - 1), lambda = par1)
    } else if (dist == 'nbinom') {
        distklim <- dnbinom(0:(k - 1), size = par1, prob = par2)
    } else if (dist == 'dweibull') {
        distklim <- DiscreteWeibull::ddweibull(1:k, q = par1,
                                               beta = par2,
                                               zero = FALSE)
        # We keep zero = false because of the parametric estimation.
    }
    if ((sumd <- sum(distklim)) < 1) {
        distklim[k] <- 1 - sumd
    }
    distklim
}

valid_estimation_fit <- function(estimation, degree, states, s, k_max,
                                 p_is_drifting, f_is_drifting, fdist, pdist) {
    # '''
    #    This function is acting on the OBJECT `dsmm_fit_nonparametric` and
    #    `dsmm_fit_parametric`. Specifically, it is used when printing through
    #    `check_attributes()`, in order to prevent someone from defining
    #    their own object in a `vile` manner.
    # '''
    if (estimation != 'nonparametric' && estimation != 'parametric') {
        stop("`estimation` attribute should be either 'nonparametric' or ",
             " 'parametric'.")
    }
    stopifnot(valid_p_dist(p_dist = pdist,
                           states = states,
                           s = s, degree = degree,
                           p_is_drifting = p_is_drifting))
    if (estimation == 'nonparametric') {
        stopifnot(valid_fdist_nonparametric(f_dist = fdist,
                                            states = states,
                                            s = s,
                                            degree = degree,
                                            f_is_drifting = f_is_drifting,
                                            k_max = k_max))
    } else {
        stopifnot(valid_estimation(estimation = estimation,
                                   fpar = fdist, s = s,
                                   f_is_drifting = f_is_drifting,
                                   degree = degree, states = states),
                  valid_fdist_parametric(fdist = fdist,
                                         params = NULL,
                                         s = s, degree = degree,
                                         f_is_drifting = f_is_drifting))
    }
    TRUE
}


# ------------------------------------------------------------------------------
# Regarding the drifting kernels.
poly_coeff <- function(degree) {
    # '''
    #    Returns a matrix with the coefficients of the polynomials `A_i`.
    # '''
    x <- (s0 <- seq(0, degree)) / degree
    x <- sapply(x, function(v) sapply(s0, function(i) v ** i))
    ans <- solve(x)
    colnames(ans) <- paste0("t^", s0)
    rownames(ans) <- paste0("A_", s0)
    ans
}

names_i_d <- function(d, kernel_name = "q") {
    # '''
    #    Generate names for the drifting cases. This is not exported.
    # '''
    kernel_name <- paste0(kernel_name, "_")
    kernel_0 <- paste0(kernel_name, "0")
    kernel_1 <- paste0(kernel_name, "1")
    if (is_integer(d) && d > 1) {
        return(c(paste0(kernel_name, "0"),
                 paste0(kernel_name, "(", seq(1, d-1), "/", d, ")"),
                 paste0(kernel_name, "1"))
        )
    } else if (d == 1) {
        return(c(kernel_0, kernel_1))
    } else if (d == 0) {
        return(kernel_name)
    }
}

seq_states_to_id <- function(seq, states) {
    # Transforming a sequence of states into numbers.
    list_which <- sapply(states, function(u) which(seq == u) )
    for (i in seq_along(states)) {
        seq[list_which[[i]]] <- i
    }
    seq
}

# Get different values needed to compute the different distributions. ----------
get_1_u <- function(id_seq, l, n, id_states, s) {
    # '''
    #    Used in `fit_dsmm()` to get the indicator function 1_u.
    #    `id_seq` is the sequence of the embedded Markov chain where we
    #    substituted the states for positive integers.
    #    `id_states` is seq(1, s).
    # '''
    zero_v <- rep(0, n)
    vector_1_u <- sapply(id_states,
                         function(u, seq_u, zero_v) {
                             zero_v[which(u == seq_u)] <- 1
                             return(zero_v)},
                         seq_u = id_seq[-l], # minus the last state
                         zero_v = zero_v)
    rownames(vector_1_u) <- paste0("t = ", 1:n)
    vector_1_u
}

get_1_uvl <- function(id_seq, l, n, X, k_max, id_states, s) {
    # '''
    #    Used in `fit_dsmm()` to get the indicator function 1_uvl.
    #    `id_seq` is the sequence of the embedded Markov chain where we
    #    substituted the states for positive integers.
    #    `id_states` is seq(1, s).
    # '''
    seq_uvl <- paste(id_seq[-l], c(id_seq[-1]), X[-l])
    Es <- rep(id_states, s)
    pEs <- paste(Es, sort(Es))
    possible_uvl_states <- paste(rep(pEs, k_max),
                                 sort(rep(1:k_max, s*s)))
    zero_v <- rep(0, n)
    vector_1_uvl <- sapply(X = possible_uvl_states,
                           FUN = function(uvl, seq_uvl, zero_v) {
                               zero_v[which(uvl == seq_uvl)] <- 1
                               return(zero_v)},
                           seq_uvl = seq_uvl, zero_v = zero_v)
    rownames(vector_1_uvl) <- paste0("t = ", 1:n)
    vector_1_uvl
}

get_model <- function(p_is_drifting, f_is_drifting) {
    # '''
    #    Should only be used AFTER `valid_model(...)` function has been used!
    #    Used in `fit_dsmm`, `parametric_dsmm` and `nonparametric_dsmm`.
    # '''
    Model_1 <- p_is_drifting && f_is_drifting
    Model_2 <- p_is_drifting && !f_is_drifting
    Model_3 <- !p_is_drifting && f_is_drifting
    model_id <- which(c(Model_1, Model_2, Model_3) == TRUE)
    model <- c("Model_1", "Model_2", "Model_3")[model_id]
    model
}

get_vl <- function(t, u, vl_names, kernel, states) {
    # '''
    #    This function is intended to be used in the generic function
    #    `simulate.dsmm`.
    #    It returns a vector of two character values:
    #    The first is the name of the next state `v`, and the second is the
    #    sojourn time `spend on previous state
    # '''
    vl <- sample(vl_names, prob = kernel[u, , , t], size = 1)
    #  ## Because we have `size = 1`,
    #  ## the `replace` argument is irrelevant
    #  ## See `sample()` for more details.
    ans <- strsplit(vl, " ")[[1]] # Splitting the space added in `paste()`.
    if ((l_ans <- length(ans)) > 2) {
        # If the state contains spaces, we take all but the last values of
        # `ans` and paste them together. The last value of `ans` is the
        # corresponding sojourn time of the state.
        v <- paste(ans[1:(l_ans-1)], collapse = " ")
        l <- ans[l_ans]
        ans <- c(v, l)
    }
    ans
}

get_A_i <- function(degree, model_size) {
    # '''
    #    For intended use in `fit_dsmm` and `get_kernel`.
    # '''
    n <- model_size
    coeff <- poly_coeff(degree = degree)
    t_n <- sapply((t <- 0:n)/n , function(tn) tn**(0:degree))
    res <- coeff %*% t_n
    colnames(res) <- paste0("t = ", t)
    return(res)
}

# Normalization of f.
get_f <- function(v) (a <- abs(v)) / sum(a)

# Check the validity of the computed kernel and return it.
get_valid_kernel <- function(Ji, Ai, s, n, k_max, states) {
    N <- n + 1
    ans <- Ji %*% Ai # ans has dimensions of (s*s*k_max, n+1)
    kernel <- array(ans, dim = c(s, s, k_max, N))
    if (!all_equal_numeric(apply(kernel, c(1,4), sum), 1)) {
        stop("\nThe sums of the kernel are not equal to 1.")
    }
    if (!is_prob(kernel)) {
        kernel_pos <- apply(kernel, c(1,4), get_f)
        dim(kernel_pos) <- c(s, k_max, s, N)
        kernel <- aperm(kernel_pos, c(3,1,2,4))
    }
    dimnames(kernel) <- list(
        as.list(states),
        as.list(states),
        as.list(paste0("l = ", 1:k_max)),
        as.list(paste0("t = ", 0:n))
    )
    kernel
}


# Random sequence creation
#' @title Simulate a sequence for states of choice.
#' @description This is a wrapper function around \code{sample()}.
#'
#' @param states Character vector of unique values. If given the value
#' "DNA" the values c("a", "c", "g", "t") are given instead.
#' @param len Optional. Positive integer with the default value
#' equal to 5000.
#' @param probs Optional. Numeric vector with values interpreted as
#' probabilities for each of the states in \code{states}. Default value
#' is equal to 1 over the number of states given, for every state.
#' @param seed Optional. Object specifying the initialization of the random
#' number generator (see more in \code{\link[base:set.seed]{set.seed}}).
#'
#' @seealso
#' For the simulation of a sequence with a Drifting semi-Markov kernel:
#' \link{simulate.dsmm}.
#'
#' The original function: \code{\link[base:sample]{sample}}.
#'
#' About random number generation in R: \code{\link[base:RNG]{RNG}}.
#'
#' For the theoretical background of Drifting semi-Markov models: \link{dsmmR}.
#'
#' @return A character sequence of length \code{len}.
#' @export
#' @examples
#' # This is equal to having the argument `probs = c(1/4, 1/4, 1/4, 1/4)`.
#' rand_dna_seq <- create_sequence(states = "DNA")
#' table(rand_dna_seq)
#'
#' random_letters <- sample(letters, size = 5, replace = FALSE)
#' rand_dna_seq2 <- create_sequence(
#'     states = random_letters,
#'     probs = c(0.6, 0.3, 0.05, 0.025, 0.025),
#'     len = 10000)
#' table(rand_dna_seq2)
create_sequence <- function(states, len = 5000, probs = NULL, seed = NULL) {
    # '''
    #    This is a convenience function intended for fast sequence simulation.
    # '''
    if (missing(states)) {
        stop("\nPlease input a character vector of states.")
    } else if (identical(states, "DNA")) {
        states <- c("a", "c", "g", "t")
    }
    if (!is.null(seed)) {
        set.seed(seed)
    }
    stopifnot(valid_states(states), is_integer(len))
    s <- length(states)
    if (is.null(probs)) {
        overs <- 1/s
        probs <- rep(overs, s)
    }
    stopifnot(is_prob(probs))
    seq <- sample(x = states, size = len, replace = TRUE, prob = probs)
    seq
}

