sequence <- create_sequence("DNA", len = 10000, probs = NULL, seed = 1)
states <- sort(unique(sequence))
degree <- 3

s <- length(states)
f_dist_1 <- matrix(c(NA,         "unif",     "geom",     "pois",
                     "pois",      NA,        "pois",     "geom",
                     "geom",     "pois",      NA,        "geom",
                     "geom",     'geom',     "pois",      NA),
                   nrow = s, ncol = s, byrow = TRUE)
# f_dist_1 <- matrix(c(NA,         "unif",     "dweibull", "nbinom",
#                      "pois",      NA,        "pois",     "dweibull",
#                      "geom",     "pois",      NA,        "geom",
#                      "dweibull", 'geom',     "pois",      NA),
#                    nrow = s, ncol = s, byrow = TRUE)
f_dist <- array(f_dist_1, dim = c(s, s, degree + 1))

test_that("fit_dsmm() non parametric estimation; p and f are drifting", {
    # ===========================================================================
    # Nonparametric Estimation.
    # Fitting a random sequence under distributions of unknown shape.
    # ===========================================================================
    
    
    # ---------------------------------------------------------------------------
    # Both p and f are drifting - Model 1.
    # ---------------------------------------------------------------------------
    
    # This should return without error or warning
    expect_no_condition(
        obj_model_1 <- fit_dsmm(
            sequence = sequence,
            degree = degree,
            f_is_drifting = TRUE,
            p_is_drifting = TRUE,
            states = states,
            initial_dist = "freq",
            estimation = "nonparametric", # default value
            f_dist = NULL # default value
        )
    )
    
    # Sums over v, l should be 1
    expect_equal(
        model_sum <- apply(obj_model_1$J_i, c(1, 4), sum),
        array(1, dim = dim(model_sum), dimnames = dimnames(model_sum))
    )
    
    expect_snapshot(obj_model_1)
    
    
})

test_that("fit_dsmm() non parametric estimation; p is drifting and f is not drifting", {
    
    # ---------------------------------------------------------------------------
    # Fitting the sequence when p is drifting and f is not drifting - Model 2.
    # ---------------------------------------------------------------------------
    
    expect_no_condition(
        obj_model_2 <- fit_dsmm(
            sequence = sequence,
            degree = degree,
            f_is_drifting = FALSE,
            p_is_drifting = TRUE
        )
    )
    
    # Sums over v, l should be 1
    expect_equal(
        model_sum <- apply(obj_model_2$J_i, c(1, 4), sum),
        array(1, dim = dim(model_sum), dimnames = dimnames(model_sum))
    )
    
    expect_snapshot(obj_model_2)
    
})

test_that("fit_dsmm() non parametric estimation; p is not drifting and f is drifting", {
    
    # ---------------------------------------------------------------------------
    # Fitting the sequence when f is drifting and p is not drifting - Model 3.
    # ---------------------------------------------------------------------------
    
    expect_no_condition(
        obj_model_3 <- fit_dsmm(
            sequence = sequence,
            degree = degree,
            f_is_drifting = TRUE,
            p_is_drifting = FALSE
        )
    )
    
    # Sums over v, l should be 1
    expect_equal(
        model_sum <- apply(obj_model_3$J_i, c(1, 4), sum),
        array(1, dim = dim(model_sum), dimnames = dimnames(model_sum))
    )
    
    expect_snapshot(obj_model_3)
    
})

test_that("fit_dsmm() parametric estimation; p and f are drifting", {
    
    # ===========================================================================
    # Parametric Estimation
    # Fitting a random sequence under distributions of known shape.
    # ===========================================================================
    
    # ---------------------------------------------------------------------------
    # Both p and f are drifting - Model 1.
    # ---------------------------------------------------------------------------
    
    set.seed(1)
    expect_no_condition(
        obj_fit_parametric_1 <- fit_dsmm(sequence = sequence,
                                         degree = degree,
                                         f_is_drifting = TRUE,
                                         p_is_drifting = TRUE,
                                         states = states,
                                         initial_dist = 'unif',
                                         estimation = 'parametric',
                                         f_dist = f_dist
        )
    )
    
    # Sums over v, l should be 1
    expect_equal(
        model_sum <- apply(obj_fit_parametric_1$J_i, c(1, 4), sum),
        array(1, dim = dim(model_sum), dimnames = dimnames(model_sum))
    )
    
    expect_snapshot(obj_fit_parametric_1)
    
})

test_that("fit_dsmm() parametric estimation; p is drifting, f is not drifting", {
    
    # ---------------------------------------------------------------------------
    # f is not drifting, only p is drifting - Model 2.
    # ---------------------------------------------------------------------------
    
    set.seed(1)
    expect_no_condition(
        obj_fit_parametric_2 <- fit_dsmm(
            sequence = sequence,
            degree = degree,
            f_is_drifting = FALSE,
            p_is_drifting = TRUE,
            initial_dist = 'unif',
            estimation = 'parametric',
            f_dist = f_dist_1
        )
    )
    
    # Sums over v, l should be 1
    expect_equal(
        model_sum <- apply(obj_fit_parametric_2$J_i, c(1, 4), sum),
        array(1, dim = dim(model_sum), dimnames = dimnames(model_sum))
    )
    
    expect_snapshot(obj_fit_parametric_2)
    
})

test_that("fit_dsmm() parametric estimation; p is not drifting, f is drifting", {
    
    # ---------------------------------------------------------------------------
    # p is not drifting, only f is drifting - Model 3.
    # ---------------------------------------------------------------------------
    
    set.seed(1)
    expect_no_condition(
        obj_fit_parametric_3 <- fit_dsmm(
            sequence = sequence,
            degree = degree,
            f_is_drifting = TRUE,
            p_is_drifting = FALSE,
            initial_dist = 'unif',
            estimation = 'parametric',
            f_dist = f_dist
        )
    )
    
    # Sums over v, l should be 1
    expect_equal(
        model_sum <- apply(obj_fit_parametric_3$J_i, c(1, 4), sum),
        array(1, dim = dim(model_sum), dimnames = dimnames(model_sum))
    )
    
    expect_snapshot(obj_fit_parametric_3)
    
})
