# We can also define states in a flexible way, including spaces.
states <- c("Dollar $", " /1'2'3/ ", " Z E T A ", "O_M_E_G_A")
s <- length(states)
d <- 1


# `p_dist` has dimensions of: (s, s, d + 1).
# Sums over v must be 1 for all u and i = 0, ..., d.

# First matrix.
p_dist_1 <- matrix(c(0,   0.1, 0.4, 0.5,
                     0.5, 0,   0.3, 0.2,
                     0.3, 0.4, 0,   0.3,
                     0.8, 0.1, 0.1, 0),
                   ncol = s, byrow = TRUE)

# Second matrix.
p_dist_2 <- matrix(c(0,   0.3, 0.6, 0.1,
                     0.3, 0,   0.4, 0.3,
                     0.5, 0.3, 0,   0.2,
                     0.2, 0.3, 0.5, 0),
                   ncol = s, byrow = TRUE)

# `f_dist` has dimensions of: (s, s, d + 1).
# First matrix.
f_dist_1 <- matrix(c(NA,         "unif", "dweibull", "nbinom",
                     "geom",      NA,    "pois",     "dweibull",
                     "dweibull", "pois",  NA,        "geom",
                     "pois",      NA,    "geom",      NA),
                   nrow = s, ncol = s, byrow = TRUE)


# Second matrix.
f_dist_2 <- matrix(c(NA,     "pois", "geom", "nbinom",
                     "geom",  NA,    "pois", "dweibull",
                     "unif", "geom",  NA,    "geom",
                     "pois", "pois", "geom",  NA),
                   nrow = s, ncol = s, byrow = TRUE)

# `f_dist_pars` has dimensions of: (s, s, 2, d + 1).
# First array of coefficients, corresponding to `f_dist_1`.
# First matrix.
f_dist_1_pars_1 <- matrix(c(NA,  5,  0.4, 4,
                            0.7, NA, 5,   0.6,
                            0.2, 3,  NA,  0.6,
                            4,   NA, 0.4, NA),
                          nrow = s, ncol = s, byrow = TRUE)
# Second matrix.
f_dist_1_pars_2 <- matrix(c(NA,  NA, 0.2, 0.6,
                            NA,  NA, NA,  0.8,
                            0.6, NA, NA,  NA,
                            NA,  NA, NA,  NA),
                          nrow = s, ncol = s, byrow = TRUE)

# Second array of coefficients, corresponding to `f_dist_2`.
# First matrix.
f_dist_2_pars_1 <- matrix(c(NA,  6,   0.4, 3,
                            0.7, NA,  2,   0.5,
                            3,   0.6, NA,  0.7,
                            6,   0.2, 0.7, NA),
                          nrow = s, ncol = s, byrow = TRUE)
# Second matrix.
f_dist_2_pars_2 <- matrix(c(NA, NA, NA, 0.6,
                            NA, NA, NA, 0.8,
                            NA, NA, NA, NA,
                            NA, NA, NA, NA),
                          nrow = s, ncol = s, byrow = TRUE)

test_that("parametric_dsmm(); p and f are drifting.", {
    # ---------------------------------------------------------------------------
    # Parametric object for Model 1.
    # ---------------------------------------------------------------------------
    # get `p_dist` as an array of p_dist_1 and p_dist_2.
    p_dist_model_1 <- array(c(p_dist_1, p_dist_2), dim = c(s, s, d + 1))
    # get `f_dist` as an array of `f_dist_1` and `f_dist_2`
    f_dist_model_1 <- array(c(f_dist_1, f_dist_2), dim = c(s, s, d + 1))
    
    f_dist_pars_model_1 <- array(c(f_dist_1_pars_1, f_dist_1_pars_2,
                                   f_dist_2_pars_1, f_dist_2_pars_2),
                                 dim = c(s, s, 2, d + 1))
    
    expect_no_condition(
        obj_par_model_1 <- parametric_dsmm(
            model_size = 10000,
            states = states,
            initial_dist = c(0.8, 0.1, 0.1, 0),
            degree = d,
            p_dist = p_dist_model_1,
            f_dist = f_dist_model_1,
            f_dist_pars = f_dist_pars_model_1,
            p_is_drifting = TRUE,
            f_is_drifting = TRUE
        )
    )
    expect_snapshot(obj_par_model_1)
    
})

test_that("parametric_dsmm(); p is drifting, f is not drifting.", {
    # `p_dist` has the same dimensions as in Model 1: (s, s, d + 1).
    p_dist_model_2 <- array(c(p_dist_1, p_dist_2), dim = c(s, s, d + 1))
    
    # `f_dist` has dimensions of: (s, s).
    f_dist_model_2 <- matrix(c( NA,       "pois",  NA,       "nbinom",
                                "geom",    NA,    "geom",    "dweibull",
                                "unif",   "geom",  NA,       "geom",
                                "nbinom", "unif", "dweibull", NA),
                             nrow = s, ncol = s, byrow = TRUE)
    
    # `f_dist_pars` has dimensions of: (s, s, 2),
    #  corresponding to `f_dist_model_2`.
    
    # First matrix.
    f_dist_pars_1_model_2 <- matrix(c(NA,  0.2, NA,  3,
                                      0.2, NA,  0.2, 0.5,
                                      3,   0.4, NA,  0.7,
                                      2,   3,   0.7, NA),
                                    nrow = s, ncol = s, byrow = TRUE)
    
    # Second matrix.
    f_dist_pars_2_model_2 <- matrix(c(NA,  NA, NA,  0.6,
                                      NA,  NA, NA,  0.8,
                                      NA,  NA, NA,  NA,
                                      0.2, NA, 0.3, NA),
                                    nrow = s, ncol = s, byrow = TRUE)
    
    
    # Get `f_dist_pars`.
    f_dist_pars_model_2 <- array(c(f_dist_pars_1_model_2,
                                   f_dist_pars_2_model_2),
                                 dim = c(s, s, 2))
    
    
    # ---------------------------------------------------------------------------
    # Parametric object for Model 2.
    # ---------------------------------------------------------------------------
    
    expect_no_condition(
        obj_par_model_2 <- parametric_dsmm(
            model_size = 10000,
            states = states,
            initial_dist = c(0.8, 0.1, 0.1, 0),
            degree = d,
            p_dist = p_dist_model_2,
            f_dist = f_dist_model_2,
            f_dist_pars = f_dist_pars_model_2,
            p_is_drifting = TRUE,
            f_is_drifting = FALSE
        )
    )
    expect_snapshot(obj_par_model_2)

})


test_that("parametric_dsmm(); p is not drifting, f is drifting.", {
    # `p_dist` has dimensions of: (s, s).
    p_dist_model_3 <- matrix(c(0,   0.1,  0.3,  0.6,
                               0.4, 0,    0.1,  0.5,
                               0.4, 0.3,  0,    0.3,
                               0.9, 0.01, 0.09, 0),
                             ncol = s, byrow = TRUE)
    
    # `f_dist` has the same dimensions as in Model 1: (s, s, d + 1).
    f_dist_model_3 <- array(c(f_dist_1, f_dist_2), dim = c(s, s, d + 1))
    
    
    # `f_dist_pars` has the same dimensions as in Model 1: (s, s, 2, d + 1).
    f_dist_pars_model_3 <- array(c(f_dist_1_pars_1, f_dist_1_pars_2,
                                   f_dist_2_pars_1, f_dist_2_pars_2),
                                 dim = c(s, s, 2, d + 1))
    
    # ---------------------------------------------------------------------------
    # Parametric object for Model 3.
    # ---------------------------------------------------------------------------
    expect_no_condition(
        obj_par_model_3 <- parametric_dsmm(
            model_size = 10000,
            states = states,
            initial_dist = c(0.3, 0.2, 0.2, 0.3),
            degree = d,
            p_dist = p_dist_model_3,
            f_dist = f_dist_model_3,
            f_dist_pars = f_dist_pars_model_3,
            p_is_drifting = FALSE,
            f_is_drifting = TRUE
        )
    )
    expect_snapshot(obj_par_model_3)
    
})
