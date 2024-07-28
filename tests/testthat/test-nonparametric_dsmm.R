# Setup.
states <- c("AA", "AC", "CC")
s <- length(states)
d <- 2
k_max <- 3

# `p_dist` has dimensions of: (s, s, d + 1).
# Sums over v must be 1 for all u and i = 0, ..., d.
p_dist_1 <- matrix(c(0,   0.1, 0.9,
                     0.5, 0,   0.5,
                     0.3, 0.7, 0),
                   ncol = s, byrow = TRUE)

p_dist_2 <- matrix(c(0,   0.6, 0.4,
                     0.7, 0,   0.3,
                     0.6, 0.4, 0),
                   ncol = s, byrow = TRUE)

p_dist_3 <- matrix(c(0,   0.2, 0.8,
                     0.6, 0,   0.4,
                     0.7, 0.3, 0),
                   ncol = s, byrow = TRUE)

# `f_dist` has dimensions of: (s, s, k_max, d + 1).
# First f distribution. Dimensions: (s, s, k_max).
# Sums over l must be 1, for every u, v and i = 0, ..., d.
f_dist_1_l_1 <- matrix(c(0,   0.2, 0.7,
                         0.3, 0,   0.4,
                         0.2, 0.8, 0),
                       ncol = s, byrow = TRUE)

f_dist_1_l_2 <- matrix(c(0,   0.3,  0.2,
                         0.2, 0,    0.5,
                         0.1, 0.15, 0),
                       ncol = s, byrow = TRUE)

f_dist_1_l_3 <- matrix(c(0,   0.5,  0.1,
                         0.5, 0,    0.1,
                         0.7, 0.05, 0),
                       ncol = s, byrow = TRUE)
# Get f_dist_1
f_dist_1 <- array(c(f_dist_1_l_1, f_dist_1_l_2, f_dist_1_l_3),
                  dim = c(s, s, k_max))

# Second f distribution. Dimensions: (s, s, k_max)
f_dist_2_l_1 <- matrix(c(0,   1/3, 0.4,
                         0.3, 0,   0.4,
                         0.2, 0.1, 0),
                       ncol = s, byrow = TRUE)

f_dist_2_l_2 <- matrix(c(0,   1/3, 0.4,
                         0.4, 0,   0.2,
                         0.3, 0.4, 0),
                       ncol = s, byrow = TRUE)

f_dist_2_l_3 <- matrix(c(0,   1/3, 0.2,
                         0.3, 0,   0.4,
                         0.5, 0.5, 0),
                       ncol = s, byrow = TRUE)

# Get f_dist_2
f_dist_2 <- array(c(f_dist_2_l_1, f_dist_2_l_2, f_dist_2_l_3),
                  dim = c(s, s, k_max))

# Third f distribution. Dimensions: (s, s, k_max)
f_dist_3_l_1 <- matrix(c(0,    0.3, 0.3,
                         0.3,  0,   0.5,
                         0.05, 0.1, 0),
                       ncol = s, byrow = TRUE)

f_dist_3_l_2 <- matrix(c(0,   0.2, 0.6,
                         0.3, 0,   0.35,
                         0.9, 0.2, 0),
                       ncol = s, byrow = TRUE)

f_dist_3_l_3 <- matrix(c(0,    0.5, 0.1,
                         0.4,  0,   0.15,
                         0.05, 0.7, 0),
                       ncol = s, byrow = TRUE)

# Get f_dist_3
f_dist_3 <- array(c(f_dist_3_l_1, f_dist_3_l_2, f_dist_3_l_3),
                  dim = c(s, s, k_max))

test_that("nonparametric_dsmm(); p and f are drifting", {
    
    # ---------------------------------------------------------------------------
    # Defining distributions for Model 1 - both p and f are drifting.
    # ---------------------------------------------------------------------------
    # Get `p_dist` as an array of p_dist_1, p_dist_2 and p_dist_3.
    p_dist <- array(c(p_dist_1, p_dist_2, p_dist_3),
                    dim = c(s, s, d + 1))
    
    # Get f_dist as an array of f_dist_1, f_dist_2 and f_dist_3.
    f_dist <- array(c(f_dist_1, f_dist_2, f_dist_3),
                    dim = c(s, s, k_max, d + 1))
    
    # ---------------------------------------------------------------------------
    # Non-Parametric object for Model 1.
    # ---------------------------------------------------------------------------
    expect_no_condition(
        obj_nonpar_model_1 <- nonparametric_dsmm(
            model_size = 8000,
            states = states,
            initial_dist = c(0.3, 0.5, 0.2),
            degree = d,
            k_max = k_max,
            p_dist = p_dist,
            f_dist = f_dist,
            p_is_drifting = TRUE,
            f_is_drifting = TRUE
        )
    )
    
    expect_snapshot(obj_nonpar_model_1)
    
})

test_that("nonparametric_dsmm(); p is drifting, f is not drifting.", {
    # ---------------------------------------------------------------------------
    # Defining Model 2 - p is drifting, f is not drifting.
    # -----------------------i----------------------------------------------------
    
    # p_dist has the same dimensions as in Model 1: (s, s, d + 1).
    p_dist_model_2 <- array(c(p_dist_1, p_dist_2, p_dist_3),
                            dim = c(s, s, d + 1))
    
    # f_dist has dimensions of: (s,s,k_{max}).
    f_dist_model_2 <- f_dist_2
    
    # ---------------------------------------------------------------------------
    # Non-Parametric object for Model 2.
    # ---------------------------------------------------------------------------
    
    expect_no_condition(
        obj_nonpar_model_2 <- nonparametric_dsmm(
            model_size = 10000,
            states = states,
            initial_dist = c(0.7, 0.1, 0.2),
            degree = d,
            k_max = k_max,
            p_dist = p_dist_model_2,
            f_dist = f_dist_model_2,
            p_is_drifting = TRUE,
            f_is_drifting = FALSE
        )
    )
    
    expect_snapshot(obj_nonpar_model_2)
    
})

test_that("nonparametric_dsmm(); f is drifting, p is not drifting.", {
    # ---------------------------------------------------------------------------
    # Defining Model 3 - f is drifting, p is not drifting.
    # ---------------------------------------------------------------------------
    
    # `p_dist` has dimensions of: (s, s, d + 1).
    p_dist_model_3 <- p_dist_3
    
    # `f_dist` has the same dimensions as in Model 1: (s, s, d + 1).
    f_dist_model_3 <- array(c(f_dist_1, f_dist_2, f_dist_3),
                            dim = c(s, s, k_max, d + 1))
    
    
    # ---------------------------------------------------------------------------
    # Non-Parametric object for Model 3.
    # ---------------------------------------------------------------------------
    expect_no_condition(
        obj_nonpar_model_3 <- nonparametric_dsmm(
            model_size = 10000,
            states = states,
            initial_dist = c(0.3, 0.4, 0.3),
            degree = d,
            k_max = k_max,
            p_dist = p_dist_model_3,
            f_dist = f_dist_model_3,
            p_is_drifting = FALSE,
            f_is_drifting = TRUE
        )
    )
    
    expect_snapshot(obj_nonpar_model_3)
    
})
#' # ===========================================================================
#' # Using methods for non-parametric objects.
#' # ===========================================================================
#'
#' kernel_parametric <- get_kernel(obj = obj_nonpar_model_3)
#' str(kernel_parametric)
#'
#' sim_seq_par <- simulate(obj_nonpar_model_3, nsim = 50)
#' str(sim_seq_par)
