states <- c("Rouen", "Bucharest", "Samos", "Aigio", "Marseille")
sequence <- create_sequence(states, probs = c(0.3, 0.1, 0.1, 0.3, 0.2))
obj_model <- fit_dsmm(
    sequence = sequence,
    states = states,
    degree = 3,
    f_is_drifting = FALSE,
    p_is_drifting = TRUE
)

test_that("get_kernel()", {
    
    klim <- 10
    kernel <- get_kernel(obj_model, klim = klim)
    
    expect_type(kernel, "double")
    
    # Each of the first two dimensions should be equal to the length of states.
    expect_identical(
        dim(kernel)[1:2], rep(length(states), 2)
    )
    
    # Third dimension of kernel should be less or equal to `klim`.
    expect_lte(
        dim(kernel)[3], klim
    )
    
    # Fourth dimension of kernel should be less or equal to `length(sequence)`
    # (Equal is extremely unlikely)
    expect_lte(
        dim(kernel)[4], length(sequence)
    )
    
    # Probabilities should sum to 1 over v, l.
    expect_equal(
        kernel_sum <- apply(kernel, c(1, 4), sum),
        array(1, dim = dim(kernel_sum), dimnames = dimnames(kernel_sum))
    )
    
})

test_that("simulate()", {
    
    max_seq_length <- 10
    sim_seq <- simulate(obj_model, max_seq_length = max_seq_length)
    
    # `length(sim_seq)` should be less or equal to `max_seq_length`
    expect_lte(length(sim_seq), max_seq_length)
    
    # `sim_seq` should be a subset of `states`
    expect_in(sim_seq, states)
    
})




















