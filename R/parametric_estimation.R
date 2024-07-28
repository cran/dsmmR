# '''
#    This file concerns itself with the helper function for the parametric
#    estimation used in `fit_dsmm()`.
# '''

parametric_estimation <- function(lprobs, dist, kmax, i, j, d, degree, states) {
    # '''
    #     This function is used in `fit_dsmm()` when estimation = "parametric".
    #   * `lprobs` is a vector of probabilities corresponding to the sojourn
    #     times of the sequence.
    #   * `dist` is either NA or one of c('unif', 'pois', 'geom', 'nbinom',
    #     'dweibull').
    #   * `kmax` is the maximum sojourn time of the sequence.
    #   * `i`, `j`, `d` are values corresponding to the current u, v and
    #      sojourn time distribution f_(i/d).
    #   * `degree` is the polynomial degree of the drift.
    # '''
    if (dist == 'unif') {
        posprobs <- sum(lprobs != 0)
        ratio <- c(lprobs[-1], NA) / lprobs
        similar_consec_ratios <- which(ratio <= posprobs & ratio >= 1/posprobs)
        theta <- max(sapply(split(similar_consec_ratios,
                                  cumsum(c(1, diff(similar_consec_ratios) != 1))),
                            length)) + 1
        return(c(theta, NA))
    } else if (dist == 'pois') {
        theta <- sum(0:(kmax - 1) * lprobs)
        return(c(theta, NA))
    } else if (dist == 'geom') {
        theta <- 1 / sum(1:kmax * lprobs)
        return(c(theta, NA))
    } else if (dist == 'nbinom') {
        zerokmax <- 0:(kmax - 1)
        expectation <- sum(zerokmax * lprobs)
        expectationx2 <- sum((zerokmax**2) * lprobs)
        expectationsquared <- expectation**2
        variance <- expectationx2 - expectationsquared
        if (expectation >= variance) {
            stop("The Negative Binomial distribution is not appropriate\n",
                 "for modeling the conditional sojourn time distribution\n",
                 "associated with:\n the current state u = ", states[i],
                 ",\n the next state v = ", states[j], "\n and the distribution ",
                 names_i_d(degree, 'f')[d],
                 ",\nbecause variance = ", variance,
                 " >= ", expectation, " = expectation.")
        }
        phat <- expectation / variance
        alphahat <- expectationsquared / (variance - expectation)
        return(c(alphahat, phat))
    } else if (dist == 'dweibull') {
        if (length(unique(lprobs)) == 1) {
            stop("The Discrete Weibull is not appropriate for modeling\n",
                 "the conditional sojourn time distribution",
                 " describing:\n the previous state u = ", states[i],
                 ",\n the next state v = ", states[j],
                 "\n for the sojourn time distribution ",
                 names_i_d(degree, 'f')[d], ",\nsince the estimation of",
                 " the second parameter beta is not possible.\n\n",
                 "This happens because we only have 1 non-negative value,\n",
                 "and the estimation of beta requires at least 2 non-negative ",
                 "values,\nwhich makes the estimation of beta impossible in",
                 " this case.")
        }
        qhat <- 1 - lprobs[1]
        cumsumf <- 1 - cumsum(lprobs)
        beta_i <- sapply(2:kmax, function(i)
            log(log(cumsumf[i], base = qhat), base = i))
        betahat <- mean(beta_i[is.finite(beta_i)])
        if (is.nan(betahat)) {
            stop("The Discrete Weibull is not appropriate for modeling\n",
                 "the conditional sojourn time distribution",
                 " describing\n the previous state u = ", states[i],
                 ",\n the next state v = ", states[j],
                 "\n for the sojourn time distribution ",
                 names_i_d(degree, 'f')[d], ",\n since the estimation of",
                 " the second parameter beta is not possible.\n\n",
                 "This behaviour is perhaps accounted to the fact that",
                 " beta > 1,\nwhich is generally hard to estimate.")
        }
        return(c(qhat, betahat))
    }
}

