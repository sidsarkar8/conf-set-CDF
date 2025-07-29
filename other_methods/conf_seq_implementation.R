# --- confseq.misc functionality ---

# Return the running intersection of lower and upper confidence bounds.
# l: numeric vector of lower bounds
# u: numeric vector of upper bounds
get_running_intersection <- function(l, u) {
  list(lower = cummax(l), upper = cummin(u))
}

# Get a sequence of confidence intervals.
# x: numeric vector (RealArray) of observations between 0 and 1.
# ci_fn: a function that takes a RealArray of observations and returns a 
#        two-element numeric vector c(l, u) for the lower and upper CI.
# times: vector of positive integers at which to compute the CI.
# parallel: if TRUE, computation is parallelized using mclapply.
get_ci_seq <- function(x, ci_fn, times, parallel = FALSE) {
  x <- as.numeric(x)
  l <- rep(0.0, length(times))
  u <- rep(1.0, length(times))
  
  if (parallel) {
    if (!requireNamespace("parallel", quietly = TRUE)) {
      stop("The 'parallel' package is required for parallel execution.")
    }
    n_cores <- parallel::detectCores()
    cat("Using", n_cores, "cores\n")
    results <- parallel::mclapply(times, function(time) {
      ci_fn(x[1:time])
    }, mc.cores = n_cores)
    results <- do.call(rbind, results)
    l <- results[, 1]
    u <- results[, 2]
  } else {
    for (i in seq_along(times)) {
      time <- times[i]
      ci <- ci_fn(x[1:time])
      l[i] <- ci[1]
      u[i] <- ci[2]
    }
  }
  
  return(list(lower = l, upper = u))
}

# --- End confseq.misc functionality ---

# --- conf_seq.betting_strategies functionality ---

# lambda_predmix_eb: Predictable mixture lambda values ("bets")
# Parameters:
#   x: RealArray (numeric vector) of observations (bounded between 0 and 1)
#   truncation: positive real (or Inf) at which to truncate lambda
#   alpha: significance level in (0, 1)
#   fixed_n: positive integer or NULL; if not NULL, lambda is optimized for sample size fixed_n
#   prior_mean: prior mean for regularization (default 0.5)
#   prior_variance: prior variance for regularization (default 0.25)
#   fake_obs: number of "fake observations" to add (default 1)
#   scale: multiplier for final lambda (default 1)
lambda_predmix_eb <- function(x, truncation = Inf, alpha = 0.05, fixed_n = NULL,
                              prior_mean = 1/2, prior_variance = 1/4, fake_obs = 1, scale = 1) {
  t <- seq_along(x)
  # Regularized running mean (capped at 1)
  mu_hat_t <- pmin((fake_obs * prior_mean + cumsum(x)) / (t + fake_obs), 1)
  # Compute running variance based on squared deviation from mu_hat_t
  sigma2_t <- (fake_obs * prior_variance + cumsum((x - mu_hat_t)^2)) / (t + fake_obs)
  # Regularized "lagged" variance: use prior_variance for time 1 and then previous sigma2_t values
  sigma2_tminus1 <- c(prior_variance, sigma2_t[1:(length(x) - 1)])
  
  if (is.null(fixed_n)) {
    lambdas <- sqrt(2 * log(1 / alpha) / (t * log(1 + t) * sigma2_tminus1))
  } else {
    lambdas <- sqrt(2 * log(1 / alpha) / (fixed_n * sigma2_tminus1))
  }
  
  # Replace any NaN with 0
  lambdas[is.nan(lambdas)] <- 0
  
  # Truncate lambda values elementwise
  lambdas <- pmin(truncation, lambdas)
  
  return(lambdas * scale)
}

# --- End conf_seq.betting_strategies functionality ---

# --- Confidence Sequence Functions ---

# Predictable mixture lower confidence sequence.
# x: RealArray of observations.
# v: RealArray (variance increments).
# lambdas_fn: function mapping a RealArray to lambda values (RealArray).
# psi_fn: function mapping a RealArray of lambda values to another RealArray.
# alpha: significance level.
# running_intersection: if TRUE, returns the running intersection (cumulative max).
# N: population size if sampling without replacement.
predmix_lower_cs <- function(x, v, lambdas_fn, psi_fn, alpha = 0.05,
                             running_intersection = FALSE, N = NULL) {
  x <- as.numeric(x)
  t <- seq_along(x)
  
  S_t <- cumsum(x)
  S_tminus1 <- c(0, S_t[-length(x)])
  
  if (!is.null(N)) {
    Zstar <- S_tminus1 / (N - t + 1)
    Wstar <- (t - 1) / (N - t + 1)
  } else {
    Zstar <- 0
    Wstar <- 0
  }
  
  lambdas <- lambdas_fn(x)
  psi <- psi_fn(lambdas)
  
  margin <- (log(1 / alpha) + cumsum(v * psi)) / cumsum(lambdas * (1 + Wstar))
  weighted_mu_hat_t <- cumsum(lambdas * (x + Zstar)) / cumsum(lambdas * (1 + Wstar))
  
  l <- weighted_mu_hat_t - margin
  l <- pmax(l, 0)
  
  if (running_intersection) {
    l <- cummax(l)
  }
  
  return(l)
}

# Predictable mixture empirical Bernstein lower confidence sequence.
# x: RealArray of observations (assumed nonnegative).
# alpha: significance level.
# truncation: level at which to truncate lambda.
# running_intersection: if TRUE, returns running intersection.
# N: population size if sampling without replacement.
# fixed_n: fixed time for optimizing the bound.
# truncate_mean: if TRUE, cap the running mean at 1.
predmix_empbern_lower_cs <- function(x, alpha = 0.05, truncation = 1/2,
                                     running_intersection = FALSE, N = NULL,
                                     fixed_n = NULL, truncate_mean = FALSE) {
  t <- seq_along(x)
  if (truncate_mean) {
    mu_hat_t <- pmin(cumsum(x) / t, 1)
  } else {
    mu_hat_t <- cumsum(x) / t
  }
  mu_hat_tminus1 <- c(0, mu_hat_t[-length(x)])
  v <- (x - mu_hat_tminus1)^2
  
  # Use the lambda_predmix_eb defined above.
  lambdas_fn <- function(y) {
    lambda_predmix_eb(y, truncation = truncation, alpha = alpha, fixed_n = fixed_n)
  }
  psi_fn <- function(lambdas) {
    -log(1 - lambdas) - lambdas
  }
  
  return(predmix_lower_cs(x, v, lambdas_fn, psi_fn, alpha, running_intersection, N))
}

# Predictable mixture Hoeffding lower confidence sequence.
# x: RealArray of observations.
# alpha: significance level.
# truncation: truncation level for lambda.
# running_intersection: if TRUE, return running intersection.
# N: population size if sampling without replacement.
# fixed_n: fixed time for optimizing the bound.
predmix_hoeffding_lower_cs <- function(x, alpha = 0.05, truncation = 1,
                                       running_intersection = FALSE, N = NULL,
                                       fixed_n = NULL) {
  t <- seq_along(x)
  if (!is.null(fixed_n)) {
    lambdas_fn <- function(y) {
      rep(sqrt(8 * log(1 / alpha) / fixed_n), length(x))
    }
  } else {
    lambdas_fn <- function(y) {
      pmin(sqrt(8 * log(1 / alpha) / (t * log(1 + t))), truncation)
    }
  }
  psi_fn <- function(lambdas) {
    lambdas^2 / 8
  }
  
  v <- rep(1, length(x))
  return(predmix_lower_cs(x, v, lambdas_fn, psi_fn, alpha, running_intersection, N))
}

# Predictable mixture Hoeffding confidence sequence (two-sided).
# x: RealArray of observations in [0, 1].
# Returns a list with elements 'lower' and 'upper'.
predmix_hoeffding_cs <- function(x, alpha = 0.05, truncation = 1,
                                 running_intersection = FALSE, N = NULL,
                                 fixed_n = NULL) {
  lower_cs <- predmix_hoeffding_lower_cs(x, alpha = alpha / 2,
                                         truncation = truncation,
                                         running_intersection = running_intersection,
                                         N = N, fixed_n = fixed_n)
  upper_cs <- 1 - predmix_hoeffding_lower_cs(1 - x, alpha = alpha / 2,
                                             truncation = truncation,
                                             running_intersection = running_intersection,
                                             N = N, fixed_n = fixed_n)
  return(list(lower = lower_cs, upper = upper_cs))
}

# Predictable mixture empirical Bernstein confidence sequence (two-sided).
# x: RealArray of observations in [0, 1].
predmix_empbern_twosided_cs <- function(x, alpha = 0.05, truncation = 1/2,
                                        running_intersection = FALSE, N = NULL,
                                        fixed_n = NULL) {
  l <- predmix_empbern_lower_cs(x = x, alpha = alpha / 2,
                                truncation = truncation,
                                running_intersection = running_intersection,
                                N = N, fixed_n = fixed_n)
  u <- 1 - predmix_empbern_lower_cs(x = 1 - x, alpha = alpha / 2,
                                    truncation = truncation,
                                    running_intersection = running_intersection,
                                    N = N, fixed_n = fixed_n)
  return(list(lower = l, upper = u))
}

# Deprecated: Predictable mixture empirical Bernstein confidence sequence.
predmix_empbern_cs <- function(x, alpha = 0.05, truncation = 1/2,
                               running_intersection = FALSE, N = NULL,
                               fixed_n = NULL) {
  warning("predmix_empbern_cs is deprecated. Please use predmix_empbern_twosided_cs instead.",
          call. = FALSE)
  return(predmix_empbern_twosided_cs(x, alpha = alpha, truncation = truncation,
                                     running_intersection = running_intersection,
                                     N = N, fixed_n = fixed_n))
}

# Predictable mixture Hoeffding confidence interval (final CI at time n).
# x: RealArray of observations.
predmix_hoeffding_ci <- function(x, alpha = 0.05, N = NULL, running_intersection = TRUE) {
  cs <- predmix_hoeffding_cs(x, alpha = alpha, truncation = Inf,
                             running_intersection = running_intersection,
                             N = N, fixed_n = length(x))
  lower_cs <- cs$lower
  upper_cs <- cs$upper
  return(c(lower = tail(lower_cs, 1), upper = tail(upper_cs, 1)))
}

# Sequence of Hoeffding CIs computed at specified times.
# x: RealArray of observations.
# times: vector of positive integers indicating time points.
predmix_hoeffding_ci_seq <- function(x, times, alpha = 0.05, N = NULL,
                                     running_intersection = TRUE, parallel = FALSE) {
  ci_fn <- function(x_sub) {
    predmix_hoeffding_ci(x_sub, alpha = alpha, N = N, running_intersection = running_intersection)
  }
  return(get_ci_seq(x, ci_fn, times, parallel))
}

# Predictable mixture empirical Bernstein confidence interval (final CI at time n).
# x: RealArray of observations.
predmix_empbern_ci <- function(x, alpha = 0.05, truncation = 1/2, N = NULL,
                               running_intersection = TRUE) {
  cs <- predmix_empbern_cs(x, alpha = alpha, truncation = truncation,
                           N = N, fixed_n = length(x))
  lower_cs <- cs$lower
  upper_cs <- cs$upper
  return(c(lower = tail(lower_cs, 1), upper = tail(upper_cs, 1)))
}

# Sequence of empirical Bernstein CIs computed at specified times.
# x: RealArray of observations.
predmix_empbern_ci_seq <- function(x, times, alpha = 0.05, truncation = 1/2,
                                   N = NULL, running_intersection = TRUE,
                                   parallel = FALSE) {
  ci_fn <- function(x_sub) {
    predmix_empbern_ci(x_sub, alpha = alpha, truncation = truncation,
                       N = N, running_intersect on = running_intersection)
  }
  return(get_ci_seq(x, ci_fn, times, parallel))
}

# --- Example usage ---

# Generate example data and compute a two-sided Hoeffding confidence sequence.
set.seed(123)
x <- runif(100)  # example data in [0, 1]
cs <- predmix_hoeffding_cs(x, alpha = 0.05)
plot(cs$lower, type = "l", col = "blue", ylim = c(0, 1),
     ylab = "Confidence Bound", xlab = "Time")
lines(cs$upper, col = "red")
legend("topright", legend = c("Lower Bound", "Upper Bound"),
       col = c("blue", "red"), lty = 1)
