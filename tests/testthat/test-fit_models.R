test_that("fit_epifilter runs", {
  rt_fun <- function(t) {
    if (t <= 60) {
      R <- 2
    } else if (t <= 90) {
      R <- 0.5
    } else {
      R <- 1
    }
    R
  }

  # simulation parameters
  nt <- 200
  mean_si <- 6.5
  sd_si <- 4.03
  i_0 <- 10

  # data frame of outputs
  epidemic_df <- simulate_renewal_epidemic(rt_fun, nt, mean_si, sd_si, i_0)

  # test optimisation runs
  fit <- fit_epifilter(
    N = length(epidemic_df$i_t),
    C = epidemic_df$i_t,
    w = epidemic_df$w_dist,
    is_sampling = FALSE,
    as_vector = FALSE
  )
  expect_true(is.numeric(fit$par$sigma))

  # test sampling runs
  fit <- suppressWarnings({ # here because Stan outputs sampling warnings
    fit_epifilter(
      N = length(epidemic_df$i_t),
      C = epidemic_df$i_t,
      w = epidemic_df$w_dist,
      is_sampling = TRUE,
      iter = 10,
      chains = 1,
      refresh = 0
    )
  })
  expect_s4_class(fit, "stanfit")
})
