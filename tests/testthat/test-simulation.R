test_that("Output is a numeric vector of correct length", {
  nt <- 10
  mean_si <- 5
  sd_si <- 2
  result <- generate_vector_serial(nt, mean_si, sd_si)
  expect_type(result, "double")
  expect_length(result, nt)
})

# Test 2: Verify the function's output for known inputs
test_that("Output is correct for known inputs", {
  nt <- 1
  mean_si <- 1
  sd_si <- 1
  expected_result <- c(pgamma(1, shape = 1, scale = 1) - pgamma(0, shape = 1, scale = 1))
  result <- generate_vector_serial(nt, mean_si, sd_si)
  expect_equal(result, expected_result)
})

# Test 4: Test the function with very small values for mean_si and sd_si
test_that("Output is correct for small mean_si and sd_si", {
  nt <- 10
  mean_si <- 0.1
  sd_si <- 0.1
  result <- generate_vector_serial(nt, mean_si, sd_si)
  expect_true(all(result >= 0))
  expect_equal(sum(result), 1)
})

test_that("Function handles nt = 0", {
  expect_error(generate_vector_serial(0, 5, 2), "Parameter 'nt' should be a positive integer.")
})

test_that("Function handles non-integer nt", {
  expect_error(generate_vector_serial(1.5, 5, 2), "Parameter 'nt' should be a positive integer.")
})

# Test 5: Test the function with negative nt (invalid input)
test_that("Function handles negative nt", {
  expect_error(generate_vector_serial(-5, 5, 2), "Parameter 'nt' should be a positive integer.")
})

# Test 6: Test the function with non-numeric nt (invalid input)
test_that("Function handles non-numeric nt", {
  expect_error(generate_vector_serial("ten", 5, 2), "Parameter 'nt' should be a positive integer.")
})

# Test 7: Test the function with mean_si <= 0 (invalid input)
test_that("Function handles mean_si <= 0", {
  expect_error(generate_vector_serial(10, 0, 2), "Parameter 'mean_si' should be a positive numeric value.")
  expect_error(generate_vector_serial(10, -1, 2), "Parameter 'mean_si' should be a positive numeric value.")
})

# Test 8: Test the function with non-numeric mean_si (invalid input)
test_that("Function handles non-numeric mean_si", {
  expect_error(generate_vector_serial(10, "five", 2), "Parameter 'mean_si' should be a positive numeric value.")
})

# Test 9: Test the function with sd_si <= 0 (invalid input)
test_that("Function handles sd_si <= 0", {
  expect_error(generate_vector_serial(10, 5, 0), "Parameter 'sd_si' should be a positive numeric value.")
  expect_error(generate_vector_serial(10, 5, -2), "Parameter 'sd_si' should be a positive numeric value.")
})

# Test 10: Test the function with non-numeric sd_si (invalid input)
test_that("Function handles non-numeric sd_si", {
  expect_error(generate_vector_serial(10, 5, "two"), "Parameter 'sd_si' should be a positive numeric value.")
})

test_that("simulate_renewal_epidemic returns a data frame", {
  rt_fun <- function(t) {
    1.5 * exp(-0.05 * t)
  }
  result <- simulate_renewal_epidemic(rt_fun, 100, 5, 2, 10)
  expect_type(result, "list")
  expect_s3_class(result, "data.frame")
})

test_that("simulate_renewal_epidemic returns correct dimensions", {
  rt_fun <- function(t) {
    1.5 * exp(-0.05 * t)
  }
  nt <- 100
  result <- simulate_renewal_epidemic(rt_fun, nt, 5, 2, 10)
  expect_equal(nrow(result), nt)
  expect_equal(ncol(result), 5)
})

test_that("simulate_renewal_epidemic returns a data frame", {
  rt_fun <- function(t) {
    1.5 * exp(-0.05 * t)
  }
  result <- simulate_renewal_epidemic(rt_fun, 100, 5, 2, 10)
  expect_type(result, "list")
  expect_s3_class(result, "data.frame")
})

test_that("simulate_renewal_epidemic returns correct dimensions", {
  rt_fun <- function(t) {
    1.5 * exp(-0.05 * t)
  }
  nt <- 100
  result <- simulate_renewal_epidemic(rt_fun, nt, 5, 2, 10)
  expect_equal(nrow(result), nt)
  expect_equal(ncol(result), 5)
})

test_that("simulate_renewal_epidemic handles invalid nt", {
  rt_fun <- function(t) {
    1.5 * exp(-0.05 * t)
  }
  expect_error(simulate_renewal_epidemic(rt_fun, -10, 5, 2, 10), "Parameter 'nt' should be a positive integer.")
  expect_error(simulate_renewal_epidemic(rt_fun, 0, 5, 2, 10), "Parameter 'nt' should be a positive integer.")
  expect_error(simulate_renewal_epidemic(rt_fun, 10.5, 5, 2, 10), "Parameter 'nt' should be a positive integer.")
})

test_that("simulate_renewal_epidemic handles invalid mean_si", {
  rt_fun <- function(t) {
    1.5 * exp(-0.05 * t)
  }
  expect_error(simulate_renewal_epidemic(rt_fun, 100, -5, 2, 10), "Parameter 'mean_si' should be a positive numeric value.")
  expect_error(simulate_renewal_epidemic(rt_fun, 100, 0, 2, 10), "Parameter 'mean_si' should be a positive numeric value.")
  expect_error(simulate_renewal_epidemic(rt_fun, 100, "five", 2, 10), "Parameter 'mean_si' should be a positive numeric value.")
})

test_that("simulate_renewal_epidemic handles invalid sd_si", {
  rt_fun <- function(t) {
    1.5 * exp(-0.05 * t)
  }
  expect_error(simulate_renewal_epidemic(rt_fun, 100, 5, -2, 10), "Parameter 'sd_si' should be a positive numeric value.")
  expect_error(simulate_renewal_epidemic(rt_fun, 100, 5, 0, 10), "Parameter 'sd_si' should be a positive numeric value.")
  expect_error(simulate_renewal_epidemic(rt_fun, 100, 5, "two", 10), "Parameter 'sd_si' should be a positive numeric value.")
})

test_that("simulate_renewal_epidemic handles invalid i_0", {
  rt_fun <- function(t) {
    1.5 * exp(-0.05 * t)
  }
  expect_error(simulate_renewal_epidemic(rt_fun, 100, 5, 2, -10), "Parameter 'i_0' should be a positive integer.")
  expect_error(simulate_renewal_epidemic(rt_fun, 100, 5, 2, 0), "Parameter 'i_0' should be a positive integer.")
  expect_error(simulate_renewal_epidemic(rt_fun, 100, 5, 2, 10.5), "Parameter 'i_0' should be a positive integer.")
})

test_that("simulate_renewal_epidemic handles invalid Rt_fun", {
  expect_error(simulate_renewal_epidemic(NULL, 100, 5, 2, 10), "Parameter 'Rt_fun' should be a function.")
  expect_error(simulate_renewal_epidemic(5, 100, 5, 2, 10), "Parameter 'Rt_fun' should be a function.")
})

test_that("simulate_renewal_epidemic returns consistent results for fixed Rt_fun", {
  rt_fun <- function(t) {
    rep(2, length(t))
  }
  set.seed(123)
  result1 <- simulate_renewal_epidemic(rt_fun, 10, 2, 1, 1)
  set.seed(123)
  result2 <- simulate_renewal_epidemic(rt_fun, 10, 2, 1, 1)
  expect_equal(result1, result2)
})

test_that("simulate_renewal_epidemic produces non-negative incidence", {
  rt_fun <- function(t) {
    1.5 * exp(-0.05 * t)
  }
  result <- simulate_renewal_epidemic(rt_fun, 100, 5, 2, 10)
  expect_true(all(result$i_t >= 0))
})
