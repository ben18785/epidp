

#' Generate a Discretized Serial Interval Distribution
#'
#' This function generates a discretized serial interval distribution based on
#' the gamma distribution characterized by the specified mean and standard deviation.
#'
#' @param nt An integer specifying the number of time points (length of the distribution).
#' @param mean_si A numeric value representing the mean of the serial interval distribution.
#' @param sd_si A numeric value representing the standard deviation of the serial interval distribution.
#'
#' @return A numeric vector of length `nt` representing the discretized serial interval distribution.
#'
#' @details
#' The serial interval distribution is discretized using a gamma distribution
#' with shape and scale parameters derived from the specified mean and standard deviation.
#' The discretization is performed by computing the difference in the cumulative distribution
#' function (CDF) of the gamma distribution at successive integer time points.
#'
#' @examples
#' \dontrun{
#' generate_vector_serial(10, 5, 2)
#' }
#' @export
generate_vector_serial <- function(nt, mean_si, sd_si) {

  if (!is.numeric(nt) || nt <= 0 || nt != as.integer(nt)) {
    stop("Parameter 'nt' should be a positive integer.")
  }
  if (!is.numeric(mean_si) || mean_si <= 0) {
    stop("Parameter 'mean_si' should be a positive numeric value.")
  }
  if (!is.numeric(sd_si) || sd_si <= 0) {
    stop("Parameter 'sd_si' should be a positive numeric value.")
  }

  # Shape-scale parameters of gamma serial interval
  shape <- mean_si^2 / sd_si^2
  scale <- sd_si^2 / mean_si

  # Time points and cumulative distribution values
  tdist <- 0:nt
  cdf_vals <- pgamma(tdist, shape = shape, scale = scale)

  # Compute differences to get the discretized distribution
  w_dist <- diff(cdf_vals)

  w_dist
}

#' Simulate a Renewal Epidemic Model
#'
#' This function simulates an epidemic using a renewal model based on the
#' specified time-varying reproduction number (Rt), the serial interval distribution
#' characterized by a gamma distribution, and an initial number of cases.
#'
#' @param rt_fun A function that takes a vector of time points and returns a vector
#'   of reproduction numbers (Rt) corresponding to those time points.
#' @param nt An integer specifying the number of time points (duration of the epidemic simulation).
#' @param mean_si A numeric value representing the mean of the serial interval distribution.
#' @param sd_si A numeric value representing the standard deviation of the serial interval distribution.
#' @param i_0 An integer specifying the initial number of cases at time point 1.
#' @param X An optional matrix of covariates with `nt` rows.
#'
#' @return A list containing the following components:
#' \describe{
#'   \item{i_t}{A numeric vector of the incidence of new cases at each time point.}
#'   \item{lambda_t}{A numeric vector of the total infectiousness at each time point.}
#'   \item{w_dist}{A numeric vector of the discretized serial interval distribution.}
#'   \item{r_t}{A numeric vector of the reproduction number (Rt) at each time point.}
#'   \item{t}{A numeric vector of the time points from 1 to nt.}
#' }
#'
#' @details
#' The renewal model is a common approach for simulating the spread of infectious diseases.
#' The total infectiousness at each time point is computed as the convolution of past incidence
#' with the serial interval distribution. The incidence at each time point is then simulated
#' as a Poisson random variable with a mean equal to the product of the total infectiousness and
#' the time-varying reproduction number.
#'
#' @examples
#' \dontrun{
#' rt_fun <- function(t) { 1.5 * exp(-0.05 * t) }
#' simulate_renewal_epidemic(rt_fun, 100, 5, 2, 10)
#' }
#' @export
simulate_renewal_epidemic <- function(Rt_fun, nt, mean_si, sd_si, i_0, X=NULL){

  # Input validation
  if (!is.numeric(nt) || nt <= 0 || nt != as.integer(nt)) {
    stop("Parameter 'nt' should be a positive integer.")
  }
  if (!is.numeric(mean_si) || mean_si <= 0) {
    stop("Parameter 'mean_si' should be a positive numeric value.")
  }
  if (!is.numeric(sd_si) || sd_si <= 0) {
    stop("Parameter 'sd_si' should be a positive numeric value.")
  }
  if (!is.numeric(i_0) || i_0 <= 0 || i_0 != as.integer(i_0)) {
    stop("Parameter 'i_0' should be a positive integer.")
  }
  if (!is.function(Rt_fun)) {
    stop("Parameter 'Rt_fun' should be a function.")
  }
  if (!is.null(X) && (!is.matrix(X) || nrow(X) != nt)) {
    stop("Parameter 'X' should be a matrix with 'nt' rows.")
  }

  # Time series and Rt
  t = 1:nt
  Rt <- vector(length = nt)
  for(i in seq_along(Rt)) {
    if(is.null(X))
      Rt[i] = Rt_fun(t[i])
    else
      Rt[i] <- Rt_fun(t[i], X[i, ])
  }

  # Total infectiousness and incidence with initial imports
  Lt = rep(0, nt); It = Lt; It[1] = i_0; Lt[1] = It[1]

  w_dist <- generate_vector_serial(nt, mean_si, sd_si)

  # Simulate from standard renewal model
  for(i in 2:nt){
    # Total infectiousness is a convolution
    Lt[i] = sum(It[seq(i-1, 1, -1)] * w_dist[1:(i-1)])
    # Poisson renewal model
    It[i] = rpois(1, Lt[i] * Rt[i])
  }

  data.frame(t=t, i_t=It, R_t=Rt, lambda_t=Lt, w_dist=w_dist)
}
