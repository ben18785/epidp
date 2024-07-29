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
simulate_renewal_epidemic <- function(rt_fun, nt, mean_si, sd_si, i_0){

  # Time series and Rt
  t = 1:nt; Rt = rt_fun(t)

  # Shape-scale parameters of gamma serial interval
  pms = c(0,0); pms[1] = mean_si^2/sd_si^2; pms[2] = sd_si^2/mean_si

  # Discretise serial interval distribution
  tdist = c(0, t); wdist = rep(0, nt)
  for (i in 1:nt){
    wdist[i] = pgamma(tdist[i+1], shape = pms[1], scale = pms[2]) -
      pgamma(tdist[i], shape = pms[1], scale = pms[2])
  }

  # Total infectiousness and incidence with initial imports
  Lt = rep(0, nt); It = Lt; It[1] = i_0; Lt[1] = It[1]

  # Simulate from standard renewal model
  for(i in 2:nt){
    # Total infectiousness is a convolution
    Lt[i] = sum(It[seq(i-1, 1, -1)]*wdist[1:(i-1)])
    # Poisson renewal model
    It[i] = rpois(1, Lt[i]*Rt[i])
  }
  
  list(i_t=i_t, lambda_t=lambda_t, w_dist=w_dist, r_t=r_t, t=t)
}
