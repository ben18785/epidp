
#' Estimate the time-varying reproduction number from incidence data
#'
#' @param N number of data points
#' @param C case counts
#' @param w discretised serial interval
#' @param ...
#'
#' @return a stanfit object
#' @export
fit_epifilter <- function(N, C, w, ...) {
  wmax <- length(w)
  standata <- list(N=N,
                   C=C,
                   wmax=wmax,
                   w=w)
  out <- rstan::sampling(stanmodels$epifilter, data = standata, ...)
  return(out)
}

#' Estimate the time-varying reproduction number from incidence data using covariate
#' data
#'
#' @param N number of data points
#' @param C case counts
#' @param w discretised serial interval
#' @param X matrix of covariates (which likely should include a first column of 1s)
#' of dimensions N x (N_covariates + 1)
#' @param ...
#'
#' @return a stanfit object
#' @export
fit_epifilter_covariates <- function(N, C, w, X, ...) {
  wmax <- length(w)
  N_covariates <- dim(X)[2]
  standata <- list(N=N,
                   C=C,
                   w=w,
                   wmax=wmax,
                   N_covariates=N_covariates,
                   X=X)
  out <- rstan::sampling(stanmodels$epifilter, data = standata, ...)
  return(out)
}
