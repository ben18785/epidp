
#' Estimate the time-varying reproduction number from incidence data
#'
#' @param N number of data points
#' @param C case counts
#' @param w discretised serial interval
#' @param is_sampling Boolean indicating if sampling is to be used (if TRUE) or
#' optimisation (if FALSE). Defaults to TRUE
#' @param ... Arguments passed to [rstan::sampling()] (if `is_sampling=TRUE`) or
#' [rstan::optimizing()] (if `is_sampling=FALSE`)
#'
#' @return a stanfit object
#' @export
fit_epifilter <- function(N, C, w, is_sampling=TRUE, ...) {
  wmax <- length(w)
  standata <- list(N=N,
                   C=C,
                   wmax=wmax,
                   w=w)
  if(is_sampling)
    out <- rstan::sampling(stanmodels$epifilter, data = standata, ...)
  else
    out <- rstan::optimizing(stanmodels$epifilter, data = standata, ...)
  return(out)
}

#' Estimate the time-varying reproduction number from incidence data using covariate
#' data
#'
#' @inheritParams fit_epifilter
#' @param X matrix of covariates (which likely should include a first column of 1s)
#' of dimensions N x (N_covariates + 1)
#'
#' @return a stanfit object
#' @export
fit_epifilter_covariates <- function(N, C, w, X, is_sampling=TRUE, ...) {
  wmax <- length(w)
  N_covariates <- dim(X)[2]
  standata <- list(N=N,
                   C=C,
                   w=w,
                   wmax=wmax,
                   N_covariates=N_covariates,
                   X=X)
  if(is_sampling)
    out <- rstan::sampling(stanmodels$epifilter_covariates, data = standata, ...)
  else
    out <- rstan::optimizing(stanmodels$epifilter_covariates, data = standata, ...)
  return(out)
}
