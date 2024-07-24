simRenewalEpidemic <- function(RtFun, nt, mean_si, sd_si, I0){

  # Time series and Rt
  t = 1:nt; Rt = RtFun(t)

  # Shape-scale parameters of gamma serial interval
  pms = c(0,0); pms[1] = mean_si^2/sd_si^2; pms[2] = sd_si^2/mean_si

  # Discretise serial interval distribution
  tdist = c(0, t); wdist = rep(0, nt)
  for (i in 1:nt){
    wdist[i] = pgamma(tdist[i+1], shape = pms[1], scale = pms[2]) -
      pgamma(tdist[i], shape = pms[1], scale = pms[2])
  }

  # Total infectiousness and incidence with initial imports
  Lt = rep(0, nt); It = Lt; It[1] = I0; Lt[1] = It[1]

  # Simulate from standard renewal model
  for(i in 2:nt){
    # Total infectiousness is a convolution
    Lt[i] = sum(It[seq(i-1, 1, -1)]*wdist[1:(i-1)])
    # Poisson renewal model
    It[i] = rpois(1, Lt[i]*Rt[i])
  }
  # Main outputs
  out = list(It, Lt, wdist, Rt, t)
  names(out) = c('It', 'Lt', 'wdist', 'Rt', 't')
  simRenewalEpidemic = out
}
