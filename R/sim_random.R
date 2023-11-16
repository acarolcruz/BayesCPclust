sim_random <- function(S, seed, scenario, casenumber, M, N, w, d, as, bs, al, bl, a, b, alpha0, lambda, maxIter,
                          startpoint_sampling, par.values, sigma2_equal = FALSE, emp = FALSE, initial_val){
  if(sigma2_equal == TRUE){
    source("Gibbs_commonsigma2_v2.R")
  } else{
    source("Gibbs_new_v6_all_clusters.R")
  }

  dir.create(scenario, showWarnings = FALSE)
  dir.create(paste0(scenario, "/case", casenumber), showWarnings = FALSE)

  cls <- parellel::makeForkCluster(parellel::detectCores())
  outt <- parSapply(cls, 1:S, save_sim_random, seed = seed, scenario = scenario,
                    casenumber = casenumber, M = M, N = N, w = w,
                    d = d, as = as, bs = bs, al = al, bl = bl, a = a, b = b, alpha0 = alpha0, lambda = lambda,
                    maxIter = maxIter, startpoint_sampling = startpoint_sampling, par.values = par.values,
                    emp = emp, initial_val = initial_val)
  parellel::stopCluster(cls)


  return(outt)
}
