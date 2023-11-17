#' Run and saves the Gibbs sampler algorithm results for simulated scenarios when
#'    the initial values are set to be equal to the true parameters
#'
#' @param S A scalar representing the number of simulations to run.
#' @param scenario A string with the name of the folder you want to create to save the results
#' @param casenumber A scalar with the id for the simulation scenario you want to run.
#'     If only one simulation will be performed set this argument to 1.
#' @param par.values A list with lists with parameters for each cluster.
#'    The first argument in each list is the number of change points, then the positions for the change points,
#'    where $T_1 = 1, T_{last} = M + 1$, and for each interval between change points
#'    you need to specify a value for the constant level.
#' @param d A scalar representing the number of clusters.
#' @param M A scalar representing the number of points available for each observation
#' @param w A scalar representing the minimum number of points in each interval between two change points
#' @param N A scalar representing the number of observations
#' @param as The hyperparameter value for the shape parameter in the inverse-gamma prior for the variance component
#' @param bs The hyperparameter value for the scale parameter in the inverse-gamma prior for the variance component
#' @param al The hyperparameter value for the shape parameter in the gamma prior for lambda
#' @param bl The hyperparameter value for the scale parameter in the gamma prior for lambda
#' @param a The hyperparameter value for the shape parameter in the gamma prior for alpha0
#' @param b The hyperparameter value for the scale parameter in the gamma prior for alpha0
#' @param alpha0 A scalar defining the parameter for the Dirichlet process prior
#'    that controls the number of clusters
#' @param lambda A scalar defining the parameter for the Truncate Poisson distribution
#'    that controls the number of change points
#' @param maxIter A scalar for the number of iteration to run in the Gibbs sampler
#' @param startpoint_sampling A scalar defining the step value to use when selecting a systematic sample
#'    after the burn-in.
#' @param emp A boolean value (TRUE/FALSE) for the usage of empirical variance. Default is set to `FALSE`
#' @param sigma2_equal A boolean value (TRUE/FALSE) to determine if the observations should have the same variance.
#'    Default is set to `FALSE`
#' @param seed A scalar representing the seed
#'
#'
#' @return It saves the posterior estimates in the folder provided and returns the processing time.
#' @export
#'
sim_true <- function(S, seed, scenario, casenumber, M, N, w, d, as, bs, al, bl,
                          a, b, alpha0, lambda, maxIter, startpoint_sampling, par.values, sigma2_equal = FALSE, emp = FALSE){
  if(sigma2_equal == TRUE){
    #source("Gibbs_commonsigma2_v2.R")
  } else{
    #source("gibbs_alg.R")
  }

  dir.create(scenario, showWarnings = FALSE)
  dir.create(paste0(scenario, "/case", casenumber), showWarnings = FALSE)


  cls <- parallel::makeForkCluster(parallel::detectCores())
  outt <- parallel::parSapply(cls, 1:S, save_sim_true, seed = seed, scenario = scenario,
                    casenumber = casenumber, M = M, N = N, w = w,
                    d = d, as = as, bs = bs, al = al, bl = bl, a = a, b = b, alpha0 = alpha0, lambda = lambda,
                    maxIter = maxIter, startpoint_sampling = startpoint_sampling, par.values = par.values, emp = emp)
  parallel::stopCluster(cls)
  return(outt)
}
