#' Simulates a dataset using `simdata` function and runs the Gibbs sampler algorithm
#'    for simulated scenarios using initial values
#'
#' @param par.values A list with lists with parameters for each cluster.
#'    The first argument in each list is the number of change points, then the positions for the change points,
#'    where $T_1 = 1, T_{last} = M + 1$, and for each interval between change points
#'    you need to specify a value for the constant level.
#' @param d A scalar representing the number of clusters.
#' @param M A scalar representing the number of points available for each observation
#' @param w A scalar representing the minimum number of points in each interval between two change points
#' @param N A scalar representing the number of observations
#' @param as The hyperparameter value for the shape parameter in the inverse-gamma prior for the variance
#'    component
#' @param bs The hyperparameter value for the scale parameter in the inverse-gamma prior for the variance
#'    component
#' @param al The hyperparameter value for the shape parameter in the gamma prior for lambda
#' @param bl The hyperparameter value for the scale parameter in the gamma prior for lambda
#' @param a The hyperparameter value for the shape parameter in the gamma prior for alpha0
#' @param b The hyperparameter value for the scale parameter in the gamma prior for alpha0
#' @param alpha0 A scalar defining the parameter for the Dirichlet process prior
#'    that controls the number of clusters (or its initial values)
#' @param lambda A scalar defining the parameter for the Truncate Poisson distribution
#'    that controls the number of change points (or its initial values)
#' @param maxIter A scalar for the number of iteration to run in the Gibbs sampler
#' @param emp A boolean value (TRUE/FALSE) for the usage of empirical variance. Default is set to `FALSE`
#' @param initial_val A list with the same structure as the `par.values` argument specifying the
#'    initial values for each parameter to estimate
#' @param seed A scalar representing the seed
#'
#'
#'
#' @return A list with the output from the `simdata` function and the estimates for each iteration
#'     of the Gibbs sampler for each parameter
#' @export
#'
run_gibbs_random <- function(M, N, w, d, as = 2, bs = 100, al = 2, bl = 1000, a = 2, b = 1000, alpha0 = 1/100, lambda = 2,
                             maxIter = 10000, par.values, emp, initial_val, seed){

  # maximum number of change-points
  kstar <- ((M - 1)/w) - 1

  # simulate observations accordingly to scenario 1
  data <- simdata(par.values =  par.values, d = d, M = M, w = w, N = N, bs = bs, as = as, emp = emp, seed = seed)

  print(data)

  cluster_miss <- data[[2]]

  if(cluster_miss[1] == 1){
    cluster_miss[1] <- 2
  } else{
    cluster_miss[1] <- 1
  }

  if(cluster_miss[N] == 1){
    cluster_miss[N] <- 2
  } else{
    cluster_miss[N] <- 1
  }

  if(N == 50){

    if(cluster_miss[2] == 1){
      cluster_miss[2] <- 2
    } else{
      cluster_miss[2] <- 1
    }

    if(cluster_miss[10] == 1){
      cluster_miss[10] <- 2
    } else{
      cluster_miss[10] <- 1
    }

    if(cluster_miss[25] == 1){
      cluster_miss[25] <- 2
    } else{
      cluster_miss[25] <- 1
    }
  }

  K <- initial_val[[1]]

  Tl <- initial_val[[2]]

  cluster <- cluster_miss

  alpha <- initial_val[[3]]

  sigma2 <- initial_val[[4]]


  res <- gibbs_alg(alpha0 = alpha0, N = N, w = w, M = M,
                   K = K,
                   Tl = Tl,
                   cluster = cluster,
                   alpha = alpha,
                   sigma2 = sigma2,
                   bs = bs, as = as, al = al, bl = bl, a = a, b = b, kstar = kstar,lambda = lambda,
                   Y = data[[1]], d = d, maxIter = maxIter)


  return(list(data, res))
}

