#' Simulate change-point data
#'
#' @param par.values A list with a list with parameter for each cluster in the data vector with, at most, one element.
#' @param d A scalar representing the number of clusters.
#' @param M A scalar representing the number of points available for each observation
#' @param w A scalar representing the minimum number of points in each interval between two change points
#' @param N A scalar representing the number of observations
#' @param as The hyperparameter value for the shape parameter in the inverse-gamma prior for the variance component
#' @param bs The hyperparameter value for the scale parameter in the inverse-gamma prior for the variance component
#' @param emp A boolean value (TRUE/FALSE) for the usage of empirical variance
#' @param seed A scalar representing the seed
#'
#'
#' @return A list composed of a matrix with the data, a vector with the cluster assignments for each observation,
#'    and a vector with the variance for each observation.
#' @export
#'
#' @examples
#' set.seed(1238)
#' sigma2 <- round(extraDistr::rinvgamma(10, 21, 0.05*20),4)
#' param1 <- list(list(2, c(1, 19, 34, 51), c(5, 20, 10)),
#'           list(2, c(1, 15, 32, 51), c(17, 10, 2)),
#'           sigma2)
#' simdata(par.values = param1, d = 2, M = 50, w = 10, N = 10, as = 2, bs = 1000,
#'         emp = FALSE, seed = 1252)
simdata <- function(par.values, d, M, w, N, as, bs, emp, seed){

  # cluster assignment (may not produce balanced clusters)
  set.seed(seed)
  cluster <- sample(c(1:d), N, replace = T)

  lcount <- 1
  while(length(unique(cluster)) == 1){
    set.seed(seed + lcount)
    cluster <- sample(c(1:d), N, replace = T)
    lcount <- lcount + 1
  }

  #fix sigma2
  sigma2 <- par.values[[3]]

  # sample error term
  set.seed(seed)
  error <- MASS::mvrnorm(M, rep(0, N), Sigma = diag(sigma2), empirical = emp)

  y <- matrix(0, M, N)

  # obs in columns and bins in rows
  for(i in 1:d){

    if(par.values[[i]][[1]] > 0){
      Tl <- par.values[[i]][[2]]
      alphal <- par.values[[i]][[3]]
      obs <- which(cluster == i)
      nl <- diff(Tl)
      alpha <- rep(alphal, nl)
      y[,obs] <- alpha + error[,obs]
    } else{
      Tl <- c(1, par.values[[i]][[2]] + 1)
      alphal <- par.values[[i]][[3]]
      obs <- which(cluster == i)
      nl <- diff(Tl)
      alpha <- rep(alphal, nl)
      y[,obs] <- alpha + error[,obs]
    }
  }
  return(list(data = y, cluster_true = cluster, sigma2_true = sigma2))
}

