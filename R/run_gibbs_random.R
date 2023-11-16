# computing posterior estimates
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

