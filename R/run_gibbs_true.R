# computing posterior estimates
run_gibbs_true <- function(M, N, w, d, as = 2, bs = 100, al = 2, bl = 1000, a = 2, b = 1000, alpha0 = 1/100, lambda = 2,
                           maxIter = 10000, par.values, emp, seed){

  # maximum number of change-points
  kstar <- ((M - 1)/w) - 1

  # simulate observations accordingly to model
  data <- simdata(par.values =  par.values, d = d, M = M, w = w, N = N, bs = bs, as = as, emp = emp, seed = seed)

  print(data)

  K <- c(par.values[[1]][[1]], par.values[[2]][[1]])

  Tl <- list(par.values[[1]][[2]], par.values[[2]][[2]])

  res <- gibbs_alg(alpha0 = alpha0, N = N, w = w, M = M,
                  K = K,
                  Tl = Tl,
                  cluster = data[[2]],
                  alpha = list(par.values[[1]][[3]], par.values[[2]][[3]]),
                  sigma2 = par.values[[3]],
                  bs = bs, as = as, al = al, bl = bl, a = a, b = b, kstar = kstar,
                  lambda = lambda, Y = data[[1]], d = d, maxIter = maxIter)


  return(list(data, res))
}




