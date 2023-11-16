#posterior_functions

pk <- function(k, kstar, lambda){
  (lambda^k/factorial(k)) / sum(sapply(0:kstar,FUN=function(l){lambda^l/factorial(l)}))
}


qn0_mk <- function(w, m0, bs, as, M, km, lambda, mk, Yn, kstar){

  Tl <- c()
  if(km == 0){
    Tl <- M
    Xn <- rep(1, Tl)
  } else{
    for(i in 2:(km+1)){
      Tl[1] <- 1
      Tl[i] <- mk[i-1] + Tl[i-1] + w
      #print(i)
    }
    Tl <- append(Tl, M + 1)
    Xn <- apply(diag(km + 1), 2, function(x){rep(x, diff(Tl))})

  }

  Vn <- t(Xn)%*%Xn
  B <- t(Yn)%*%Yn - t(Yn)%*%Xn%*%solve(Vn)%*%t(Xn)%*%Yn

  Hmk <- -((M-(km+1))/2)*log(2*pi) -(1/2)*log(det(Vn)) - as*log(bs) - (((M-(km+1))/2) + as)*log((B/2) + (1/bs)) -
    lgamma(as) + lgamma((((M-(km+1))/2) + as))

  res1 <- exp(Hmk)*exp((lfactorial(m0) - sum(lfactorial(mk))))*((1/(km+1))^m0)*pk(km, kstar, lambda)

  return(res1)
}


qn0 <- function(alpha0, w, N, M, bs, as, kstar, lambda, Yn){

  res2 <- c()
  for(km in 0:kstar){
    m0 <- M - 1 - (km+1)*w
    mk <- permuteGeneral(0:m0, km + 1,
                         constraintFun = "sum",
                         comparisonFun = "==", limitConstraints = m0)


    res2 <-c(res2, sum(apply(mk, 1,
                             FUN = function(x){qn0_mk(w=w, m0=m0, bs=bs, as=as, M=M, km=km,
                                                      lambda=lambda, mk=x, Yn=Yn, kstar=kstar)})))
  }
  return(alpha0*sum(res2))
}


qnj <- function(N, M, as, bs, Yn, alpha, cluster, Tl, K){
  res3 <- c()
  for(clusteri in 1:length(alpha)){
    alphan <- alpha[[clusteri]]
    kj <- K[clusteri]
    Tlj <- Tl[[clusteri]]
    if(kj == 0){
      Xn <- rep(1, Tlj)

      Nj <- sum(cluster == clusteri)

      Hj <- lgamma(M/2 + as) - (M/2 + as)*log(((t((Yn - Xn*alphan))%*%(Yn - Xn*alphan))/2) + (1/bs))
      -(M/2)*log(2*pi) - lgamma(as) - as*log(bs)
    } else{
      Xn <- apply(diag(kj+1), 2, function(x){rep(x, diff(Tlj))})

      Nj <- sum(cluster == clusteri)

      Hj <- lgamma(M/2 + as) - (M/2 + as)*log(((t(Yn - Xn%*%alphan)%*%(Yn - Xn%*%alphan))/2) + (1/bs))
      -(M/2)*log(2*pi) - lgamma(as) - as*log(bs)
    }

    res3 <- c(res3, Nj*exp(Hj))

  }
  return(res3)
}



# update sigma^2n
possigma2n <- function(as, bs, M, Yn, k, Tln, alphan){


  if(k == 0){
    Tln <- M
    Xn <- as.matrix(rep(1, Tln))
  } else{
    Xn <- apply(diag(k+1),2,function(x){rep(x,diff(Tln))})


  }

  an <- (M/2) + as
  bn <- (((t(Yn - Xn%*%alphan)%*%(Yn - Xn%*%alphan))/2) + (1/bs))#rate parameter : exp((1/B)*x) bn = 1/B

  return(rinvgamma(1, alpha = an, beta = bn))
}

postK_mk <- function(k, m0, w, M, Yn, sigma2n, cellsn, mk, Cr){
  Tl <- c()
  if(k == 0){
    Tl <- M
    Xn <- as.matrix(rep(1,Tl))

    Vn <- t(Xn)%*%Xn

    A <- sum(sapply(1:length(cellsn), function(x){(sigma2n[x]^(-1))*t(Xn)%*%Yn[,x]}))

  }else{
    for(i in 2:(k+1)){
      Tl[1] <- 1
      Tl[i] <- mk[i-1] + Tl[i-1] + w
    }
    Tl <- append(Tl, M + 1)
    Xn <- apply(diag(k + 1), 2, function(x){rep(x,diff(Tl))})


    Vn <- t(Xn)%*%Xn

    A <- as.matrix(rowSums(sapply(1:length(cellsn), function(x){(sigma2n[x]^(-1))*t(Xn)%*%Yn[,x]})))

  }



  H <- -(M/2)*sum(log(sigma2n)) - ((M*Cr - k - 1)/2)*log(2*pi) - ((k+1)/2)*log(sum(sigma2n^(-1))) - (1/2)*log(det(Vn)) - (1/2)*sum(sapply(1:length(cellsn), function(x){(sigma2n[x]^(-1))*t(Yn[,x])%*%Yn[,x]})) + (1/2)*t(A)%*%solve((sum(sigma2n^(-1)))*t(Xn)%*%Xn)%*%A



  return(H)
}


logsumexp <- function(x, min_x = Inf){
  min_x <- min(min_x, min(x))
  const <- min_x / (-740) #-740
  return(list(x = x/const, min_x = min_x))
}


postK <- function(kstar, w, M, Y, cluster, sigma2, lambda, clusteri){

  # filtering the cells in cluster i
  cellsn <- which(cluster == clusteri)

  # filtering the alpha values for cluster i
  # alphan <- alpha[[clusteri]]

  # selecting the Y values for the cells in cluster i
  Yn <- as.matrix(Y[,cellsn])

  # selecting sigma2 for cells in cluster i
  sigma2n <- sigma2[cellsn]

  Cr <- sum(cluster == clusteri)

  res2 <- c()
  min_H <- numeric(kstar+1)
  min_x <- Inf
  for(k in 0:kstar){
    m0 <- M - 1 -(k+1)*w
    mk <- permuteGeneral(0:m0, k + 1,
                         constraintFun = "sum",
                         comparisonFun = "==", limitConstraints = m0,
                         repetition = TRUE)
    Hlist <- apply(mk, 1,
                   FUN = function(x){postK_mk(k = k, m0 = m0, w = w, M = M, Yn = Yn,
                                              sigma2n = sigma2n, cellsn = cellsn, mk = x, Cr = Cr)})

    min_x <- logsumexp(Hlist, min_x)$min_x


  }


  for(k in 0:kstar){
    m0 <- M - 1 -(k+1)*w
    mk <- permuteGeneral(0:m0, k + 1,
                         constraintFun = "sum",
                         comparisonFun = "==", limitConstraints = m0,
                         repetition = TRUE)

    Hlist <- apply(mk, 1,
                   FUN = function(mkx){postK_mk(k = k, m0 = m0, w = w, M = M, Yn = Yn,
                                                sigma2n = sigma2n, cellsn = cellsn, mk = mkx, Cr = Cr)})

    res1 <- c()
    res1 <- exp(logsumexp(Hlist, min_x)$x)*exp(apply(mk, 1, function(z){(lfactorial(m0) - sum(lfactorial(z))) - m0*log(k+1)}))

    res2 <- c(res2, sum(res1)*pk(k = k, kstar = kstar, lambda = lambda))

  }

  return(c(0:kstar)[which.max(res2)])

}

# posterior of (m1, ..., mk+1)

postmk <- function(w, M, Y, K, cluster, sigma2, clusteri){

  # filtering the cells in cluster i
  cellsn <- which(cluster == clusteri)

  # filtering the alpha values for cluster i
  # alphan <- alpha[[clusteri]]

  # selecting the Y values for the cells in cluster i
  Yn <- as.matrix(Y[,cellsn])

  # selecting sigma2 for cells in cluster i
  sigma2n <- sigma2[cellsn]

  # number of change points
  k <- K[clusteri]

  Cr <- sum(cluster == clusteri)

  m0 <- M - 1 - (k+1)*w
  mk <- permuteGeneral(0:m0, k + 1,
                       constraintFun = "sum",
                       comparisonFun = "==", limitConstraints = m0,
                       repetition = TRUE)

  Hvalues <- apply(mk, 1,
                   function(mkx){postK_mk(k = k, m0 = m0, w = w, M = M, Yn = Yn,
                                          sigma2n = sigma2n, cellsn = cellsn, mk = mkx, Cr = Cr)})

  res1 <- exp(logsumexp(Hvalues)$x)*exp(apply(mk, 1, function(z){(lfactorial(m0) - sum(lfactorial(z))) - m0*log(k+1)}))

  xn <- sample(1:nrow(mk), 1, prob = res1/sum(res1))
  mkn <- mk[xn,]

  return(mkn)

}


# posterior alpha

postalphak <- function(M, Y, sigma2, K, Tl, cluster, clusteri){
  # filtering the cells in cluster i
  cellsn <- which(cluster == clusteri)

  # selecting the Y values for the cells in cluster i
  Yn <- as.matrix(Y[,cellsn])

  # selecting sigma2 for cells in cluster i
  sigma2n <- sigma2[cellsn]

  # number of change points for cluster i
  k <- K[clusteri]

  # Change point positions for cluster i
  Tln <- Tl[[clusteri]]

  Tlnn <- c()
  if(k == 0){
    Tlnn <- M
    Xn <- as.matrix(rep(1, Tlnn))
  } else{
    Xn <- apply(diag(k+1),2,function(x){rep(x,diff(Tln))})
  }

  Sigman <- solve((sum(sigma2n^(-1)))*t(Xn)%*%Xn)

  if(k == 0){
    mun <- Sigman*(sum(sapply(1:length(cellsn), function(x){(sigma2n[x]^(-1))*t(Xn)%*%Yn[,x]})))
    alphak <- rnorm(1,mun, sqrt(Sigman))
  } else{
    mun <- Sigman%*%as.matrix(rowSums(sapply(1:length(cellsn), function(x){(sigma2n[x]^(-1))*t(Xn)%*%Yn[,x]})))
    alphak <- sapply(1:(k+1), function(x){rnorm(1,mun[x,1],sqrt(diag(Sigman)[x]))})
  }
  return(alphak)
}


# posterior for lambda


full_cond <- function(kstar, lambda, cluster, al, bl, K, N){
  Nr <- as.vector(table(cluster))
  aln <- al + sum(Nr*K[sort(unique(cluster))])

  out <- (lambda^(aln - 1)*exp((-lambda/bl)))/(sum(sapply(0:kstar,FUN=function(l){lambda^l/factorial(l)}))^N)

}


update_lambda <- function(a = 4, b = 2, kstar, lambda, cluster, al, bl, K, N){

  lambda_prop <- rgamma(1, a, b)

  log.r <- log(full_cond(kstar=kstar, lambda=lambda_prop, cluster=cluster, al=al, bl=bl, K=K, N=N)) - log(full_cond(kstar=kstar, lambda=lambda, cluster=cluster, al=al, bl=bl, K=K, N=N)) + dgamma(lambda, a, b, log = TRUE) - dgamma(lambda_prop, a, b, log = TRUE)

  if(log(runif(1)) < log.r){lambda <- lambda_prop}

  return(lambda)
}

# posterior for alpha0

postalpha0 <- function(alpha0, a, b, N, cluster){

  d <- length(unique(cluster))

  kappa <- rbeta(1, alpha0 + 1, N)

  pik <- (a + d - 1)/(N*((1/b) - log(kappa)))

  if(pik > 1){
    pik = 1
  }

  x <- rbinom(1, 1, prob = pik)

  if(x == 0){
    alpha0n <- rgamma(1, a + d - 1, scale = ((1/b) - log(kappa))^(-1))
  } else{
    alpha0n <- rgamma(1, a + d, scale = ((1/b) - log(kappa))^(-1))
  }

  return(alpha0n)
}
