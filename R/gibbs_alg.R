# Run gibbs algorithm
gibbs_alg <- function(N, w, M, K, Tl, cluster, alpha, sigma2, bs = 1000, as = 2, al = 2, bl = 1000, a = 2, b = 1000,
                      alpha0 = 1/100, kstar, lambda, Y, d, maxIter = 10000){

  # creating structure to save output
  clustervalues <- matrix(NA, maxIter, N)
  sigma2values <- matrix(NA, maxIter, N)
  dvalues <- c()
  lambdavalues <- c()
  alpha0values <- c()
  ckvalues <- NULL

  resTk <- NULL
  resak <- NULL

  # run the Gibbs sampler
  for(iter in 1:maxIter){
    resTk[[iter]] <- list()
    resak[[iter]] <- list()
    # take a Gibbs step at each data point
    for(cell in 1:N){

      k <- K[cluster[cell]]

      yn <- Y[,cell]

      alphan <- alpha[[cluster[cell]]]

      Tln <- Tl[[cluster[cell]]]

      sigma2n <- sigma2[cell]

      m0 <- M - 1 - (k + 1)*w

      qn0.res <- qn0(alpha0=alpha0, w=w, N=N, M=M, bs=bs, as=as, kstar=kstar, lambda=lambda, Yn=yn)
      qnj.res <- qnj(N=N, M=M, as=as, bs=bs, Yn=yn, alpha=alpha, cluster=cluster, Tl=Tl, K=K)

      p0 <- qn0.res/(qn0.res + sum(qnj.res))

      if(rbinom(1, 1, 1 - p0) == 0){ #

        ccurrent <- cluster[cell]

        #generate a new cluster
        cluster[cell] <- max(cluster) + 1

        if((sum(cluster == ccurrent) == 0) & (cluster[cell] != ccurrent)){
          d <- d - 1
        }

        knew <- postK(kstar=kstar, w=w, M=M, Y=Y, cluster=cluster,
                      sigma2=sigma2, lambda=lambda, clusteri=cluster[cell])

        # update K
        K[cluster[cell]] <- knew

        mknew <- postmk(w=w, M=M, Y=Y, K=K, cluster=cluster, sigma2=sigma2, clusteri=cluster[cell])

        if(K[cluster[cell]] == 0){
          Tnew <- mknew + w + 1
          Tl[[cluster[cell]]] <- Tnew
        } else{
          Tnew <- cumsum(c(1,mknew + w))[-c(1,length(c(1,mknew)))]
          # update Tl
          Tl[[cluster[cell]]] <- c(1, Tnew, M + 1)
        }

        alphanew <- postalphak(M=M, Y=Y, sigma2=sigma2, K=K, Tl=Tl, cluster=cluster, clusteri=cluster[cell])


        # Update alpha
        alpha[[cluster[cell]]] <- alphanew

        # Update number of clusters
        d <- d + 1


      } else {
        # assigned to an existing cluster

        #current cluster
        ccurrent <- cluster[cell]

        pj <- qnj.res

        cluster[cell] <- c(1:max(cluster))[which.max(pj)] # debugged

        if((sum(cluster == ccurrent) == 0) & (cluster[cell] != ccurrent)){
          d <- d - 1
        }

      }
    }
    cat("Step 1 completed for iteration", iter, "\n")


    # step 2 update sigma2n
    for(cell in 1:N){
      yn <- Y[,cell]

      k <- K[cluster[cell]]

      alphan <- alpha[[cluster[cell]]]

      Tln <- Tl[[cluster[cell]]]

      # bs = rate parameter in the inverse Gamma
      sigma2new <- possigma2n(as=as, bs=bs, M=M, Yn=yn, k=k, Tln=Tln, alphan=alphan)

      sigma2[cell] <- sigma2new
    }
    cat("Step 2 completed for iteration", iter, "\n")

    # Step 3: Update cluster parameters
    for(cl in 1:length(unique(cluster))){

      # sample a different K, T, alpha from their posteriors

      knew <- postK(kstar=kstar, w=w, M=M, Y=Y, cluster=cluster, sigma2=sigma2, lambda=lambda,
                    clusteri=sort(unique(cluster))[cl])

      # update K
      K[sort(unique(cluster))[cl]] <- knew

      mknew <- postmk(w=w, M=M, Y=Y, K=K, cluster=cluster, sigma2=sigma2, clusteri=sort(unique(cluster))[cl])

      if(K[sort(unique(cluster))[cl]] == 0){
        Tnew <- mknew + w + 1

        Tl[[sort(unique(cluster))[cl]]] <- Tnew
      } else{
        Tnew <- cumsum(c(1,mknew + w))[-c(1,length(c(1,mknew)))]

        # update Tl
        Tl[[sort(unique(cluster))[cl]]] <- c(1, Tnew, M + 1)
      }

      alphanew <- postalphak(M=M, Y=Y, sigma2=sigma2, K=K, Tl=Tl, cluster=cluster, clusteri=sort(unique(cluster))[cl])

      # Update alpha
      alpha[[sort(unique(cluster))[cl]]] <- alphanew

    }

    cat("Step 3 completed for iteration", iter, "\n")

    # Step 4: Update lambda

    lambda_new <- update_lambda(kstar=kstar, lambda=lambda, cluster=cluster, al=al, bl=bl, K=K, N=N)
    lambda <- lambda_new

    cat("Step 4 completed for iteration", iter, "\n")

    # Step 5: Update alpha0

    alpha0new <- postalpha0(alpha0, a, b, N, cluster)

    alpha0 <- alpha0new

    cat("Step 5 completed for iteration", iter, "\n")

    clustervalues[iter,] <- cluster
    sigma2values[iter,] <- sigma2
    dvalues <- c(dvalues, d)
    lambdavalues <- c(lambdavalues, lambda)
    alpha0values <- c(alpha0values, alpha0)
    ckvalues[[iter]] <- K


    for(i in 1:length(alpha)){
      resTk[[iter]][[i]] <- NULL
      resak[[iter]][[i]] <- NULL
      #resTk[[i]] <- Tm
      #resak[[i]] <- alpham
      resTk[[iter]][[i]] <- Tl[[i]]
      resak[[iter]][[i]] <- alpha[[i]]

    }

  }
  return(list(Ck = ckvalues, Nclusters = dvalues, lambda = lambdavalues, alpha0 = alpha0values,
              clusters = clustervalues, sigma2 = sigma2values, Tl = resTk, alphal = resak))

}

