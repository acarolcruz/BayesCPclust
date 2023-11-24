save_sim_random <- function(sim, seed, scenario, casenumber, M, N, w, d, as, bs, al, bl, a, b, alpha0,
                        lambda, maxIter, startpoint_sampling, par.values, emp, initial_val){

  tt <- Sys.time()
  out <- run_gibbs_random(M = M, N = N, w = w, d = d, as = as, bs = bs, al = al, bl = bl, a = a, b = b,
                          alpha0 = alpha0, lambda = lambda, maxIter = maxIter, par.values = par.values,
                          emp = emp, initial_val = initial_val, seed = seed + sim)

  tt <- difftime(Sys.time(), tt, units = 'mins')

  data <- out[[1]]

  res <- out[[2]]

  save(res, file = paste0(scenario, "/case", casenumber, "/Sim", sim, "IterGibbs_all.RData"))

  index = seq(from = round((maxIter/2)+1), to = maxIter, by = startpoint_sampling)

  res1 <- list(Ck = res$Ck[index], Nclusters = res$Nclusters[index], lambda = res$lambda[index],
               alpha0 = res$alpha0[index], clusters = res$clusters[index,], sigma2 = res$sigma2[index,],
               Tl = res$Tl[index], alphal = res$alphal[index])

  save(res1, file = paste0(scenario, "/case", casenumber, "/Sim", sim, "IterGibbs_final.RData"))

  #computing posterior means and modes for all parameters except alpha and T

  out2 <- list(Nclusters = Mode(res1$Nclusters),
               lambda = mean(res1$lambda), alpha0 = mean(res1$alpha0),
               clusters = apply(res1$clusters, 2, Mode), sigma2 = apply(res1$sigma2, 2, mean))

  # merging Ck for selected iterations (considering we can have different sizes in each iteration)

  Ckres <- do.call("rbind", lapply(res1$Ck, "[", seq(1:max(sapply(res1$Ck, length)))))

  out2$Ck <- apply(Ckres[,sort(unique(out2$clusters))], 2, Mode)

  # merging selected iteration for Tl and alpha (considering we can have different sizes in each iteration)

  Tlist_tab <- lapply(1:length(res1$Tl),
                      function(x){do.call("list",
                                          lapply(res1$Tl[[x]][sort(unique(res1$clusters[x,]))], "[",
                                                 seq(1:max(sapply(res1$Tl[[x]][sort(unique(res1$clusters[x,]))], length)))))})



  # creating list for each cluster, considering the most frequent number of clusters
  Tlist <- lapply(1:out2$Nclusters, function(x){lapply(Tlist_tab , "[[", x)})


  Tvalues <- lapply(1:length(Tlist),
                    function(x){do.call("rbind",
                                        lapply(Tlist[[x]], "[", seq(1:max(sapply(Tlist[[x]], length)))))})

  a_tab <- lapply(1:length(res1$alphal),
                  function(x){do.call("list",
                                      lapply(res1$alphal[[x]][sort(unique(res1$clusters[x,]))], "[",
                                             seq(1:max(sapply(res1$alphal[[x]][sort(unique(res1$clusters[x,]))], length)))))})


  alist <- lapply(1:out2$Nclusters, function(x){lapply(a_tab, "[[", x)})


  avalues <- lapply(1:length(alist),
                    function(x){do.call("rbind",
                                        lapply(alist[[x]], "[", seq(1:max(sapply(alist[[x]], length)))))})



  # computing posterior means for alpha and T
  # computing posterior means for alpha and T
  res_alpha <- list()
  res_Tl <- list()
  for(i in 1:out2$Nclusters){
    ck <- out2$Ck[i]
    if(ck == 0){
      alphapos <- mean(avalues[[i]])
      Tlpos <- Mode(Tvalues[[i]])
    } else{
      Tl <- Tvalues[[i]][,-c(1,ncol(Tvalues[[i]]))]
      alphapos <- apply(avalues[[i]][,1:(ck+1)], 2, mean)
      Tlpos <- apply(Tl[,1:ck], 2, Mode)
    }
    res_alpha[[i]] <- alphapos
    res_Tl[[i]] <- Tlpos
  }

  out2$alpha <- res_alpha
  out2$Tl  <- res_Tl

  #saving data and true parameters

  for(param in names(data)){
    if(is.matrix(data[[param]])){
      write(data[[param]], file = paste0(scenario, "/case", casenumber, "/Sim", sim, param, ".txt"), ncolumns = ncol(data[[param]]))
    } else{
      write(data[[param]], file = paste0(scenario, "/case", casenumber, "/Sim", sim, param, ".txt"), ncolumns = length(data[[param]]))
    }
  }

  #saving posterior parameters

  for(g in names(out2)){
    if(is.matrix(out2[[g]])){
      write(out2[[g]], file = paste0(scenario, "/case", casenumber, "/Sim", sim, g, ".txt"), ncolumns = ncol(out2[[g]]))

    } else{
      if(is.list(out2[[g]])){
        for(j in 1:length(out2[[g]])){
          # browser()
          write(out2[[g]][[j]], file = paste0(scenario, "/case", casenumber, "/Sim", sim, g, "cluster", j, ".txt"))
        }
      } else{
        # browser()
        write(out2[[g]], file = paste0(scenario, "/case", casenumber, "/Sim", sim, g, ".txt"), ncolumns = length(out2[[g]]))
      }
    }
  }

  cat("Simulations saved in folder")
  cat("Processing time", tt)

  return(tt)
}
