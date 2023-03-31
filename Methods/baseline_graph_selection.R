##baseline estimation and inference methods
#----------------------------------------------------#
find_graph_stb_tuning <- function(X, Obs, n_total, p, kept_nodes, method, zerofill, stb_thrs = 0.05){
  #covariance for full data and subsampled data
  X_kept <- X[ , kept_nodes]; Obs_kept <- Obs[ , kept_nodes]; p_kept <- length(kept_nodes)
  Sigma_hat <- array(rep(0, p_kept^2*21), c(21, p_kept, p_kept)); N <- array(rep(0, p_kept^2*21), c(21, p_kept, p_kept));
  start_t <- Sys.time()
  Sigma_out <- calc_Sigma_PSD(X_kept, Obs_kept, p_kept, Entryest = TRUE); 
  Sigma_hat[1, , ] <- Sigma_out$PSD_Sigma; Entryest_Sigma <- Sigma_out$Entryest_Sigma
  N[1, , ] <- Sigma_out$N
  graph_est <- list(); 
  lambda_list <- find_lambda_list(Sigma_hat[1, , ], method, p_kept, n_total, N[1, , ], zerofill)
  print(Sys.time() - start_t)
  print("Found lambda range!")
  start_t <- Sys.time()
  if(method == "CLIME"){
    nsubsample <- 10
    ntuning <- 10
  }else{
    nsubsample <- 20
    ntuning <- 20
  }
  for(i in 2 : (nsubsample+1)){
      start_t_tmp <- Sys.time()
      set.seed(i - 1)
      selected_ind <- as.logical(rbinom(n_total, 1, 0.8))
      Sigma_out <- calc_Sigma_PSD(X_kept[selected_ind, ], Obs_kept[selected_ind, ], p_kept)
      Sigma_hat[i, , ] <- Sigma_out$PSD_Sigma; N[i, , ] <- Sigma_out$N
      graph_est[[i - 1]] <- fit_baseline_graph(Sigma_hat[i, , ], method = method, 
                                    lambda = lambda_list * sqrt(n_total) / sqrt(sum(selected_ind)), 
                                    p_kept, sum(selected_ind), N[i, , ], zerofill)$graph_est_list
      print(Sys.time() - start_t_tmp)
      print(sprintf("Fitted graphs for %dth subsampling!", i))
  }
  select_n <- array(rep(0, ntuning * p_kept^2), c(ntuning, p_kept, p_kept))
  for(j in 1:ntuning){
    for(i in 1:nsubsample){
      select_n[j, , ] <- select_n[j, , ] + as.matrix(graph_est[[i]][[j]])
    }
  }
  stb_vec <- apply(select_n / nsubsample * (1 - select_n / nsubsample) * 2, 1, mean)
  if(stb_vec[1] > stb_thrs){
    choice_ind <- 1
  }else if(max(stb_vec) <= stb_thrs){
    choice_ind <- length(lambda_list)
  }else{
    choice_ind <- min(which(stb_vec > stb_thrs)) - 1
  }
  print(Sys.time() - start_t)
  print("Finished Tuning!")
  graph_est_final <- matrix(rep(0, p^2), p, p); #Theta_hat <- matrix(rep(0, p^2), p, p);
  out <- fit_baseline_graph(Sigma_hat[1, , ], method = method, 
                            lambda = lambda_list[choice_ind], 
                            p_kept, n_total, N[1, , ], zerofill)
  graph_est_final[kept_nodes, kept_nodes] <- as.numeric(out$graph_est_list[[1]])

  #result is nonzero for diagonal elements in clime!!!
  if(method == "glasso"){
    #Theta_hat[kept_nodes, kept_nodes] <- as.numeric(out$Theta_est[[1]])
    Theta_hat <- as.numeric(out$Theta_est[[1]])
    return(list(graph_est = graph_est_final, Sigma_hat = Sigma_hat[1, , ], Theta_hat = Theta_hat, lambda = lambda_list[choice_ind]))
  }else{
    return(list(graph_est = graph_est_final, Sigma_hat = Sigma_hat[1, , ], lambda = lambda_list[choice_ind]))
  }
}

find_lambda_list <- function(Sigma_hat, method, p, n_total, N, zerofill){
  lambda_max <- sqrt(log(p)/n_total)
  out <- fit_baseline_graph(Sigma_hat, method, lambda = lambda_max, p, n_total, N, zerofill)$graph_est_list[[1]]
  while(sum(out) != 0){
    lambda_max <- lambda_max * 2
    out <- fit_baseline_graph(Sigma_hat, method, lambda = lambda_max, p, n_total, N, zerofill)$graph_est_list[[1]]
  }
  if(lambda_max == sqrt(log(p)/n_total)){
      while(sum(out) == 0){
        lambda_max <- lambda_max / 2
        out <- fit_baseline_graph(Sigma_hat, method, lambda = lambda_max, p, n_total, N, zerofill)$graph_est_list[[1]]
      }
      lambda_max <- lambda_max * 2
  }
  if(method == "CLIME"){
    lambda_list <- exp(log(lambda_max)-log(10) / 9 * (0 : 9))
  }else{
    lambda_list <- exp(log(lambda_max)-log(10) / 19 * (0 : 19))
  }
  return(lambda_list)
}

fit_baseline_graph <- function(Sigma_hat, method, lambda_list, p, n, N, zerofill, thrs_c = 5){
  if(zerofill){
    Sigma_hat[N <= thrs_c * log(p)] <- 0
  }
  if(method == "nblasso_and"){
    out <- huge(Sigma_hat, method = "mb", sym = "and", lambda = lambda_list, verbose = FALSE)
    return(list(graph_est_list = out$path))
  }else if(method == "nblasso_or"){
    out <- huge(Sigma_hat, method = "mb", sym = "or", lambda = lambda_list, verbose = FALSE)
    return(list(graph_est_list = out$path))
    }else if(method == "glasso"){
      out <- huge(Sigma_hat, method = "glasso", lambda = lambda_list, verbose = FALSE)
      return(list(graph_est_list = out$path, Theta_est_list = out$icov))
      }else if(method == "CLIME"){
        linsolver <- ifelse(n > p, "simplex", "primaldual"); 
        out <- clime(Sigma_hat, sigma = TRUE, lambda = lambda_list, linsolver = linsolver)$Omega
        return(list(graph_est_list  = lapply(out, function(x){diag(x)=0; x!=0}), Theta_est = out))
    }
}

#inference baseline: debiased graphical lasso (Jankova and Van de Geer, 2015); Ren et al 2015; + Bonferroni / Holm's
#inference baseline: GGM FDR control

test_glasso <- function(Sigma_hat, Theta_hat, n, p, kept_nodes, signif_pval, testing_thrs = 0){
  Theta_tilde <- 2 * Theta_hat - Theta_hat %*% Sigma_hat %*% Theta_hat
  var_est <- diag(Theta_hat) %*% t(diag(Theta_hat)) + Theta_hat^2
  test <- matrix(rep(0, p^2), p, p)
  if(testing_thrs != 0){
    test[kept_nodes, kept_nodes] <- (abs(Theta_tilde) - testing_thrs) / sqrt(var_est) * sqrt(n)
    diag(test) <- 0
    p_val <- pmin(2 * (1 - pnorm(test)), 1)
  }else{
    test[kept_nodes, kept_nodes] <- Theta_tilde / sqrt(var_est) * sqrt(n)
    diag(test) <- 0
    p_val <- 2 * (1 - pnorm(abs(test)))
  }

  
  #Bonferroni correction
  p_val_Bonferroni <- pmin(p_val*p*(p-1)/2, 1)
  signif_graph_Bonferroni <- p_val_Bonferroni < signif_pval
  
  #Holm-Bonferroni correction
  ind <- matrix(1 : (p ^ 2), p, p)
  p_val_vec <- p_val[upper.tri(p_val)]; order_p_val <-order(p_val_vec)
  k <- min(which(p_val_vec[order_p_val] > signif_pval / seq(from = p*(p-1)/2, to = 1, by = -1)))
  signif_ind_holm <- ind[upper.tri(ind)][order_p_val[1 : (k - 1)]]
  signif_graph_holm <- matrix(rep(0, p^2), p, p)
  signif_graph_holm[signif_ind_holm] <- 1; signif_graph_holm <- signif_graph_holm + t(signif_graph_holm)
  
  return(list(p_val = p_val, p_val_Bonferroni = p_val_Bonferroni, signif_graph_Bonferroni = signif_graph_Bonferroni, 
              signif_graph_holm = signif_graph_holm))
}
