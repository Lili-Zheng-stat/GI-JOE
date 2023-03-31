
##fitting functions for general missing patterns or independent missingness
#----------------------------------------------------#
library(MTS)
library(reticulate)
library(rTensor)

#----------------------------------------------------#
find_p_val_graph <- function(X, Obs, n_total, p, thrs_c, signif_pval, onetuning = FALSE, symmetry = FALSE, ...){
  #given data X, Obs, n_total, p, screening threshold constant, 
  #find p value matrix, adjusted and unadjusted; and significant graph based on signif_p_val
  #if onetuning = true, only one tuning constant is used for all pairs of nodes
  input <- list(...);
  
  #screening  
  N <- t(Obs)%*%Obs
  ind_order <- order(apply(N,1,mean), decreasing = TRUE)
  p_screened <- p
  while(sum(N[ind_order[1:p_screened], ind_order[1:p_screened]]<=thrs_c*log(p))>0){
    p_screened <- p_screened - 1
  }
  
  if(p_screened < 2){
    print("Only one node left!")
    return(list(success = FALSE, p_screened = p_screened))
  }else{
    X_screened <- X[, ind_order[1:p_screened]]
    Obs_screened <- Obs[, ind_order[1:p_screened]]
    
    #testing
    if(length(input)>0){
      test_out <- graph_testing(X_screened, Obs_screened, n_total, p_screened, input[[1]])
    }else if(!onetuning){
      test_out <- graph_testing(X_screened, Obs_screened, n_total, p_screened)
    }else{
      test_out <- graph_testing_onetuning(X_screened, Obs_screened, n_total, p_screened, symmetry)
    }
    test <- matrix(rep(0, p^2), p, p)
    
    if(!symmetry){
      test[ind_order[1:p_screened], ind_order[1:p_screened]] <- test_out$test+t(test_out$test)
      p_val <- matrix(pmin(1,(1-pnorm(abs(test)))*2), p, p)
    }else{
      test[ind_order[1:p_screened], ind_order[1:p_screened]] <- test_out$test
      p_val <- matrix(pmin(1,(1-pnorm(abs(test)))*2), p, p)
      p_val <- pmin(p_val, t(p_val))*2
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
    
    #FDR control
    k <- max(which(p_val_vec[order_p_val] <= signif_pval * (1 : (p*(p-1)/2)) / (p*(p-1)/2)))
    if(sum(p_val_vec[order_p_val] <= signif_pval * (1 : (p*(p-1)/2)) / (p*(p-1)/2))>0){
      rho <- signif_pval * k / (p*(p-1)/2)
    }else{
      rho <- signif_pval / (p*(p-1)/2)
    }
    t_max <- sqrt(2 * log(p * (p - 1) / 2) - 2 * log(log(p * (p - 1) / 2))); rho_min <- 2 * (1 - pnorm(t_max))
    if(rho < rho_min){
      rho <- 2 * (1 - pnorm(sqrt(2 * log(p * (p - 1) / 2))))
    }
    signif_ind_FDR <- ind[upper.tri(ind)][which(p_val_vec <= rho)]
    signif_graph_FDR <- matrix(rep(0, p^2), p, p)
    signif_graph_FDR[signif_ind_FDR] <- 1; signif_graph_FDR <- signif_graph_FDR + t(signif_graph_FDR)
    
    #estimated graph 
    graph_est <- matrix(rep(0, p^2), p, p)
    graph_est[ind_order[1:p_screened], ind_order[1:p_screened]] <- test_out$beta_hat != 0 
    graph_est_OR <- pmax(graph_est, t(graph_est))
    graph_est_AND <- pmin(graph_est, t(graph_est))
    
    return(list(success = TRUE, test_out = test_out, test_mat = test, p_val_mat = p_val, 
                p_val_Bonferroni = p_val_Bonferroni, signif_graph_Bonferroni = signif_graph_Bonferroni,
                signif_graph_holm = signif_graph_holm, signif_graph_FDR = signif_graph_FDR,
                graph_est_OR = graph_est_OR, graph_est_AND = graph_est_AND, kept_nodes = ind_order[1:p_screened]))
    
  }
}


#----------------------------------------------------#
graph_testing <- function(X, Obs, n_total, p, ...){
  ## Given data X, observational pattern Obs, pairwise sample sizes N, 
  ## output the test statistics and p-values for each edge (neighborhood lasso on the node with smaller index; 
  ## if wanting neighborhood lasso for the other node, can use a permutation on rows of X before input and 
  ## apply the reversed permutation on computed p-values).
  # tuning_c can be an additional input; otherwise, chosen by stability
  input <- list(...);
  
  #find entrywise estimate of Sigma over full data and random subsamples; also find pairwise sample sizes for subsamples
  Sigma_hat <- array(rep(0, p^2*21), c(21, p, p)); N <- array(rep(0, p^2*21), c(21, p, p));
  start_t <- Sys.time()
  Sigma_out <- calc_Sigma_PSD(X, Obs, p, Entryest = TRUE); 
  Sigma_hat[1, , ] <- Sigma_out$PSD_Sigma; Entryest_Sigma <- Sigma_out$Entryest_Sigma
  N[1, , ] <- Sigma_out$N
  for(i in 2 : 21){
    set.seed(i - 1)
    selected_ind <- as.logical(rbinom(n_total, 1, 0.8))
    Sigma_out <- calc_Sigma_PSD(X[selected_ind, ], Obs[selected_ind, ], p)
    Sigma_hat[i, , ] <- Sigma_out$PSD_Sigma; N[i, , ] <- Sigma_out$N
  }
  print(Sys.time() - start_t)
  print("Found covariance estimates for all subsampling!")
  
  #find N4 (exhausting memory for p = 200!)
  #N4 <- m4_sum(Obs)
  
  #loop over all pairs (a, b): neighborhood lasso for a and debiasing vector for (a,b)
  eta <- 1;tol <- 0.00001; lambda_scale <- sqrt(log(p)/apply(N[1, , ], 1, min))
  beta_hat <- matrix(rep(0, p^2), p, p); theta_hat <- array(rep(0, p^3), c(p, p, p)); beta_tilde <- matrix(rep(0, p^2), p, p)
  var_est <- matrix(rep(0, p^2), p, p); tuning_c_vec <- rep(0, p); tuning_c_mat <- matrix(rep(0, p^2), p, p); test <- matrix(rep(0, p^2), p, p)
  start_t <- Sys.time()
  for(node_a in 1 : (p - 1)){
    if(length(input) > 0){
      fit_out1 <- fit_nb(Sigma_hat[1, , ], p, node_a, lambda_scale[-node_a] * input[[1]], eta, tol)
    }else{
      fit_out1 <- est_nb_tuning(Sigma_hat, N, p, stb_thrs = 0.05, node_a, lambda_scale[-node_a], eta, tol)
      tuning_c_vec[node_a] <- fit_out1$tuning_c
    }
    beta_hat[node_a, -node_a] <- fit_out1$beta_hat
    ind <- (1 : p)[-node_a]
    for(node_b in (node_a + 1) : p){
      node_b_ind <- which(ind == node_b)
      lambda_scale_tmp <- sqrt(log(p - 1)/apply(N[1, -node_a, -node_a], 1, min))
      if(length(input) > 0){
        fit_out2 <- fit_nb(Sigma_hat[1, -node_a, -node_a], p, node_b_ind, lambda_scale_tmp[-node_b_ind] * input[[1]], eta, tol)
      }else{
        fit_out2 <- est_nb_tuning(Sigma_hat[, -node_a, -node_a], N[, -node_a, -node_a], p - 1, stb_thrs = 0.05, node_b_ind, lambda_scale_tmp[-node_b_ind], eta, tol)
        tuning_c_mat[node_a, node_b] <- fit_out2$tuning_c
      }
      theta_hat[node_a, node_b, -c(node_a, node_b)] <- fit_out2$beta_hat
      tau <- as.numeric(1 / (Sigma_hat[1, node_b, node_b] - Sigma_hat[1, node_b, ] %*% theta_hat[node_a, node_b, ]))
      theta_hat[node_a, node_b, ] <- -theta_hat[node_a, node_b,] * tau
      theta_hat[node_a, node_b, node_b] <- tau
      #debiased neighborhood lasso (node_a, node_b)
      beta_tilde[node_a, node_b] = as.numeric(beta_hat[node_a, node_b] - t(theta_hat[node_a, node_b, ]) %*% (Entryest_Sigma %*% beta_hat[node_a, ] - Entryest_Sigma[, node_a]))
      #variance estimation
      beta_hat_tmp <- -beta_hat[node_a, ];
      beta_hat_tmp[node_a] <- 1;
      var_est[node_a, node_b] <- find_var(Entryest_Sigma, N[1, , ], Obs, beta_hat_tmp, theta_hat[node_a, node_b, ])
      #test statistic
      test[node_a, node_b] <- beta_tilde[node_a, node_b]/sqrt(var_est[node_a, node_b])
    }
  }
  print(Sys.time() - start_t)
  print("Graph testing using nblasso is done!")
  return(list(beta_hat = beta_hat, theta_hat = theta_hat, beta_tilde = beta_tilde, tuning_c1 = tuning_c_vec, tuning_c2 = tuning_c_mat, var_est = var_est, test = test))
}

#----------------------------------------------------#
graph_testing_onetuning <- function(X, Obs, n_total, p, symmetry){
  ## Choose one fixed tuning_c for all pairs of node; 
  ## if symmetry = true, for each pair of nodes, perform two tests in a symmetric way;
  ## otherwise, perform one test based on the order
  ## Given data X, observational pattern Obs, pairwise sample sizes N, 
  ## output the test statistics and p-values for each edge (neighborhood lasso on the node with smaller index; 
  ## if wanting neighborhood lasso for the other node, can use a permutation on rows of X before input and 
  ## apply the reversed permutation on computed p-values).
  
  #find entrywise estimate of Sigma over full data and random subsamples; also find pairwise sample sizes for subsamples
  Sigma_hat <- array(rep(0, p^2*21), c(21, p, p)); N <- array(rep(0, p^2*21), c(21, p, p));
  start_t <- Sys.time()
  Sigma_out <- calc_Sigma_PSD(X, Obs, p, Entryest = TRUE); 
  Sigma_hat[1, , ] <- Sigma_out$PSD_Sigma; Entryest_Sigma <- Sigma_out$Entryest_Sigma
  N[1, , ] <- Sigma_out$N
  for(i in 2 : 21){
    set.seed(i - 1)
    selected_ind <- as.logical(rbinom(n_total, 1, 0.8))
    Sigma_out <- calc_Sigma_PSD(X[selected_ind, ], Obs[selected_ind, ], p)
    Sigma_hat[i, , ] <- Sigma_out$PSD_Sigma; N[i, , ] <- Sigma_out$N
  }
  print(Sys.time() - start_t)
  print("Found covariance estimates for all subsampling!")
  
  #perform tuning and output neighborhood lasso result for each node
  eta <- 1;tol <- 0.00001; lambda_scale <- sqrt(log(p)/apply(N[1, , ], 1, min))
  tuning_fitting_results <- est_graph_tuning(Sigma_hat, N, p, stb_thrs = 0.05, lambda_scale, eta, tol)
  beta_hat <- tuning_fitting_results$nblasso_combined
  tuning_c <- tuning_fitting_results$tuning_c
  print("Finished tuning!")
  
  #loop over all pairs (a, b): debiasing vector for (a,b); 
  #if symmetry = true, (a, b) and (b, a) are both considered; otherwise, consider (a, b) for a < b
  theta_hat <- array(rep(0, p^3), c(p, p, p)); beta_tilde <- matrix(rep(0, p^2), p, p)
  var_est <- matrix(rep(0, p^2), p, p); test <- matrix(rep(0, p^2), p, p)
  start_t <- Sys.time()
  for(node_a in 1 : (p - 1)){
    ind <- (1 : p)[-node_a]
    for(node_b in ind){
      if(node_b > node_a || symmetry){
        node_b_ind <- which(ind == node_b)
        lambda_scale_tmp <- sqrt(log(p - 1)/apply(N[1, -node_a, -node_a], 1, min))
        fit_out2 <- fit_nb(Sigma_hat[1, -node_a, -node_a], p - 1, node_b_ind, lambda_scale_tmp[-node_b_ind] * tuning_c, eta, tol)
        theta_hat[node_a, node_b, -c(node_a, node_b)] <- fit_out2$beta_hat
        tau <- as.numeric(1 / (Sigma_hat[1, node_b, node_b] - Sigma_hat[1, node_b, ] %*% theta_hat[node_a, node_b, ]))
        theta_hat[node_a, node_b, ] <- -theta_hat[node_a, node_b,] * tau
        theta_hat[node_a, node_b, node_b] <- tau
        #debiased neighborhood lasso (node_a, node_b)
        beta_tilde[node_a, node_b] = as.numeric(beta_hat[node_a, node_b] - t(theta_hat[node_a, node_b, ]) %*% (Entryest_Sigma %*% beta_hat[node_a, ] - Entryest_Sigma[, node_a]))
      }
    }
  }
  print(Sys.time() - start_t)
  print("Graph estimation using neighborhood lasso is done!")
  #variance estimation
  #var_results <- foreach(a = 1 : (p - 1)) %:% 
  #  foreach(b = (a + 1) : p, .combine='c') %dopar% {
  #    wd <- "/Applications/Files/research/graph-quilting/graph-learning-varying-sample-size/simulations/inference_simulations"
  #    setwd(wd)
  #    source("fitting_functions_missing.R")
  #    beta_hat_tmp <- -beta_hat[a, ];
  #    beta_hat_tmp[a] <- 1;
  #    find_var(Entryest_Sigma, N[1, , ], Obs, beta_hat_tmp, theta_hat[a, b, ])
  #  }
  #for(a in 1 : (p -1 )){
  #  var_est[a, (a + 1) : p] <- var_results[[a]]
  #}
  start_t <- Sys.time()
  for(a in 1 : (p - 1)){
    ind <- (1 : p)[-a]
    for(b in ind){
      if(b > a || symmetry){
        beta_hat_tmp <- -beta_hat[a, ]; beta_hat_tmp[a] <- 1; beta_hat_tmp[abs(beta_hat_tmp) / max(abs(beta_hat_tmp)) <= 10^(-5)] <- 0;
        theta_tmp <- theta_hat[a, b, ]; theta_tmp[abs(theta_tmp) / max(abs(theta_tmp)) <= 10^(-5)] <- 0;
        var_est[a, b] <- find_var(Entryest_Sigma, N[1, , ], Obs, beta_hat_tmp, theta_tmp)
      }
    }
    print(sprintf("finished variance estimation for node %d", a))
  }
  print(Sys.time() - start_t)
  print("Variance estimation is done!")
  #test statistic
  if(symmetry){
    test <- beta_tilde/sqrt(var_est)
    diag(test) <- 0
    test <- matrix(test, p ,p)
  }else{
    test[upper.tri(test)] <- beta_tilde[upper.tri(beta_tilde)]/sqrt(var_est[upper.tri(beta_tilde)])
  }
  return(list(beta_hat = beta_hat, theta_hat = theta_hat, beta_tilde = beta_tilde, tuning_c = tuning_c, var_est = var_est, test = test))
}

#----------------------------------------------------#
est_graph_tuning <- function(Sigma_hat, N, p, stb_thrs, lambda_scale, eta, tol){
  ## perform neighborhood lasso for all nodes with tuning parameter selected based on stability
  ## Sigma_hat includes the PSD covariance estimate using the full data and 20 subsampled data;
  ## lambda_scale is the p-dimensional scale of penality parameters for the whole data set
  ## N is 21*p*p, N[1,,] is the pairwise sample size of full data; N[i,,] for i=2,...,21 is the pairwise sample size of subsampled data
  
  #find range of tuning_c
  #tuning_c_max is the approximately the smallest tuning parameter that leads to an empty graph
  tuning_c_max <- 1
  node_a <- 1
  out <- fit_nb(Sigma_hat[1, , ], p, node_a, lambda_scale[-node_a] * tuning_c_max, eta, tol)
  while(sum(out$beta_hat != 0) > 0 || node_a < p){
    if(sum(out$beta_hat != 0) > 0){
      tuning_c_max <- tuning_c_max * 2
    }else{
       node_a <- node_a + 1;
    }
    out <- fit_nb(Sigma_hat[1, , ], p, node_a, lambda_scale[-node_a] * tuning_c_max, eta, tol)
  }
  if(tuning_c_max == 1){
    while(sum(out$beta_hat != 0) == 0){
      if(node_a < p){
        node_a <- node_a + 1
      }else{
        tuning_c_max <- tuning_c_max / 2
        node_a <- 1
      }
      out <- fit_nb(Sigma_hat[1, , ], p, node_a, lambda_scale[-node_a] * tuning_c_max, eta, tol)
    }
    tuning_c_max <- tuning_c_max * 2
  }
  tuning_c_vec <- exp(log(tuning_c_max)-log(10) / 19 * (0 : 19))
  
  #Given the vector of tuning parameters, fit graphs for subsampled data and return stabilities
  tuning_out <- graph_stb_tuning(Sigma_hat[2 : 21, , ], N[2 : 21, , ], p, tuning_c_vec, eta, tol)
  
  if(tuning_out$stb_vec[1] > stb_thrs){
    choice_ind <- 1
  }else if(max(tuning_out$stb_vec) <= stb_thrs){
    choice_ind <- length(tuning_c_vec)
  }else{
    choice_ind <- min(which(tuning_out$stb_vec > stb_thrs)) - 1
  }
  
  nblasso_combined <- matrix(rep(0, p^2), p, p); out_final <- list()
  for(node_a in 1 : p){
    out_final[[node_a]] <- fit_nb(Sigma_hat[1, , ], p, node_a, lambda_scale[-node_a] * tuning_c_vec[choice_ind], eta, tol)
    nblasso_combined[node_a, -node_a] <- out_final[[node_a]]$beta_hat
  }
  return(list(out_final = out_final, nblasso_combined = nblasso_combined, stb_vec = tuning_out$stb_vec, 
              theta_stb = tuning_out$theta_stb, tuning_c_vec = tuning_c_vec, tuning_c = tuning_c_vec[choice_ind]))
}
#----------------------------------------------------#
graph_stb_tuning <- function(Sigma_hat, N, p, tuning_c_vec, eta, tol){
  #use stability tuning; 
  n_subsample <- dim(Sigma_hat)[1]
  graph_est_arr <- array(rep(0, length(tuning_c_vec) * n_subsample * p * p), dim = c(length(tuning_c_vec), n_subsample, p, p))
  for(kk in 1 : n_subsample){
    lambda_scale_tmp <- sqrt(log(p)/apply(N[kk, , ],1,min))
    for(j in 1 : length(tuning_c_vec)){
      for(node_a in 1 : p){
        out <- fit_nb(Sigma_hat[kk, , ], p, node_a, lambda_scale_tmp[-node_a] * tuning_c_vec[j], eta, tol)
        graph_est_arr[j, kk, node_a, -node_a] <- out$beta_hat != 0
      }
    }
  }
  theta_stb <- apply(graph_est_arr, c(1, 3, 4), mean)
  stb_vec <- apply(2 * theta_stb * (1 - theta_stb), 1, mean)
  return(list(theta_stb = theta_stb, stb_vec = stb_vec))
}

est_nb_tuning <- function(Sigma_hat, N, p, stb_thrs, node_a, lambda_scale, eta, tol){
  ## estimate neighborhood with tuning parameter selected based on stability
  ## Sigma_hat includes the PSD covariance estimate using the full data and 20 subsampled data;
  ## lambda_scale is the scale of penality parameters for the whole data (sqrt(log p/N))
  ## N is 21*p*p, N[1,,] is the pairwise sample size of full data; N[i,,] for i=2,...,21 is the pairwise sample size of subsampled data
  
  #find range of tuning_c
  tuning_c_max <- 1
  out <- fit_nb(Sigma_hat[1, , ], p, node_a, lambda_scale * tuning_c_max, eta, tol)
  while(sum(out$beta_hat != 0) > 0){
    tuning_c_max <- tuning_c_max * 2
    out <- fit_nb(Sigma_hat[1, , ], p, node_a, lambda_scale * tuning_c_max, eta, tol)
  }
  if(tuning_c_max == 1){
    while(sum(out$beta_hat != 0) == 0){
      tuning_c_max <- tuning_c_max / 2
      out <- fit_nb(Sigma_hat[1, , ], p, node_a, lambda_scale * tuning_c_max, eta, tol)
    }
    tuning_c_max <- tuning_c_max * 2
  }
  tuning_c_vec <- exp(log(tuning_c_max)-log(10) / 19 * (0 : 19))
  tuning_out <- stb_tuning(Sigma_hat[2 : 21, , ], N[2 : 21, , ], p, node_a, tuning_c_vec, eta, tol)
  if(tuning_out$stb_vec[1] > stb_thrs){
    choice_ind <- 1
  }else if(max(tuning_out$stb_vec) <= stb_thrs){
    choice_ind <- length(tuning_c_vec)
  }else{
    choice_ind <- min(which(tuning_out$stb_vec > stb_thrs)) - 1
  }
  
  out_final <- fit_nb(Sigma_hat[1, , ], p, node_a, lambda_scale * tuning_c_vec[choice_ind], eta, tol)
  out_final$stb_vec <- tuning_out$stb_vec
  out_final$theta_stb <- tuning_out$theta_stb
  out_final$tuning_c_vec <- tuning_c_vec
  out_final$tuning_c <- tuning_c_vec[choice_ind]
  return(out_final)
}

stb_tuning <- function(Sigma_hat, N, p, node_a, tuning_c_vec, eta, tol){
  #use stability tuning; 
  n_subsample <- dim(Sigma_hat)[1]
  nbhd_est <- array(rep(0, length(tuning_c_vec) * n_subsample * (p - 1)), dim = c(length(tuning_c_vec), n_subsample, p - 1))
  for(kk in 1 : n_subsample){
    lambda_scale_tmp <- sqrt(log(p)/apply(N[kk, , ],1,min))
    for(j in 1 : length(tuning_c_vec)){
      out <- fit_nb(Sigma_hat[kk, , ], p, node_a, lambda_scale_tmp[-node_a] * tuning_c_vec[j], eta, tol)
      nbhd_est[j, kk, ] <- out$beta_hat != 0
    }
  }
  theta_stb <- apply(nbhd_est, c(1, 3), mean)
  stb_vec <- apply(2 * theta_stb * (1 - theta_stb), 1, mean)
  return(list(theta_stb = theta_stb, stb_vec = stb_vec))
}



calc_Sigma_PSD <- function(X, Obs, p, Entryest = FALSE){
  #calculate covariate estimate (entrywise and positive-semi definite projection) and pairwise sample size matrix,
  #using data X with observational pattern Obs
  Sigma_sum <- t(X) %*% X
  N <- t(Obs) %*% Obs
  Entryest_Sigma <- Sigma_sum / N
  Entryest_Sigma[N == 0] <- 0
  proj_out <- proj_cov(Entryest_Sigma, N, 0.01, 1, matrix(rep(0, p^2), p, p), matrix(rep(1, p^2), p, p), 0.001, 500)
  PSD_Sigma  <- proj_out$projected_S
  if(Entryest){
    return(list(PSD_Sigma = PSD_Sigma, N = N, Entryest_Sigma = Entryest_Sigma))
  }else{
    return(list(PSD_Sigma = PSD_Sigma, N = N))
  }
}

#----------------------------------------------------#
fit_nb <- function(Sigma, p, a, lambda, eta, tol,...){
  ## use projected gradient descent to solve neighborhood lasso
  # lambda: the p-1 dimensional penalty parameter
  # Sigma: the estimated p by p sample covariance matrix
  # a: target node
  # eta: intial step size parameter
  # ... can be initial beta (p-1 dimensional), if not supplied, beta_init = 0
  input <- list(...);
  eta_init=eta;
  if(length(input)>1){
    iter <- input[[2]]
    beta_init <- input[[1]]
  }else{
    iter <- 1000
    if(length(input)>0){
      beta_init <- input[[1]]
    }else{
      beta_init <- rep(0,p-1)
    }
  }
  beta <- beta_init
  kk=0;
  chg_per=Inf;
  A <- Sigma[-a,-a];
  b <- Sigma[-a,a]
  
  # initialize loss 
  loss=t(beta)%*%A%*%beta/2-sum(b*beta)+sum(abs(lambda*beta));
  
  while(kk<iter&chg_per>tol){
    grad <- A%*%beta - b
    # determine initial step size^(-1)
    if(kk > 0){
      denom <- sum((beta-beta_prev)*(grad-grad_prev));
      num <- sum((beta-beta_prev)^2);
      eta = num/denom;
    }
    accept <- F;
    while(!accept){
      grad_upd <- beta-grad*eta;
      # soft-thresholding
      beta_temp <- pmax(abs(grad_upd)-eta*lambda,0)*sign(grad_upd)
      loss_temp <- t(beta_temp)%*%A%*%beta_temp/2 - sum(b*beta_temp) + sum(abs(lambda*beta_temp));
      # check acceptance: non-increase in loss
      if(kk==0){
        accept <- (loss_temp <= loss)|eta<=0.001;
      }else{
        accept <- (loss_temp<=loss[kk+1])|eta<=0.001;
      }
      eta <- eta/2;
    }
    # check convergence
    beta_prev <- beta;
    grad_prev <- grad;
    beta <- beta_temp
    nm_prev <- norm(beta_prev,type='2')
    nm_diff <- norm(beta-beta_prev,type='2')
    if(nm_prev > 0){
      chg_per <- nm_diff/nm_prev;
    }else{
      chg_per <- nm_diff;
    }
    
    kk <- kk+1;
    loss <- c(loss,loss_temp);
  }
  if(chg_per <= tol){
    cvg <- T;
  }else{
    cvg=F;
  }
  return(list(beta_hat = beta,loss_path = loss,final_grad = grad,cvg = cvg))
}

find_var <- function(Sigma, N, Obs, beta, theta){
  p <- dim(Sigma)[[1]]
  var <- 0
  supp_theta <- which(theta != 0); supp_beta <- which(beta != 0)
  nonzeroObs <- which(N > 0)
  for(j1 in supp_theta){
    for(j2 in supp_beta){
      if(N[j1, j2] > 0){
        for(j3 in supp_theta){
            for(j4 in supp_beta){
              if(N[j3, j4] > 0){
              n_quad <- sum(Obs[, j1] * Obs[, j2] * Obs[, j3] * Obs[, j4])
              var <- var + (Sigma[j1, j3] * Sigma[j2, j4] + Sigma[j1, j4] * Sigma[j2, j3]) * 
                n_quad / N[j1, j2] / N[j3, j4] * theta[j1] * theta[j3] * beta[j2] * beta[j4]
            }
          }
        }
      }
    }
  }
  return(var)
}

#block <- list(1 : 20, 21 : 40, 41:60,61:80, 81:100, 101:120, 121:140,141:160, 161:180, 181:200); N4 <- list(list(list(list())))
#for(i1 in 1:length(block)){
#  for(i2 in 1:length(block)){
#    for(i3 in 1:length(block)){
#      for(i4 in 1:length(block)){
#        a1 <- array(Obs[, block[[i1]]], c(n_total, length(block[[i1]]), length(block[[i2]]), length(block[[i3]]), 
#                                          length(block[[i4]])))
#        a2 <- aperm(array(Obs[, block[[i2]]], c(n_total, length(block[[i2]]), length(block[[i1]]), length(block[[i3]]), 
#                                          length(block[[i4]]))), c(1, 3, 2, 4, 5))
#        a3 <- aperm(array(Obs[, block[[i3]]], c(n_total, length(block[[i3]]), length(block[[i1]]), length(block[[i2]]), 
#                                                length(block[[i4]]))), c(1, 3, 4, 2, 5))
#        a4 <- aperm(array(Obs[, block[[i4]]], c(n_total, length(block[[i4]]), length(block[[i1]]), length(block[[i2]]), 
#                                                length(block[[i3]]))), c(1, 3, 4, 5, 1))
#        N4[[i1]][[i2]][[i3]][[i4]] <- apply(a1 * a2 * a3 * a4, 1, sum)
#      }
#    }
#  }
#}

#find_var <- function(Sigma, N, N4, beta, theta){
#  cT <- outer(Sigma, Sigma)
#  cT <- aperm(cT, c(1,3,2,4)) + aperm(cT, c(1,3,4,2))
#  cTn <- cT * N4 / outer(N, N) 
#  var <- seq_ttm(as.tensor(cTn), list(theta, beta, theta, beta))@data[1]
#  return(var)
#}

m4_sum <- function(X){
  ## For a n*p matrix X, compute  a order-4 tensor with dimensions all equal to p. Each entry of the tensor is the summation of the 4d tuple data over n samples
  dim <- dim(X)
  n <- dim[1]
  p <- dim[2]
  m4_sum <- array(rep(0,p^4),c(p,p,p,p))
  if(n>p^4){
    for(j1 in 1:p){
      for(k1 in 1:p){
        for(j2 in 1:p){
          for(k2 in 1:p){
            m4_sum[j1,k1,j2,k2] <- sum(X[,j1]*X[,k1]*X[,j2]*X[k2])
          }
        }
      }
    }
  }
  else{
    for(i in 1:n){
      m4_sum <- m4_sum + outer(outer(outer(X[i,],X[i,]),X[i,]),X[i,])
    }
  }
  return(m4_sum)
}

seq_ttm <- function(tensor, vec_list){
  for(i in 1:length(vec_list)){
    tensor <- ttm(tensor, t(vec_list[[i]]), i) 
  }
  return(tensor)
}

proj_cov <- function(S, N, epsilon, mu, init_B, init_Lambda, tol, maxiter){
  ## projection of sample covariance matrix upon positive semi-definite cone w.r.t. weighted l_infty norm
  # S: sample covariance matrix (entry-wise estimation)
  # N: matrix of pair-wise sample sizes
  # init_B and init_Lambda need to be symmetric to ensure the symmetry of Sigma
  B <- init_B;
  Lambda <- init_Lambda;
  original_loss <- NULL; primal_loss <- NULL; Resid_sz <- NULL
  #dual_loss <- NULL
  for(i in 1:maxiter){
    if(i>1){
      Sigma_prev <- Sigma;
    }
    Sigma <- B + S + mu * Lambda;
    eig_out <- eigen(Sigma);
    Sigma <- eig_out[[2]]%*%diag(pmax(eig_out[[1]],epsilon))%*%t(eig_out[[2]])
    #Resid <- Sigma - B - S
    #inner_loss <- c(inner_loss, primal_loss - sum(Lambda * Resid) + 1 / 2 / mu * norm(Resid, "F"))
    A <- Sigma - mu * Lambda - S; B <- A;
    B[N > 0] <- A[N > 0] - proj_wt_l1_vec(A[N > 0], 1 / sqrt(N[N > 0]), mu / 2);
    primal_loss <- c(primal_loss, max(abs(B * sqrt(N)))/2);
    Resid <- Sigma - B - S
    Resid_sz <- c(Resid_sz, norm(Resid, "F"))
    #dual_loss <- inner_loss[inner_iter]
    Lambda_prev <- Lambda;
    Lambda <- Lambda - Resid/mu;
    original_loss <- c(original_loss,max(abs((Sigma - S)*sqrt(N))))
    if(i>1){
      chg_per <- norm(Sigma - Sigma_prev,"F")/norm(Sigma,"F")
      if(chg_per <= tol){
        break;
      }
    }
  }
  return(list(projected_S=Sigma,resid=B,original_loss=original_loss,Resid_sz = Resid_sz, primal_loss <- primal_loss,iteration_num=i,chg_per=chg_per))
}
#----------------------------------------------------#
proj_wt_l1_vec <- function(v,omega,r){
  ## projection of vector v on weighted l1 ball
  # omega: weight vector
  # r: radius
  if(sum(abs(omega*v))<=r){
    return(v)
  }else{
    z <- v/omega;
    sort_ind <- order(z,decreasing = TRUE)
    a <- cumsum((abs(omega*v))[sort_ind])
    b <- cumsum((omega^2)[sort_ind])
    j <- length(z);
    for(i in 1:length(z)){
      j_new <- sum(z[sort_ind]>(a[j]-r)/b[j])
      if(j_new == j){
        break;
      }else{
        j <- j_new;
      }
    }
    c <- (a[j]-r)/b[j];
    proj_v <- pmax(abs(v) - c*omega, 0 )*sign(v)
    return(proj_v)
  }
}

