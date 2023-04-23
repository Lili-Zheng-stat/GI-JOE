##First steps of our inference methods: debiasing without normalized by the variance

#----------------------------------------------------#
GIJOE_step1 <- function(X, Obs, n_total, p, thrs_c, signif_pval, onetuning = FALSE, symmetry = FALSE, ...){
  #given data X, Obs, n_total, p, screening threshold constant, 
  #find debiased test statistics for the whole graph, without normalized by the variance estimates
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
    print("Only 1 node left!")
    return(list(success = FALSE, p_screened = p_screened))
  }else{
    X_screened <- X[, ind_order[1:p_screened]]
    Obs_screened <- Obs[, ind_order[1:p_screened]]
    
    #testing
    if(length(input)>0){
      test_out <- GIJOE_testing_step1(X_screened, Obs_screened, n_total, p_screened, input[[1]])
    }else if(!onetuning){
      test_out <- GIJOE_testing_step1(X_screened, Obs_screened, n_total, p_screened)
    }else{
      test_out <- GIJOE_testing_step1_onetuning(X_screened, Obs_screened, n_total, p_screened, symmetry)
    }
    return(list(success = TRUE, test_out = test_out, kept_nodes = ind_order[1:p_screened]))
  }
}

#----------------------------------------------------#
GIJOE_testing_step1 <- function(X, Obs, n_total, p, ...){
  ## Given data X, observational pattern Obs, pairwise sample sizes N, 
  ## output the debiased test statistics for each edge (neighborhood lasso on the node with smaller index; 
  ## if wanting neighborhood lasso for the other node, can use a permutation on rows of X before input).
  # tuning_c can be an additional input; otherwise, chosen by stability
  input <- list(...);
  
  #find entrywise estimate of Sigma over full data and random subsamples; also find pairwise sample sizes for subsamples
  Sigma_hat <- array(rep(0, p^2*21), c(21, p, p)); N <- array(rep(0, p^2*21), c(21, p, p));
  Sigma_out <- calc_Sigma_PSD(X, Obs, p, Entryest = TRUE); 
  Sigma_hat[1, , ] <- Sigma_out$PSD_Sigma; Entryest_Sigma <- Sigma_out$Entryest_Sigma
  N[1, , ] <- Sigma_out$N
  for(i in 2 : 21){
    set.seed(i - 1)
    selected_ind <- as.logical(rbinom(n_total, 1, 0.8))
    Sigma_out <- calc_Sigma_PSD(X[selected_ind, ], Obs[selected_ind, ], p)
    Sigma_hat[i, , ] <- Sigma_out$PSD_Sigma; N[i, , ] <- Sigma_out$N
  }

  #loop over all pairs (a, b): neighborhood lasso for a and debiasing vector for (a,b)
  eta <- 1;tol <- 0.00001; lambda_scale <- sqrt(log(p)/apply(N[1, , ], 1, min))
  beta_hat <- matrix(rep(0, p^2), p, p); theta_hat_debiasing <- array(rep(0, p^3), c(p, p, p)); 
  theta_hat_varest <- array(rep(0, p^3), c(p, p, p)); beta_tilde <- matrix(rep(0, p^2), p, p)
  var_est <- matrix(rep(0, p^2), p, p); tuning_c_vec <- rep(0, p); tuning_c_mat <- matrix(rep(0, p^2), p, p); test <- matrix(rep(0, p^2), p, p)
  for(node_a in 1 : (p - 1)){
    if(length(input) > 0){
      fit_out1 <- fit_nb(Sigma_hat[1, , ], node_a, lambda_scale[-node_a] * input[[1]], eta, tol)
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
        fit_out2 <- fit_nb(Sigma_hat[1, -node_a, -node_a], node_b_ind, lambda_scale_tmp[-node_b_ind] * input[[1]], eta, tol)
      }else{
        fit_out2 <- est_nb_tuning(Sigma_hat[, -node_a, -node_a], N[, -node_a, -node_a], p - 1, stb_thrs = 0.05, node_b_ind, lambda_scale_tmp[-node_b_ind], eta, tol)
        tuning_c_mat[node_a, node_b] <- fit_out2$tuning_c
      }
      theta_hat_varest[node_a, node_b, -c(node_a, node_b)] <- - fit_out2$beta_hat
      theta_hat_varest[node_a, node_b, node_b] <- 1
      tau_varest <- as.numeric(1 / (t(theta_hat_varest[node_a, node_b, ]) %*% Sigma_hat[1, , ] %*% theta_hat_varest[node_a, node_b, ]))
      tau_debiasing <- as.numeric(1 / (Sigma_hat[1, node_b, ] %*% theta_hat_varest[node_a, node_b, ]))
      theta_hat_debiasing[node_a, node_b, ] <- theta_hat_varest[node_a, node_b,] * tau_debiasing
      theta_hat_varest[node_a, node_b, ] <- theta_hat_varest[node_a, node_b,] * tau_varest
      #debiased neighborhood lasso (node_a, node_b)
      beta_tilde[node_a, node_b] = as.numeric(beta_hat[node_a, node_b] - t(theta_hat_debiasing[node_a, node_b, ]) %*% (Entryest_Sigma %*% beta_hat[node_a, ] - Entryest_Sigma[, node_a]))
    }
  }
  return(list(beta_hat = beta_hat, theta_hat_debiasing = theta_hat_debiasing, theta_hat_varest = theta_hat_varest, beta_tilde = beta_tilde, N = N[1, , ], tuning_c1 = tuning_c_vec, tuning_c2 = tuning_c_mat, 
              Entryest_Sigma = EntryestSigma, PSD_Sigma = Sigma_hat[1, , ]))
}

#----------------------------------------------------#
GIJOE_testing_step1_onetuning <- function(X, Obs, n_total, p, symmetry){
  ## Choose one fixed tuning_c for all pairs of node; 
  ## if symmetry = true, for each pair of nodes, perform two tests in a symmetric way;
  ## otherwise, perform one test based on the order
  ## Given data X, observational pattern Obs, pairwise sample sizes N, 
  ## output the debiased test statistics for each edge (neighborhood lasso on the node with smaller index; 
  ## if wanting neighborhood lasso for the other node, can use a permutation on rows of X before input).
  
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
  start_t <- Sys.time()
  eta <- 1;tol <- 0.00001; lambda_scale <- sqrt(log(p)/apply(N[1, , ], 1, min))
  tuning_fitting_results <- est_graph_tuning(Sigma_hat, N, p, stb_thrs = 0.05, lambda_scale, eta, tol)
  beta_hat <- tuning_fitting_results$nblasso_combined
  tuning_c <- tuning_fitting_results$tuning_c
  print(Sys.time() - start_t)
  print("Finished tuning!")
  
  #loop over all pairs (a, b): debiasing vector for (a,b); 
  #if symmetry = true, (a, b) and (b, a) are both considered; otherwise, consider (a, b) for a < b
  theta_hat_debiasing <- array(rep(0, p^3), c(p, p, p)); 
  theta_hat_varest <- array(rep(0, p^3), c(p, p, p)); 
  beta_tilde <- matrix(rep(0, p^2), p, p)
  var_est <- matrix(rep(0, p^2), p, p); test <- matrix(rep(0, p^2), p, p)
  start_t <- Sys.time()
  for(node_a in 1 : (p - 1)){
    ind <- (1 : p)[-node_a]
    for(node_b in ind){
      if(node_b > node_a || symmetry){
        node_b_ind <- which(ind == node_b)
        lambda_scale_tmp <- sqrt(log(p - 1)/apply(N[1, -node_a, -node_a], 1, min))
        fit_out2 <- fit_nb(Sigma_hat[1, -node_a, -node_a], node_b_ind, lambda_scale_tmp[-node_b_ind] * tuning_c, eta, tol)
        theta_hat_varest[node_a, node_b, -c(node_a, node_b)] <- - fit_out2$beta_hat
        theta_hat_varest[node_a, node_b, node_b] <- 1
        tau_varest <- as.numeric(1 / (t(theta_hat_varest[node_a, node_b, ]) %*% Sigma_hat[1, , ] %*% theta_hat_varest[node_a, node_b, ]))
        tau_debiasing <- as.numeric(1 / (Sigma_hat[1, node_b, ] %*% theta_hat_varest[node_a, node_b, ]))
        theta_hat_debiasing[node_a, node_b, ] <- theta_hat_varest[node_a, node_b,] * tau_debiasing
        theta_hat_varest[node_a, node_b, ] <- theta_hat_varest[node_a, node_b,] * tau_varest
        #debiased neighborhood lasso (node_a, node_b)
        beta_tilde[node_a, node_b] = as.numeric(beta_hat[node_a, node_b] - t(theta_hat_debiasing[node_a, node_b, ]) %*% (Entryest_Sigma %*% beta_hat[node_a, ] - Entryest_Sigma[, node_a]))
      }
    }
  }
  print(Sys.time() - start_t)
  print("Graph estimation using neighborhood lasso is done!")
  return(list(beta_hat = beta_hat, theta_hat_debiasing = theta_hat_debiasing, theta_hat_varest = theta_hat_varest, beta_tilde = beta_tilde, N = N[1, , ], tuning_c = tuning_c, Entryest_Sigma = Entryest_Sigma, PSD_Sigma = Sigma_hat[1, , ]))
}

find_variances_GIJOE <- function(nblasso_step1_results, PSD_Sigma, theta_hat, Obs, row_ind, col_ind){
  beta_hat <- nblasso_step1_results$test_out$beta_hat;
  N <- nblasso_step1_results$test_out$N
  var_est <- rep(0, length(row_ind))
  start_t <- Sys.time()
  for(i in 1 : length(row_ind)){
    a <- row_ind[i] 
    b <- col_ind[i]
    if(length(a) > 0){
      beta_hat_tmp <- -beta_hat[a, ]; beta_hat_tmp[a] <- 1; beta_hat_tmp[abs(beta_hat_tmp) / max(abs(beta_hat_tmp)) <= 10^(-5)] <- 0;
      theta_tmp <- theta_hat[a, b, ]; theta_tmp[abs(theta_tmp) / max(abs(theta_tmp)) <= 10^(-5)] <- 0;
      var_est[i] <- find_var(PSD_Sigma, N, Obs[, nblasso_step1_results$kept_nodes], beta_hat_tmp, theta_tmp)
    }else{
      var_est[i] <- NA
    }
    if(i == 20 * floor(i/20)){
      print(Sys.time() - start_t)
      print(sprintf("finished variance estimation for %d edges", i))
    }
  }
  return(var_est)
}

find_variances_GIJOE_sameN <- function(nblasso_step1_results, n, a){
  # When the sample sizes are all equal to n, given nblasso step1 results and node a, find variance estimate for node pair (a,b) for b>a.
  beta_hat <- nblasso_step1_results$test_out$beta_hat;
  p <- dim(beta_hat)[1]
  theta_hat <- nblasso_step1_results$test_out$theta_hat
  Entryest_Sigma <- nblasso_step1_results$test_out$Entryest_Sigma
  beta_tmp <- -beta_hat[a, ]; beta_tmp[a] <- 1; 
  theta_tmp <- theta_hat[a, , ];
  var1 <- diag(theta_tmp %*% Entryest_Sigma %*% t(theta_tmp)) * as.numeric(t(beta_tmp) %*% Entryest_Sigma %*% beta_tmp) / n
  var2 <- (theta_tmp %*% Entryest_Sigma %*% beta_tmp)^2 / n
  return(var1[(a + 1) : p] + var2[(a + 1) : p])
}


GIJOE_thrs_step3 <- function(step1_results, var_est_results, p, signif_pval, testing_thrs = 0){
  #given step 1 results using nblasso and variance estimate, find edges with weights significantly larger than testing_thrs (default is 0)
  var_est <- matrix(rep(0, p^2), p, p); 
  p_kept <- length(step1_results$kept_nodes);
  for(a in 1 : (p_kept - 1)){
    for(b in (a + 1) : p_kept){
      if(is.list(var_est_results)){
        var_est[step1_results$kept_nodes[a], step1_results$kept_nodes[b]] <- var_est_results$var_est[var_est_results$rows == step1_results$kept_nodes[a] & var_est_results$cols == step1_results$kept_nodes[b]]
      }else{
        var_est[step1_results$kept_nodes[a], step1_results$kept_nodes[b]] <- var_est_results[a, b]
      }
    }
  }
  var_est <- var_est + t(var_est);
  var_est_kept_nodes <- var_est[step1_results$kept_nodes, step1_results$kept_nodes]
  beta_tilde <- step1_results$test_out$beta_tilde;
  
  test <- matrix(rep(0, p^2), p, p)
  test[step1_results$kept_nodes, step1_results$kept_nodes][upper.tri(beta_tilde)] <- (abs(beta_tilde[upper.tri(beta_tilde)]) - testing_thrs)/sqrt(var_est_kept_nodes[upper.tri(var_est_kept_nodes)]); 
  test <- test + t(test)
  
  #find p-values and significant graphs
  p_val <- matrix(pmin(1,(1-pnorm(test))*2), p, p)
  #Bonferroni correction
  p_val_Bonferroni <- pmin(p_val*p_kept*(p_kept-1)/2, 1)
  signif_graph_Bonferroni <- p_val_Bonferroni < signif_pval
  
  #Holm-Bonferroni correction
  ind <- matrix(1 : (p_kept ^ 2), p_kept, p_kept)
  p_val_kept_nodes <- p_val[step1_results$kept_nodes, step1_results$kept_nodes]
  p_val_vec <- p_val_kept_nodes[upper.tri(p_val_kept_nodes)]; order_p_val <-order(p_val_vec)
  k <- min(which(p_val_vec[order_p_val] > signif_pval / seq(from = p_kept*(p_kept-1)/2, to = 1, by = -1)))
  signif_ind_holm <- ind[upper.tri(ind)][order_p_val[1 : (k - 1)]]
  signif_graph_holm <- matrix(rep(0, p^2), p, p)
  signif_graph_holm[step1_results$kept_nodes, step1_results$kept_nodes][signif_ind_holm] <- 1; signif_graph_holm <- signif_graph_holm + t(signif_graph_holm)
  
  #FDR control
  k <- max(which(p_val_vec[order_p_val] <= signif_pval * (1 : (p_kept * (p_kept - 1) / 2)) / (p_kept * (p_kept - 1) / 2)))
  rho <- signif_pval * k / (p_kept * (p_kept - 1) / 2)
  t_max <- sqrt(2 * log(p_kept * (p_kept - 1) / 2) - 2 * log(log(p_kept * (p_kept - 1) / 2))); rho_min <- 2 * (1 - pnorm(t_max))
  if(rho < rho_min){
    rho <- 2 * (1 - pnorm(sqrt(2 * log(p_kept * (p_kept - 1) / 2))))
  }
  signif_ind_FDR <- ind[upper.tri(ind)][which(p_val_vec <= rho)]
  signif_graph_FDR <- matrix(rep(0, p^2), p, p)
  signif_graph_FDR[step1_results$kept_nodes, step1_results$kept_nodes][signif_ind_FDR] <- 1; signif_graph_FDR <- signif_graph_FDR + t(signif_graph_FDR)
  
  #estimated graph 
  graph_est <- matrix(rep(0, p^2), p, p)
  graph_est[step1_results$kept_nodes, step1_results$kept_nodes] <- step1_results$test_out$beta_hat != 0;
  graph_est_OR <- pmax(graph_est, t(graph_est))
  graph_est_AND <- pmin(graph_est, t(graph_est))
  
  return(list(signif_graph_holm = signif_graph_holm, signif_graph_FDR = signif_graph_FDR, graph_est_OR = graph_est_OR, graph_est_AND = graph_est_AND, p_val = p_val,
              p_val_Bonferroni = p_val_Bonferroni))
}



GIJOE_step3 <- function(step1_results, var_est_results, p, signif_pval){
  #given step 1 results using nblasso and variance estimate, find significant graphs
  var_est <- matrix(rep(0, p^2), p, p); 
  p_kept <- length(step1_results$kept_nodes);
  for(a in 1 : (p_kept - 1)){
    for(b in (a + 1) : p_kept){
      var_est[step1_results$kept_nodes[a], step1_results$kept_nodes[b]] <- var_est_results$var_est[var_est_results$rows == step1_results$kept_nodes[a] & var_est_results$cols == step1_results$kept_nodes[b]]
    }
  }
  var_est <- var_est + t(var_est);
  var_est_kept_nodes <- var_est[step1_results$kept_nodes, step1_results$kept_nodes]
  beta_tilde <- step1_results$test_out$beta_tilde;
  
  test <- matrix(rep(0, p^2), p, p)
  test[step1_results$kept_nodes, step1_results$kept_nodes][upper.tri(beta_tilde)] <- beta_tilde[upper.tri(beta_tilde)]/sqrt(var_est_kept_nodes[upper.tri(var_est_kept_nodes)]); 
  test <- test + t(test)
  
  #find p-values and significant graphs
  p_val <- matrix(pmin(1,(1-pnorm(abs(test)))*2), p, p)
  #Bonferroni correction
  p_val_Bonferroni <- pmin(p_val*p_kept*(p_kept-1)/2, 1)
  signif_graph_Bonferroni <- p_val_Bonferroni < signif_pval
  
  #Holm-Bonferroni correction
  ind <- matrix(1 : (p_kept ^ 2), p_kept, p_kept)
  p_val_kept_nodes <- p_val[step1_results$kept_nodes, step1_results$kept_nodes]
  p_val_vec <- p_val_kept_nodes[upper.tri(p_val_kept_nodes)]; order_p_val <-order(p_val_vec)
  ind_tmp <- which(p_val_vec[order_p_val] > signif_pval / seq(from = p_kept*(p_kept-1)/2, to = 1, by = -1))
  if(length(ind_tmp) > 0){
    k <- min(ind_tmp)
  }else{
    k <- length(p_val_vec) + 1
  }
  if(k > 1){
    signif_ind_holm <- ind[upper.tri(ind)][order_p_val[1 : (k - 1)]]
  }else{
    signif_ind_holm <- NULL
  }
  signif_graph_holm <- matrix(rep(0, p^2), p, p)
  signif_graph_holm[step1_results$kept_nodes, step1_results$kept_nodes][signif_ind_holm] <- 1; signif_graph_holm <- signif_graph_holm + t(signif_graph_holm)
  
  #FDR control
  ind_tmp <- which(p_val_vec[order_p_val] <= signif_pval * (1 : (p_kept * (p_kept - 1) / 2)) / (p_kept * (p_kept - 1) / 2))
  if(length(ind_tmp) > 0){
    k <- max(ind_tmp)
    rho <- signif_pval * k / (p_kept * (p_kept - 1) / 2)
  }else{
    rho <- -Inf
  }
  t_max <- sqrt(2 * log(p_kept * (p_kept - 1) / 2) - 2 * log(log(p_kept * (p_kept - 1) / 2))); rho_min <- 2 * (1 - pnorm(t_max))
  if(rho < rho_min){
    rho <- 2 * (1 - pnorm(sqrt(2 * log(p_kept * (p_kept - 1) / 2))))
  }
  signif_ind_FDR <- ind[upper.tri(ind)][which(p_val_vec <= rho)]
  signif_graph_FDR <- matrix(rep(0, p^2), p, p)
  signif_graph_FDR[step1_results$kept_nodes, step1_results$kept_nodes][signif_ind_FDR] <- 1; signif_graph_FDR <- signif_graph_FDR + t(signif_graph_FDR)
  
  #estimated graph 
  graph_est <- matrix(rep(0, p^2), p, p)
  graph_est[step1_results$kept_nodes, step1_results$kept_nodes] <- step1_results$test_out$beta_hat != 0;
  graph_est_OR <- pmax(graph_est, t(graph_est))
  graph_est_AND <- pmin(graph_est, t(graph_est))
  
  return(list(signif_graph_Bonferroni = signif_graph_Bonferroni, signif_graph_holm = signif_graph_holm, 
              signif_graph_FDR = signif_graph_FDR, graph_est_AND = graph_est_AND, graph_est_OR = graph_est_OR, p_val = p_val,
              p_val_Bonferroni = p_val_Bonferroni))
}


p_val_to_Holm_FDR <- function(p_val_mat, signif_pval = 0.05, p_kept, p){
  #p_val_mat is the p-value matrix, without reordering the nodes. 
  #The node ordering of the returned graph selection result is the same as the p_val_vec.
  #This function can be applied to debiased graphical lasso (baseline) results
  
  p_val_vec <- p_val_mat[upper.tri(p_val_mat)]
  order_p_val <-order(p_val_vec)
  #Holm-Bonferroni correction
  ind <- matrix(1 : (p ^ 2), p, p)
  k <- min(which(p_val_vec[order_p_val[1 : (p_kept * (p_kept - 1) / 2)]] > signif_pval / seq(from = p_kept * (p_kept - 1) / 2, to = 1, by = -1)))
  if(k > 1){
    signif_ind_holm <- ind[upper.tri(ind)][order_p_val[1 : (k - 1)]]
  }else{
    signif_ind_holm <- NULL
  }
  signif_graph_holm <- matrix(rep(0, p^2), p, p)
  signif_graph_holm[signif_ind_holm] <- 1; signif_graph_holm <- signif_graph_holm + t(signif_graph_holm)
  
  #FDR control
  ind_tmp <- which(p_val_vec[order_p_val[1 : (p_kept * (p_kept - 1) / 2)]] <= signif_pval * (1 : (p_kept * (p_kept - 1) / 2)) / (p_kept * (p_kept - 1) / 2))
  if(length(ind_tmp) > 0){
    k <- max(ind_tmp)
    rho <- signif_pval * k / (p_kept * (p_kept - 1) / 2)
  }else{
    rho <- -Inf
  }
  t_max <- sqrt(2 * log(p_kept * (p_kept - 1) / 2) - 2 * log(log(p_kept * (p_kept - 1) / 2))); rho_min <- 2 * (1 - pnorm(t_max))
  if(rho < rho_min){
    rho <- 2 * (1 - pnorm(sqrt(2 * log(p_kept * (p_kept - 1) / 2))))
  }
  signif_ind_FDR <- ind[upper.tri(ind)][which(p_val_vec <= rho)]
  signif_graph_FDR <- matrix(rep(0, p^2), p, p)
  signif_graph_FDR[signif_ind_FDR] <- 1; signif_graph_FDR <- signif_graph_FDR + t(signif_graph_FDR)
  return(list(signif_graph_holm = signif_graph_holm, signif_graph_FDR = signif_graph_FDR))
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
  out <- fit_nb(Sigma_hat[1, , ], node_a, lambda_scale[-node_a] * tuning_c_max, eta, tol)
  while(sum(out$beta_hat != 0) > 0 || node_a < p){
    if(sum(out$beta_hat != 0) > 0){
      tuning_c_max <- tuning_c_max * 2
    }else{
      node_a <- node_a + 1;
    }
    out <- fit_nb(Sigma_hat[1, , ], node_a, lambda_scale[-node_a] * tuning_c_max, eta, tol)
  }
  if(tuning_c_max == 1){
    while(sum(out$beta_hat != 0) == 0){
      if(node_a < p){
        node_a <- node_a + 1
      }else{
        tuning_c_max <- tuning_c_max / 2
        node_a <- 1
      }
      out <- fit_nb(Sigma_hat[1, , ], node_a, lambda_scale[-node_a] * tuning_c_max, eta, tol)
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
    out_final[[node_a]] <- fit_nb(Sigma_hat[1, , ], node_a, lambda_scale[-node_a] * tuning_c_vec[choice_ind], eta, tol)
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
        out <- fit_nb(Sigma_hat[kk, , ], node_a, lambda_scale_tmp[-node_a] * tuning_c_vec[j], eta, tol)
        graph_est_arr[j, kk, node_a, -node_a] <- out$beta_hat != 0
      }
    }
  }
  theta_stb <- apply(graph_est_arr, c(1, 3, 4), mean)
  stb_vec <- apply(2 * theta_stb * (1 - theta_stb), 1, mean)
  return(list(theta_stb = theta_stb, stb_vec = stb_vec))
}
#----------------------------------------------------#
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

