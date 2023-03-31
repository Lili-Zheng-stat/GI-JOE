## functions for simulations
#model generation
#----------------------------------------------------#
precision_generation <- function(p,edge_val_range,eigmin,seed,g_type,target_nodes = NULL, signal = NULL, ...){
  input <- list(...)
  set.seed(seed);
  Adj_mat <- matrix(rep(0,p^2),p,p); Theta <- matrix(rep(0,p^2),p,p)
  if(g_type == "banded"){
    for(i in 1 : p){
      nb <- i + c(-3 : -1, 1 : 3); nb <- nb[nb >= 1 & nb <= p]
      Adj_mat[i, nb] <- 1
      Theta[i, nb] <- 0.8^(abs(nb - i))
    }
  }else{
    if(g_type == "chain"){
      for(i in 1:(p-1)){
        Adj_mat[i,i+1]=1
      }
      Adj_mat <- Adj_mat+t(Adj_mat)
    }
    else if(g_type == "block"){
      K <- list(...)[[1]]
      for(i in 1:(p/K)){
        Adj_mat[(K*(i-1)+1):(K*i),(K*(i-1)+1):(K*i)] <- 1
      }
      diag(Adj_mat) <- 0
    }
    else if(g_type == "star"){
      if(length(input) == 0){
        Adj_mat[1, 2 : p] <- 1
      }
      else{
        nstar <- input[[1]]
        nnodes_perstar <- floor(p / nstar)
        for(i in 1 : (nstar - 1)){
          Adj_mat[nnodes_perstar * (i - 1) + 1, (nnodes_perstar * (i - 1) + 2) : (nnodes_perstar * i)] <- 1
        }
        Adj_mat[nnodes_perstar * (nstar - 1) + 1, (nnodes_perstar * (nstar - 1) + 2) : p] <- 1
      }
      Adj_mat <- Adj_mat + t(Adj_mat)
    }else if(g_type == "ER"){
      if(length(input) > 0){
        connected_prob <- input[[1]]
      }else{
        connected_prob <- 0.1
      }
      Adj_mat[lower.tri(Adj_mat)] <- rbinom(p * (p - 1) / 2, 1, connected_prob)
      Adj_mat <- Adj_mat + t(Adj_mat)
    }else if(g_type == "spatial"){
      r <- list(...)[[1]] #radius
      nodes_loc <- matrix(runif(p * 2, 0, 1), p, 2)
      Adj_mat <- as.matrix(rdist(nodes_loc, metric = "euclidean")) < r
    }else if(g_type == "circle"){
      for(i in 1:(p-1)){
        Adj_mat[i,i+1]=1
      }
      Adj_mat[p, 1] <- 1
      Adj_mat <- Adj_mat+t(Adj_mat)
    }else if(g_type == "small-world"){
      if(length(input) > 0){
        degree <- input[[1]]
        if(length(input) > 1){
          rewiring_prob <- input[[2]]
        }else{
          rewiring_prob <- 0.5
        }
      }else{
        degree <- floor(0.1 * p)
        rewriring_prob <- 0.5
      }
      g <- sample_smallworld(1, p, floor(degree / 2), rewiring_prob, loops = FALSE, multiple = FALSE)
      Adj_mat <- as.matrix(get.adjacency(g))
    }else if(g_type == "scale-free"){
      if(length(input) > 0){
        g <- sample_pa(n = p, m = input[[1]], directed = FALSE)
      }else{
        g <- sample_pa(n = p, m = floor(0.1 * p), directed = FALSE)
      }
      Adj_mat <- as.matrix(get.adjacency(g))
    }
    for(i in 1:(p-1)){
      for(j in (i+1):p){
        if(Adj_mat[i,j]==1){
          Theta[i,j] <- runif(1,edge_val_range[1],edge_val_range[2])
        }
      }
    }
    Theta <- Theta + t(Theta)
  }
  if(!is.null(target_nodes) & !is.null(signal)){
    Theta[target_nodes[1], target_nodes[2]] <- signal
    Theta[target_nodes[2], target_nodes[1]] <- signal
  }
  Theta <- Theta + (eigmin-min(eigen(Theta)[[1]]))*diag(p)
  return(list(Theta,Adj_mat))
}
#----------------------------------------------------#
#data generation
data_generation_pairwise <- function(Sigma, N, seed){
  ## generate data, observation pattern and sample covariance from pairwise measurements from N(0,Sigma), each pair measured N[i,j] times
  p <- dim(Sigma)[1]
  ind <- upper.tri(N, diag = FALSE)
  total_n <- sum(N[ind]);
  set.seed(seed)
  X <- matrix(rep(0, total_n * p), total_n, p);
  Obs <- matrix(rep(0, total_n * p), total_n, p);
  start_n <- 1;
  diag_n <- rep(0, p); diag_Sigma <- rep(0, p);
  Sigma_hat <- matrix(rep(0, p^2), p, p)
  Ind <- list(); 
  for(i in 1:(p-1)){
    Ind[[i]] <- list(); 
    for(j in (i+1):p){
      end_n <- start_n + N[i,j] - 1;
      Obs[start_n:end_n, c(i,j)] <- 1;
      sqrt_Sigma_tmp <- msqrt(Sigma[c(i,j), c(i,j)])[[1]];
      X[start_n : end_n, c(i, j)] <- matrix(rnorm(N[i, j] * 2, 0, 1), N[i, j], 2) %*% sqrt_Sigma_tmp; 
      Ind[[i]][[j]] <- start_n : end_n #Ind only includes the sample sizes for i < j
      Sigma_hat[i,j] <- sum(X[start_n : end_n, i] * X[start_n : end_n, j])/N[i, j];
      diag_Sigma[i] <- diag_Sigma[i] + sum(X[start_n : end_n, i]^2);
      diag_n[i] <- diag_n[i] + N[i, j];
      diag_Sigma[j] <- diag_Sigma[j] + sum(X[start_n : end_n,j]^2);
      diag_n[j] <- diag_n[j] + N[i, j];
      start_n <- end_n + 1;
    }
  }
  Sigma_hat <- Sigma_hat + t(Sigma_hat);
  diag(Sigma_hat) <- diag_Sigma / diag_n; 
  X <- X * Obs;
  return(list(X = X, Obs = Obs, Sigma_hat = Sigma_hat, Ind = Ind))
}
#output cross validation indices when generating data
data_generation_pairwise_cv_ind <- function(Sigma, N, N_cv_test, seed){
  ## generate data, observation pattern and sample covariance from pairwise measurements from N(0,Sigma), each pair measured N[i,j] times
  p <- dim(Sigma)[1]
  ind <- upper.tri(N, diag = FALSE)
  total_n <- sum(N[ind]);
  set.seed(seed)
  X <- matrix(rep(0, total_n * p), total_n, p);
  Obs <- matrix(rep(0, total_n * p), total_n, p);
  start_n <- 1;
  diag_n <- rep(0, p); diag_Sigma <- rep(0, p);
  Sigma_hat <- matrix(rep(0, p^2), p, p)
  Ind <- list(); Ind_cv_test <- list(); Ind_cv_train <- list(); 
  for(i in 1:(p-1)){
    Ind[[i]] <- list(); Ind_cv_test[[i]] <- list(); Ind_cv_train[[i]] <- list();
    for(j in (i+1):p){
      end_n <- start_n + N[i,j] - 1;
      Obs[start_n:end_n, c(i,j)] <- 1;
      sqrt_Sigma_tmp <- msqrt(Sigma[c(i,j), c(i,j)])[[1]];
      X[start_n : end_n, c(i, j)] <- matrix(rnorm(N[i, j] * 2, 0, 1), N[i, j], 2) %*% sqrt_Sigma_tmp; 
      Ind[[i]][[j]] <- start_n : end_n #Ind only includes the sample sizes for i < j
      Ind_cv_test[[i]][[j]] <- list(); Ind_cv_train[[i]][[j]] <- list();
      for(k in 1:5){
        Ind_cv_test[[i]][[j]][[k]] <- (start_n + (k - 1) * N_cv_test[i, j]) : (start_n + k * N_cv_test[i, j] - 1)
        Ind_cv_train[[i]][[j]][[k]] <- Ind[[i]][[j]][-Ind_cv_test[[i]][[j]][[k]]]
      }
      Sigma_hat[i,j] <- sum(X[start_n : end_n, i] * X[start_n : end_n, j])/N[i, j];
      diag_Sigma[i] <- diag_Sigma[i] + sum(X[start_n : end_n, i]^2);
      diag_n[i] <- diag_n[i] + N[i, j];
      diag_Sigma[j] <- diag_Sigma[j] + sum(X[start_n : end_n,j]^2);
      diag_n[j] <- diag_n[j] + N[i, j];
      start_n <- end_n + 1;
    }
  }
  Sigma_hat <- Sigma_hat + t(Sigma_hat);
  diag(Sigma_hat) <- diag_Sigma / diag_n; 
  X <- X * Obs;
  return(list(X = X, Obs = Obs, Sigma_hat = Sigma_hat, Ind = Ind, Ind_cv_test = Ind_cv_test, Ind_cv_train = Ind_cv_train))
}
#----------------------------------------------------#
#find key neighborhood for assigning pairwise sample sizes
find_key_ind <- function(Adj_mat, Theta, p, target_nodes){
  nb_a1 <- which(Adj_mat[target_nodes[1],] != 0)
  nb_a2 <- c(target_nodes[1], nb_a1)
  Theta_a <- matrix(rep(0, p^2), p, p);
  Theta_a[-target_nodes[1], -target_nodes[1]] <- 
    Theta[-target_nodes[1], -target_nodes[1]] - Theta[-target_nodes[1], target_nodes[1]] %*% t(Theta[target_nodes[1], -target_nodes[1]])/Theta[target_nodes[1],target_nodes[1]]
  nb_a_b2 <- which(Theta_a[target_nodes[2],] != 0)
  pair_ind1 <- unique(c(nb_a1, nb_a_b2)) #n1 is the minimum sample size for i in pair_ind1 and arbitrary j
  pair_ind2_1 <- nb_a2; pair_ind2_2 <- nb_a_b2 #n2 is the minimum sample size for i in pair_ind2_1 and j in pair_ind2_2
  return(list(Theta_a = Theta_a, pair_ind1 = pair_ind1, pair_ind2 = list(pair_ind2_1 = pair_ind2_1, pair_ind2_2 = pair_ind2_2)))
}
#----------------------------------------------------#
#perform certain experiment for one replicate; 
test_pairwise_experiment <- function(p, target_nodes, n1, n2, seed_data, tuning_c = NULL, ...){
  ##given dimension p, graph type or a particular graph (Theta, Sigma, Adj_mat; or graph_type, seed_model), 
  # generate pairwise measurements data (setting seed_data) with sample sizes n1, n2, and other N[i,j];
  # n1 is scalar, n2 can be scalar or matrix of size (d_a+1)*(d_b^{(a)}+1).
  # then using the data to test whether there is an edge between target_nodes[1] and target_nodes[2].
  input <- list(...)
  #generate model
  if(sum(is.na(input$Theta))||sum(is.na(input$Sigma))||sum(is.na(input$Adj_mat))){
      if(is.na(input$graph_type)||is.na(input$model_seed)){
        print("Graph type and model_seed or graph parameters are missing!")
      }else{
        edge_val_range <- c(0.6,0.8)
        eigmin <- 0.25
        if(input$graph_type == "block"){
          model_param <- precision_generation(p,edge_val_range,eigmin,input$model_seed,input$graph_type,input$blocksize) 
        }else{
          model_param <- precision_generation(p,edge_val_range,eigmin,input$model_seed,input$graph_type) 
        }
        Theta <- model_param[[1]]
        Adj_mat <- model_param[[2]]
        diag(Adj_mat) <- 0
        Sigma <- solve(Theta)
      }
  }else{
    Theta <- input$Theta
    Sigma <- input$Sigma
    Adj_mat <- input$Adj_mat 
  }
  
  #setting sample sizes
  N <- matrix(rep(50, p^2), p, p)
  key_ind_list = find_key_ind(Adj_mat, Theta, p, target_nodes)
  N[key_ind_list$pair_ind1,] = n1; N[, key_ind_list$pair_ind1] = n1
  N[key_ind_list$pair_ind2$pair_ind2_1,key_ind_list$pair_ind2$pair_ind2_2] = n2; N[key_ind_list$pair_ind2$pair_ind2_2,key_ind_list$pair_ind2$pair_ind2_1] = t(n2)
  N_cv_train <- floor(N / 5 * 4); N_cv_test <- N - N_cv_train;
  #generate data
  data <- data_generation_pairwise(Sigma, N, seed_data)
  
  #test
  #update diagonal elements of N
  diag(N) <- apply(N, 1, sum); diag(N_cv_train) <- apply(N_cv_train, 1, sum); diag(N_cv_test) <- apply(N_cv_test, 1, sum);
  
  #N4 <- array(rep(0,p^4),rep(p,4))
  #for(j in 1:(p-1)){
  #  for(k in ((j+1):p)){
  #   N4[c(j,k),c(j,k),c(j,k),c(j,k)] <- N[j,k]
  #  }
  #}
  #for(j in 1:p){
  #  N4[j,j,j,j] <- N[j,j]
  #}
  if(is.null(tuning_c)){
    test_out <- edge_testing_saveSigma(data$X, dim(data$X)[1], data$Ind, data$Sigma_hat, 
                             p, N, target_nodes[1], target_nodes[2])
  }else{
    test_out <- edge_testing_saveSigma(data$X, dim(data$X)[1], data$Ind, data$Sigma_hat, 
                             p, N, target_nodes[1], target_nodes[2], tuning_c)
  }
  beta <- Theta[, target_nodes[1]] / Theta[target_nodes[1], target_nodes[1]];
  
  true_var <- find_var(Sigma, N, beta, key_ind_list$Theta_a[target_nodes[2],])
  E_term <- key_ind_list$Theta_a[target_nodes[2], ] %*% (data$Sigma_hat - Sigma) %*% beta
  z_score <- E_term / sqrt(true_var)
  test_out$z_score <- z_score
  test_out$true_var <- true_var
  test_out$test_knownvar <- test_out$test * sqrt(test_out$var_est/true_var)
  test_out$Entryest_Sigma <- data$Sigma_hat
  return(test_out)
}
test_pairwise_experiment_tilde2 <- function(p, target_nodes, n1, n2, seed_data, tuning_c = NULL, ...){
  ##given dimension p, graph type or a particular graph (Theta, Sigma, Adj_mat; or graph_type, seed_model), 
  # generate pairwise measurements data (setting seed_data) with sample sizes n1, n2, and other N[i,j];
  # n1 is scalar, n2 can be scalar or matrix of size (d_a+1)*(d_b^{(a)}+1).
  # then using the data to test whether there is an edge between target_nodes[1] and target_nodes[2].
  # difference to test_pairwise_experiment_tilde: this function saves Entryest_Sigma
  input <- list(...)
  #generate model
  if(sum(is.na(input$Theta))||sum(is.na(input$Sigma))||sum(is.na(input$Adj_mat))){
    if(is.na(input$graph_type)||is.na(input$model_seed)){
      print("Graph type and model_seed or graph parameters are missing!")
    }else{
      edge_val_range <- c(0.6,0.8)
      eigmin <- 0.25
      if(input$graph_type == "block"){
        model_param <- precision_generation(p,edge_val_range,eigmin,input$model_seed,input$graph_type,input$blocksize) 
      }else{
        model_param <- precision_generation(p,edge_val_range,eigmin,input$model_seed,input$graph_type) 
      }
      Theta <- model_param[[1]]
      Adj_mat <- model_param[[2]]
      diag(Adj_mat) <- 0
      Sigma <- solve(Theta)
    }
  }else{
    Theta <- input$Theta
    Sigma <- input$Sigma
    Adj_mat <- input$Adj_mat 
  }
  
  #setting sample sizes
  N <- matrix(rep(50, p^2), p, p)
  key_ind_list = find_key_ind(Adj_mat, Theta, p, target_nodes)
  N[key_ind_list$pair_ind1,] = n1; N[, key_ind_list$pair_ind1] = n1
  N[key_ind_list$pair_ind2$pair_ind2_1,key_ind_list$pair_ind2$pair_ind2_2] = n2; N[key_ind_list$pair_ind2$pair_ind2_2,key_ind_list$pair_ind2$pair_ind2_1] = t(n2)
  N_cv_train <- floor(N / 5 * 4); N_cv_test <- N - N_cv_train;
  #generate data
  data <- data_generation_pairwise(Sigma, N, seed_data)
  
  #test
  #update diagonal elements of N
  diag(N) <- apply(N, 1, sum); diag(N_cv_train) <- apply(N_cv_train, 1, sum); diag(N_cv_test) <- apply(N_cv_test, 1, sum);
  
  #N4 <- array(rep(0,p^4),rep(p,4))
  #for(j in 1:(p-1)){
  #  for(k in ((j+1):p)){
  #   N4[c(j,k),c(j,k),c(j,k),c(j,k)] <- N[j,k]
  #  }
  #}
  #for(j in 1:p){
  #  N4[j,j,j,j] <- N[j,j]
  #}
  if(is.null(tuning_c)){
    test_out <- edge_testing_tildevar(data$X, dim(data$X)[1], data$Ind, data$Sigma_hat, 
                                      p, N, target_nodes[1], target_nodes[2])
  }else{
    test_out <- edge_testing_tildevar(data$X, dim(data$X)[1], data$Ind, data$Sigma_hat, 
                                      p, N, target_nodes[1], target_nodes[2], tuning_c)
  }
  
  beta <- Theta[, target_nodes[1]] / Theta[target_nodes[1], target_nodes[1]];
  
  true_var <- find_var(Sigma, N, beta, key_ind_list$Theta_a[target_nodes[2],])
  E_term <- key_ind_list$Theta_a[target_nodes[2], ] %*% (data$Sigma_hat - Sigma) %*% beta
  z_score <- E_term / sqrt(true_var)
  test_out$z_score <- z_score
  test_out$true_var <- true_var
  test_out$test_knownvar <- test_out$test * sqrt(test_out$var_est/true_var)
  test_out$Entryest_Sigma <- data$Sigma_hat
  return(test_out)
}
test_pairwise_experiment_tilde <- function(p, target_nodes, n1, n2, seed_data, tuning_c = NULL, ...){
  ##given dimension p, graph type or a particular graph (Theta, Sigma, Adj_mat; or graph_type, seed_model), 
  # generate pairwise measurements data (setting seed_data) with sample sizes n1, n2, and other N[i,j];
  # n1 is scalar, n2 can be scalar or matrix of size (d_a+1)*(d_b^{(a)}+1).
  # then using the data to test whether there is an edge between target_nodes[1] and target_nodes[2].
  input <- list(...)
  #generate model
  if(sum(is.na(input$Theta))||sum(is.na(input$Sigma))||sum(is.na(input$Adj_mat))){
    if(is.na(input$graph_type)||is.na(input$model_seed)){
      print("Graph type and model_seed or graph parameters are missing!")
    }else{
      edge_val_range <- c(0.6,0.8)
      eigmin <- 0.25
      if(input$graph_type == "block"){
        model_param <- precision_generation(p,edge_val_range,eigmin,input$model_seed,input$graph_type,input$blocksize) 
      }else{
        model_param <- precision_generation(p,edge_val_range,eigmin,input$model_seed,input$graph_type) 
      }
      Theta <- model_param[[1]]
      Adj_mat <- model_param[[2]]
      diag(Adj_mat) <- 0
      Sigma <- solve(Theta)
    }
  }else{
    Theta <- input$Theta
    Sigma <- input$Sigma
    Adj_mat <- input$Adj_mat 
  }
  
  #setting sample sizes
  N <- matrix(rep(50, p^2), p, p)
  key_ind_list = find_key_ind(Adj_mat, Theta, p, target_nodes)
  N[key_ind_list$pair_ind1,] = n1; N[, key_ind_list$pair_ind1] = n1
  N[key_ind_list$pair_ind2$pair_ind2_1,key_ind_list$pair_ind2$pair_ind2_2] = n2; N[key_ind_list$pair_ind2$pair_ind2_2,key_ind_list$pair_ind2$pair_ind2_1] = t(n2)
  N_cv_train <- floor(N / 5 * 4); N_cv_test <- N - N_cv_train;
  #generate data
  data <- data_generation_pairwise(Sigma, N, seed_data)
  
  #test
  #update diagonal elements of N
  diag(N) <- apply(N, 1, sum); diag(N_cv_train) <- apply(N_cv_train, 1, sum); diag(N_cv_test) <- apply(N_cv_test, 1, sum);
  
  #N4 <- array(rep(0,p^4),rep(p,4))
  #for(j in 1:(p-1)){
  #  for(k in ((j+1):p)){
  #   N4[c(j,k),c(j,k),c(j,k),c(j,k)] <- N[j,k]
  #  }
  #}
  #for(j in 1:p){
  #  N4[j,j,j,j] <- N[j,j]
  #}
  if(is.null(tuning_c)){
    test_out <- edge_testing_tildevar(data$X, dim(data$X)[1], data$Ind, data$Sigma_hat, 
                             p, N, target_nodes[1], target_nodes[2])
  }else{
    test_out <- edge_testing_tildevar(data$X, dim(data$X)[1], data$Ind, data$Sigma_hat, 
                             p, N, target_nodes[1], target_nodes[2], tuning_c)
  }
  
  beta <- Theta[, target_nodes[1]] / Theta[target_nodes[1], target_nodes[1]];
  
  true_var <- find_var(Sigma, N, beta, key_ind_list$Theta_a[target_nodes[2],])
  E_term <- key_ind_list$Theta_a[target_nodes[2], ] %*% (data$Sigma_hat - Sigma) %*% beta
  z_score <- E_term / sqrt(true_var)
  test_out$z_score <- z_score
  test_out$true_var <- true_var
  test_out$test_knownvar <- test_out$test * sqrt(test_out$var_est/true_var)
  return(test_out)
}
#memory saving functions
data_generation_pairwise_savememory <- function(Sigma, N, seed){
  ## only return sample covarance
  p <- dim(Sigma)[1]
  set.seed(seed)
  diag_n <- rep(0, p); diag_Sigma <- rep(0, p);
  Sigma_hat <- matrix(rep(0, p^2), p, p)
  for(i in 1:(p-1)){
    for(j in (i+1):p){
      sqrt_Sigma_tmp <- msqrt(Sigma[c(i,j), c(i,j)])[[1]];
      X_tmp <- matrix(rnorm(N[i, j] * 2, 0, 1), N[i, j], 2) %*% sqrt_Sigma_tmp; 
      Sigma_hat[i,j] <- sum(X_tmp[, 1] * X_tmp[, 2])/N[i, j];
      diag_Sigma[i] <- diag_Sigma[i] + sum(X_tmp[, 1]^2);
      diag_n[i] <- diag_n[i] + N[i, j];
      diag_Sigma[j] <- diag_Sigma[j] + sum(X_tmp[,2]^2);
      diag_n[j] <- diag_n[j] + N[i, j];
    }
  }
  Sigma_hat <- Sigma_hat + t(Sigma_hat);
  diag(Sigma_hat) <- diag_Sigma / diag_n; 
  return(Sigma_hat)
}
test_pairwise_experiment_notuning <- function(p, target_nodes, n1, n2, seed_data, tuning_c, ...){
  ##given dimension p, graph type or a particular graph (Theta, Sigma, Adj_mat; or graph_type, seed_model), 
  # generate pairwise measurements data (setting seed_data) with sample sizes n1, n2, and other N[i,j];
  # n1 is scalar, n2 can be scalar or matrix of size (d_a+1)*(d_b^{(a)}+1).
  # then using the data to test whether there is an edge between target_nodes[1] and target_nodes[2].
  input <- list(...)
  #generate model
  if(is.null(input$Theta)||is.null(input$Sigma)||is.null(input$Adj_mat)){
    if(is.na(input$graph_type)||is.na(input$model_seed)){
      print("Graph type and model_seed or graph parameters are missing!")
    }else{
      edge_val_range <- c(0.6,0.8)
      eigmin <- 0.25
      if(input$graph_type == "block"){
        model_param <- precision_generation(p,edge_val_range,eigmin,input$model_seed,input$graph_type,input$blocksize) 
      }else{
        model_param <- precision_generation(p,edge_val_range,eigmin,input$model_seed,input$graph_type) 
      }
      Theta <- model_param[[1]]
      Adj_mat <- model_param[[2]]
      diag(Adj_mat) <- 0
      Sigma <- solve(Theta)
    }
  }else{
    Theta <- input$Theta
    Sigma <- input$Sigma
    Adj_mat <- input$Adj_mat 
  }
  
  #setting sample sizes
  N <- matrix(rep(50, p^2), p, p)
  key_ind_list = find_key_ind(Adj_mat, Theta, p, target_nodes)
  N[key_ind_list$pair_ind1,] = n1; N[, key_ind_list$pair_ind1] = n1
  N[key_ind_list$pair_ind2$pair_ind2_1,key_ind_list$pair_ind2$pair_ind2_2] = n2; N[key_ind_list$pair_ind2$pair_ind2_2,key_ind_list$pair_ind2$pair_ind2_1] = t(n2)
  ind <- upper.tri(N, diag = FALSE)
  n_total <- sum(N[ind]);
  #generate data
  Entryest_Sigma <- data_generation_pairwise_savememory(Sigma, N, seed_data)

  #test
  #update diagonal elements of N
  diag(N) <- apply(N, 1, sum); 
  
  #N4 <- array(rep(0,p^4),rep(p,4))
  #for(j in 1:(p-1)){
  #  for(k in ((j+1):p)){
  #   N4[c(j,k),c(j,k),c(j,k),c(j,k)] <- N[j,k]
  #  }
  #}
  #for(j in 1:p){
  #  N4[j,j,j,j] <- N[j,j]
  #}
  test_out <- edge_testing_notuning(Entryest_Sigma, n_total,
                             p, N, target_nodes[1], target_nodes[2], tuning_c)
  beta <- Theta[, target_nodes[1]] / Theta[target_nodes[1], target_nodes[1]];
  
  true_var <- find_var(Sigma, N, beta, key_ind_list$Theta_a[target_nodes[2],])
  E_term <- key_ind_list$Theta_a[target_nodes[2], ] %*% (Entryest_Sigma - Sigma) %*% beta
  z_score <- E_term / sqrt(true_var)
  test_out$z_score <- z_score
  test_out$true_var <- true_var
  test_out$test_knownvar <- test_out$test * sqrt(test_out$var_est/true_var)
  return(test_out)
}

find_E_term_experiment  <- function(p, target_nodes, n1, n2, seed_data, ...){
  ##debug for the conservative issue in the test: does the normal approximation for the E term work?
  # generate pairwise measurements data (setting seed_data) with sample sizes n1, n2, and other N[i,j];
  # n1 is scalar, n2 can be scalar or matrix of size (d_a+1)*(d_b^{(a)}+1).
  # then using the data to find Sigma_hat and E term test whether there is an edge between target_nodes[1] and target_nodes[2].
  input <- list(...)
  #generate model
  if(is.na(input$Theta)||is.na(input$Sigma)||is.na(input$Adj_mat)){
    if(is.na(input$graph_type)||is.na(input$model_seed)){
      print("Graph type and model_seed or graph parameters are missing!")
    }else{
      edge_val_range <- c(0.6,0.8)
      eigmin <- 0.25
      if(input$graph_type == "block"){
        model_param <- precision_generation(p,edge_val_range,eigmin,input$model_seed,input$graph_type,input$blocksize) 
      }else{
        model_param <- precision_generation(p,edge_val_range,eigmin,input$model_seed,input$graph_type) 
      }
      Theta <- model_param[[1]]
      Adj_mat <- model_param[[2]]
      diag(Adj_mat) <- 0
      Sigma <- solve(Theta)
    }
  }else{
    Theta <- input$Theta
    Sigma <- input$Sigma
    Adj_mat <- input$Adj_mat 
  }
  
  N <- matrix(rep(5, p^2), p, p)
  key_ind_list = find_key_ind(Adj_mat, Theta, p, target_nodes)
  N[key_ind_list$pair_ind1,] = n1; N[, key_ind_list$pair_ind1] = n1
  N[key_ind_list$pair_ind2$pair_ind2_1,key_ind_list$pair_ind2$pair_ind2_2] = n2; N[key_ind_list$pair_ind2$pair_ind2_2,key_ind_list$pair_ind2$pair_ind2_1] = t(n2)
  N_cv_train <- floor(N / 5 * 4); N_cv_test <- N - N_cv_train;
  #generate data
  data <- data_generation_pairwise(Sigma, N, N_cv_test, seed_data)
  
  #test
  #update diagonal elements of N
  diag(N) <- apply(N, 1, sum); diag(N_cv_train) <- apply(N_cv_train, 1, sum); diag(N_cv_test) <- apply(N_cv_test, 1, sum);
  
  N4 <- array(rep(0,p^4),rep(p,4))
  for(j in 1:(p-1)){
    for(k in ((j+1):p)){
      N4[c(j,k),c(j,k),c(j,k),c(j,k)] <- N[j,k]
    }
  }
  for(j in 1:p){
    N4[j,j,j,j] <- N[j,j]
  }
  beta <- Theta[, target_nodes[1]] / Theta[target_nodes[1], target_nodes[1]];
  E_term <- key_ind_list$Theta_a[target_nodes[2], ] %*% (data$Sigma_hat - Sigma) %*% beta
  
  true_var <- find_var(Sigma, N, N4, beta, key_ind_list$Theta_a[target_nodes[2],])
  z_score <- E_term / sqrt(true_var)
  return(list(E_term = E_term, true_var = true_var, z_score = z_score, p_val = 2*(1-pnorm(abs(z_score)))))
}
find_F1 <- function(est_graph, true_graph){
  confusion <- c(FP=sum(est_graph!=0&true_graph==0),FN=sum(est_graph==0&true_graph!=0),
                     TP=sum(est_graph!=0&true_graph!=0),TN=sum(est_graph==0&true_graph==0))
  F1 <- 2*confusion[3]/(2*confusion[3]+confusion[2]+confusion[1])
  confusion["F1"] = F1
  return(confusion)
}


find_true_var <- function(graph_type, p, n1, n2 = n1 * alpha, model_seed = 2021, edge_val_range = c(0.6,0.8), eigmin = 0.25, target_nodes = c(2,4)){
  model_param <- precision_generation(p,edge_val_range,eigmin,model_seed,graph_type)
  Theta <- model_param[[1]]
  Adj_mat <- model_param[[2]]
  diag(Adj_mat) <- 0
  Sigma <- solve(Theta)
  n2 <- n1 * alpha
  N <- matrix(rep(5, p^2), p, p)
  key_ind_list = find_key_ind(Adj_mat, Theta, p, target_nodes)
  N[key_ind_list$pair_ind1,] = n1; N[, key_ind_list$pair_ind1] = n1
  N[key_ind_list$pair_ind2$pair_ind2_1,key_ind_list$pair_ind2$pair_ind2_2] = n2; N[key_ind_list$pair_ind2$pair_ind2_2,key_ind_list$pair_ind2$pair_ind2_1] = t(n2)
  diag(N) <- apply(N, 1, sum);
       
  N4 <- array(rep(0,p^4),rep(p,4))
  for(j in 1:(p-1)){
   for(k in ((j+1):p)){
     N4[c(j,k),c(j,k),c(j,k),c(j,k)] <- N[j,k]
   }
  }
  for(j in 1:p){
    N4[j,j,j,j] <- N[j,j]
  }
  beta <- Theta[, target_nodes[1]] / Theta[target_nodes[1], target_nodes[1]];
     
  true_var <- find_var(Sigma, N, N4, beta, key_ind_list$Theta_a[target_nodes[2],])
}

