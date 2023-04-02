## functions for simulations
#----------------------------------------------------#
#perform edge-wise testing simulation for one replicate
test_pairwise_experiment <- function(p, target_nodes, n1, n2, seed_data, tuning_c = NULL, ...){
  ## given dimension p, graph type or a particular graph (Theta, Sigma, Adj_mat; or graph_type, seed_model), 
  # generate pairwise measurements data (setting seed_data) with sample sizes n1, n2, and other N[i,j];
  # n1 is scalar, n2 can be scalar or matrix of size (d_a+1)*(d_b^{(a)}+1).
  # then using the data to test whether there is an edge between target_nodes[1] and target_nodes[2].
  ## return: a list of output, including test statistic, estimated variance, true variance, neighborhood regression results, tuning parameters, etc.
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
  #generate data
  data <- data_generation_pairwise(Sigma, N, seed_data)
  
  #test
  #update diagonal elements of N
  diag(N) <- 0; diag(N) <- apply(N, 1, sum);
  if(is.null(tuning_c)){
    test_out <- edge_testing(data$X, dim(data$X)[1], data$Ind, data$Sigma_hat, 
                                       p, N, target_nodes[1], target_nodes[2])
  }else{
    test_out <- edge_testing(data$X, dim(data$X)[1], data$Ind, data$Sigma_hat, 
                                       p, N, target_nodes[1], target_nodes[2], tuning_c)
  }
  beta <- Theta[, target_nodes[1]] / Theta[target_nodes[1], target_nodes[1]];
  
  true_var <- find_var(Sigma, N, beta, key_ind_list$Theta_a[target_nodes[2],])
  test_out$true_var <- true_var
  test_out$test_knownvar <- test_out$test * sqrt(test_out$var_est/true_var)
  test_out$Entryest_Sigma <- data$Sigma_hat
  return(test_out)
}

#find key neighborhood for assigning pairwise sample sizes in edge-wise simulations
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
find_F1 <- function(est_graph, true_graph){
  ## find the FP, FN, TP, TN, F1 scpre of an estiamted graph (est_graph is its adjacency matrix) compared with true_graph
  confusion <- c(FP=sum(est_graph!=0&true_graph==0),FN=sum(est_graph==0&true_graph!=0),
                     TP=sum(est_graph!=0&true_graph!=0),TN=sum(est_graph==0&true_graph==0))
  F1 <- 2*confusion[3]/(2*confusion[3]+confusion[2]+confusion[1])
  confusion["F1"] = F1
  return(confusion)
}

record_graphs_GIJOE_step3 <- function(step1_results, var_est_results, p, signif_pval, screenSigma, graph_type, 
                                        measurement_scenario, n_total, seed, Adj_mat, csv_path){
  selected_graphs <- GIJOE_step3(step1_results, var_est_results, p, signif_pval, screenSigma)
  #record F1 score for tested and estimated graphs
  record_graph(selected_graphs$signif_graph_Bonferroni, method = "GIJOE_Bonf", graph_type, 
                   measurement_scenario, n_total, seed, Adj_mat, csv_path)
  record_graph(selected_graphs$signif_graph_holm, method = "GIJOE_Holm", graph_type, 
                   measurement_scenario, n_total, seed, Adj_mat, csv_path)
  record_graph(selected_graphs$signif_graph_FDR, method = "GIJOE_FDR", graph_type, 
                   measurement_scenario, n_total, seed, Adj_mat, csv_path)
  record_graph(selected_graphs$graph_est_AND, method = "Nlasso_JOE(AND)", graph_type, 
                   measurement_scenario, n_total, seed, Adj_mat, csv_path)
  record_graph(selected_graphs$graph_est_OR, method = "Nlasso_JOE(OR)", graph_type, 
                   measurement_scenario, n_total, seed, Adj_mat, csv_path)
}


record_graph <- function(graph_est, method, graph_type, measurement_scenario, n_total, seed, true_graph, csv_path){
  #write its F1 score and confusion to a csv file
  p <- dim(graph_est)[1]
  graph_est_true <- matrix(rep("0", p^2), p, p)
  graph_est_true[graph_est!=0&true_graph==0] <- "FP"; graph_est_true[graph_est!=0&true_graph!=0] <- "TP";
  graph_est_true[graph_est==0&true_graph==0] <- "TN"; graph_est_true[graph_est==0&true_graph!=0] <- "FN";
  
  F1_out <- find_F1(graph_est, true_graph)
  F1_result <- data.frame(graph_type = graph_type, ms_sc = measurement_scenario, n = n_total, seed = seed, est_method = method, FP = F1_out["FP"], 
                          FN = F1_out["FN"], TP = F1_out["TP"], TN = F1_out["TN"], F1 = F1_out['F1'])
  if(file.exists(csv_path)){
    write.table(F1_result, file = csv_path, sep = " ", row.names=FALSE, append = TRUE, col.names = FALSE)
  }else{
    write.table(F1_result, file = csv_path, sep = " ", row.names=FALSE)
  }
}


