## simulation 1.2 with a different variance estimate (Theta_{b,b}^{(a)} estimation)
## PSD_Sigma instead of Entryest_Sigma is used for variance estimate

mylibpath <- "/Users/lilybean/Downloads/software/R_packages"
.libPaths(mylibpath)

wd <- "/Applications/Files/research/graph-quilting/graph-learning-varying-sample-size/simulations/inference_simulations"
setwd(wd)
RData_path <- "sim1_2/sim1_2_notuning/RData"
SigmaData_path <- "sim1_2/sim1_2_notuning/SigmaData_corrected"
csv_path <- "sim1_2/sim1_2_notuning/simulation1_2_results.csv"
source("fitting_functions.R")
source("sim_functions.R")
sim1_2_results <- read.table(file = "sim1_2/sim1_2_notuning/simulation1_2_results.csv", sep = " ", header = TRUE)
nrow(sim1_2_results)
ncol(sim1_2_results)
colnames(sim1_2_results)
sim1_2_results <- sim1_2_results[,c(1:13,17:18)]
sim1_2_results$test_varnewpsd <- rep(0, nrow(sim1_2_results));sim1_2_results$test_allnewpsd <- rep(0, nrow(sim1_2_results));sim1_2_results$var_est_newpsd <- rep(0, nrow(sim1_2_results))

for(graph_type in c("chain", "star", "ER")){
  for(alpha in c(1, 1.2, 1.5)){
    for(n2 in c(125, 250, 500)){
      for(signal_strength_ind in 1:8){
        signal_strength_vec <- c(0, 0.2, 0.4, 0.6, 0.8, 1.2, 2, 3) * sqrt(125 / n2)
        signal_strength <- signal_strength_vec[signal_strength_ind]
        p <- 200; seed_vec <- 1 : 200
        target_nodes <- c(2, 3); 
        model_seed <- 2021;
        edge_val_range <- c(0.6,0.8)
        eigmin <- 0.25
        
        if(graph_type == "ER"){
          #Have checked that with model_seed = 2021, p = 50 or 100, the generated graph does not have an edge (2,4)
          model_param <- precision_generation(p,edge_val_range,eigmin,model_seed,graph_type, 3/(p - 1), target_nodes = target_nodes, signal = signal_strength)
        }else if(graph_type == "star"){
          model_param <- precision_generation(p,edge_val_range,eigmin,model_seed,graph_type, 3, target_nodes = target_nodes, signal = signal_strength)
        }else{
          model_param <- precision_generation(p,edge_val_range,eigmin,model_seed,graph_type, target_nodes = target_nodes, signal = signal_strength)
        }
        Theta <- model_param[[1]];
        Adj_mat <- model_param[[2]]
        diag(Adj_mat) <- 0
        Sigma <- solve(Theta)
        
        n1 <- floor(n2 / alpha)
        var_est_new <- rep(0, length(seed_vec))
        test_varnew <- rep(0, length(seed_vec))
        test_allnew <- rep(0, length(seed_vec))
        load(paste(RData_path, sprintf("/%s_alpha_%d_n2_%d_signal_%d.RData", graph_type, 10*alpha, n2, signal_strength_ind), sep = ""))
        for(kk in 1 : length(seed_vec)){
          seed <- seed_vec[kk]
          ind <- which(sim1_2_results$graph == graph_type&sim1_2_results$n2==n2
                       &sim1_2_results$alpha==alpha&sim1_2_results$signal_ind==signal_strength_ind&sim1_2_results$seed == seed)
          Sigma_file <- paste(SigmaData_path, sprintf("/Sigma_%s_alpha_%d_n2_%d_signal_%d_seed_%d.RData", graph_type, 10*alpha, n2, signal_strength_ind, seed), sep = "")
          if(file.exists(Sigma_file)&sim1_2_results$var_est_newpsd[ind[1]]==0){
            N <- matrix(rep(50, p^2), p, p)
            key_ind_list = find_key_ind(Adj_mat, Theta, p, target_nodes)
            N[key_ind_list$pair_ind1,] = n1; N[, key_ind_list$pair_ind1] = n1
            N[key_ind_list$pair_ind2$pair_ind2_1,key_ind_list$pair_ind2$pair_ind2_2] = n2; N[key_ind_list$pair_ind2$pair_ind2_2,key_ind_list$pair_ind2$pair_ind2_1] = t(n2)
            #update diagonal elements of N
            diag(N) <- 0; diag(N) <- apply(N, 1, sum);
            load(Sigma_file)
            beta_hat2 <- -test_out[[kk]]$beta_hat;
            beta_hat2[target_nodes[1]] <- 1;
            
            theta_hat0 <- test_out[[kk]]$theta_hat
            theta_hat0 <- theta_hat0 / theta_hat0[target_nodes[2]]
            tau_new <- as.numeric(1 / (t(theta_hat0) %*% PSD_Sigma %*% theta_hat0))
            theta_hat_new <- theta_hat0 * tau_new
            beta_tilde_new = as.numeric(test_out[[kk]]$beta_hat[target_nodes[2]] - t(theta_hat_new) %*% (Entryest_Sigma %*% test_out[[kk]]$beta_hat - Entryest_Sigma[, target_nodes[1]]))
            
            var_est_new[kk] <- find_var(PSD_Sigma, N, beta_hat2, theta_hat_new)
            test_varnew[kk] <- test_out[[kk]]$beta_tilde/sqrt(var_est_new[kk])
            test_allnew[kk] <- beta_tilde_new/sqrt(var_est_new[kk])
            
            sim1_2_results$test_varnewpsd[ind] <- test_varnew[kk]
            sim1_2_results$test_allnewpsd[ind] <- test_allnew[kk]
            sim1_2_results$var_est_newpsd[ind] <- var_est_new[kk]
            rm(PSD_Sigma)          
            }
        }
        #system("say Yay!")
      }
    }
  }
}
sum(sim1_2_results$var_est_newpsd!=0)
nrow(sim1_2_results)
mean(pnorm(abs(sim1_2_results$test_varnewpsd))>0.975)
write.table(sim1_2_results, file = csv_path, sep = " ", row.names=FALSE)
