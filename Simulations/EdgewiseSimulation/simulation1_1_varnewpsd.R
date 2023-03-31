## simulation 1.1 with a different variance estimate (Theta_{b,b}^{(a)} estimation)
## PSD_Sigma instead of Entryest_Sigma is used for variance estimate

mylibpath <- "/Users/lilybean/Downloads/software/R_packages"
.libPaths(mylibpath)

wd <- "/Applications/Files/research/graph-quilting/graph-learning-varying-sample-size/simulations/inference_simulations"
setwd(wd)
RData_path <- "sim1_1/sim1_1_notuning/RData_tilde"
SigmaData_path <- "sim1_1/sim1_1_notuning/SigmaData_corrected"
csv_path <- "sim1_1/sim1_1_notuning/simulation1_1_results_tilde.csv"
source("fitting_functions.R")
source("sim_functions.R")
sim1_1_results_tilde <- read.table(file = "sim1_1/sim1_1_notuning/simulation1_1_results_tilde.csv", sep = " ", header = TRUE)
ncol(sim1_1_results_tilde)
sim1_1_results_tilde$test_varnewpsd <- rep(0, nrow(sim1_1_results_tilde));sim1_1_results_tilde$test_allnewpsd <- rep(0, nrow(sim1_1_results_tilde));
sim1_1_results_tilde$var_estnewpsd <- rep(0, nrow(sim1_1_results_tilde))

for(graph_type in c("chain", "star", "ER")){
  for(alpha in c(1, 1.2, 1.5)){
    for(p in c(50, 100, 200)){
      for(n1 in c(250, 500, 1000)){
        target_nodes <- c(2,4); 
        model_seed <- 2021;
        edge_val_range <- c(0.6,0.8)
        eigmin <- 0.25
        seed_vec <- 1 : 1000
        
        if(graph_type == "ER"){
          #Have checked that with model_seed = 2021, p = 50 or 100, the generated graph does not have an edge (2,4)
          model_param <- precision_generation(p,edge_val_range,eigmin,model_seed,graph_type, prob = 3/(p - 1))
        }else if(graph_type == "star"){
          model_param <- precision_generation(p,edge_val_range,eigmin,model_seed,graph_type, nstar = 3)
        }else{
          model_param <- precision_generation(p,edge_val_range,eigmin,model_seed,graph_type)
        }
        Theta <- model_param[[1]]
        Adj_mat <- model_param[[2]]
        diag(Adj_mat) <- 0
        Sigma <- solve(Theta)
        
        n2 <- floor(n1 * alpha); 
        var_est_new <- rep(0, length(seed_vec))
        test_varnew <- rep(0, length(seed_vec))
        test_allnew <- rep(0, length(seed_vec))
        load(paste(RData_path, sprintf("/tilde_%s_alpha_%d_p_%d_n1_%d.RData", graph_type, 10*alpha, p, n1), sep = ""))
        for(kk in 1 : length(seed_vec)){
          seed <- seed_vec[kk]
          ind <- which(sim1_1_results_tilde$graph == graph_type&sim1_1_results_tilde$p==p&sim1_1_results_tilde$n1==n1
                       &sim1_1_results_tilde$alpha==alpha&sim1_1_results_tilde$seed == seed)
          if(length(ind)==0){
            test_tmp <- test_out[[kk]]
            result <- data.frame(graph = graph_type, p = p, n1 = n1, alpha = alpha, seed = seed, test_stat = test_tmp$test, var_est = test_tmp$var_est, var_est_tilde = test_tmp$var_est_tilde,
                                 test_tilde = test_tmp$test_tilde, test_knownvar = test_tmp$test_knownvar, E_term_score = test_tmp$z_score, time = 0,
                                 test_varnew = 0, test_allnew = 0, var_est_new = 0)
            sim1_1_results_tilde <- rbind(sim1_1_results_tilde, result)
            ind <- which(sim1_1_results_tilde$graph == graph_type&sim1_1_results_tilde$p==p&sim1_1_results_tilde$n1==n1
                         &sim1_1_results_tilde$alpha==alpha&sim1_1_results_tilde$seed == seed)
          }
          Sigma_file <- paste(SigmaData_path, sprintf("/Sigma_%s_alpha_%d_p_%d_n1_%d_seed_%d.RData", graph_type, 10*alpha, p, n1, seed), sep = "")
          if(file.exists(Sigma_file)&sim1_1_results_tilde$var_estnewpsd[ind[1]]==0){
            N <- matrix(rep(50, p^2), p, p)
            key_ind_list = find_key_ind(Adj_mat, Theta, p, target_nodes)
            N[key_ind_list$pair_ind1,] = n1; N[, key_ind_list$pair_ind1] = n1
            N[key_ind_list$pair_ind2$pair_ind2_1,key_ind_list$pair_ind2$pair_ind2_2] = n2; N[key_ind_list$pair_ind2$pair_ind2_2,key_ind_list$pair_ind2$pair_ind2_1] = t(n2)
            #update diagonal elements of N
            diag(N) <- 0; diag(N) <- apply(N, 1, sum);
            load(Sigma_file)
            beta_hat2 <- -test_out[[kk]]$beta_hat;
            beta_hat2[target_nodes[1]] <- 1;
            
            theta_hat0 <- test_out[[kk]]$theta_hat[target_nodes[2],]
            theta_hat0 <- theta_hat0 / theta_hat0[target_nodes[2]]
            tau_new <- as.numeric(1 / (t(theta_hat0) %*% PSD_Sigma %*% theta_hat0))
            theta_hat_new <- theta_hat0 * tau_new
            beta_tilde_new = as.numeric(test_out[[kk]]$beta_hat[target_nodes[2]] - t(theta_hat_new) %*% (Entryest_Sigma %*% test_out[[kk]]$beta_hat - Entryest_Sigma[, target_nodes[1]]))
            
            var_est_new[kk] <- find_var(PSD_Sigma, N, beta_hat2, theta_hat_new)
            test_varnew[kk] <- test_out[[kk]]$beta_tilde[target_nodes[2]]/sqrt(var_est_new[kk])
            test_allnew[kk] <- beta_tilde_new/sqrt(var_est_new[kk])
            
            sim1_1_results_tilde$test_varnewpsd[ind] <- test_varnew[kk]
            sim1_1_results_tilde$test_allnewpsd[ind] <- test_allnew[kk]
            sim1_1_results_tilde$var_estnewpsd[ind] <- var_est_new[kk]
            rm(PSD_Sigma)          
          }
        }
        print(sprintf("Finished graph %s, alpha = %.2f, p = %d, n1 = %d", graph_type, alpha, p, n1))
      }
    }
  }
}
sum(sim1_1_results_tilde$var_estnewpsd!=0)
nrow(sim1_1_results_tilde)
mean(pnorm(abs(sim1_1_results_tilde$test_varnewpsd))>0.975)
write.table(sim1_1_results_tilde, file = csv_path, sep = " ", row.names=FALSE)
