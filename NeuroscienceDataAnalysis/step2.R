## step 2: for our inference methods, finding the variance estimation for the batch edges 
wd <- "/Applications/Files/research/graph-quilting/graph-learning-varying-sample-size/simulations/inference_simulations"
setwd(wd)
mylibpath <- "/Users/lilybean/Downloads/software/R_packages"
.libPaths(mylibpath)
RData_path <- "ABA_experiments/RData/"
csv_path <- "ABA_experiments/varest_results_updated/"
source("fitting_functions_missing.R")
source("fitting_functions_glasso.R")
source("sim2_functions.R")
source("sim_functions.R")
source("baseline_fitting.R")
source("sim2/step1/Inference_methods_step1.R")
source("sim2/step2/variance_estimate_methods.R")
library(huge)

args <- commandArgs(trailingOnly=TRUE);
print(length(args))
exprID <- as.numeric(args[[1]])
seed <- as.numeric(args[[2]])
batch_edges <- as.numeric(args[[3]])
batch_size <- 200

print(sprintf("start running experiment %d, seed %d, with %dth batch of edges", exprID, seed, batch_edges))
load(paste(RData_path, sprintf("raw_thrs_expr%d_seed%d.RData", exprID, seed), sep = ""))
p_kept <- length(nblasso_step1_results$kept_nodes)

## project Sigma
Entryest_Sigma <- nblasso_step1_results$test_out$Entryest_Sigma
N <- nblasso_step1_results$test_out$N
p_kept <- length(nblasso_step1_results$kept_nodes)
proj_out <- proj_cov(Entryest_Sigma, N, 0.01, 1, matrix(rep(0, p_kept^2), p_kept, p_kept), matrix(rep(1, p_kept^2), p_kept, p_kept), 0.001, 500)
PSD_Sigma  <- proj_out$projected_S

## find new theta
theta_hat <- nblasso_step1_results$test_out$theta_hat
theta_hat_updated <- theta_hat
for(node_a in 1 : (p_kept-1)){
  ind <- (1 : p_kept)[-node_a]
  for(node_b in ind){
    if(node_b > node_a){
      theta_hat0 <- theta_hat[node_a, node_b, ] / theta_hat[node_a, node_b, node_b]
      theta_hat_updated[node_a, node_b, ] <- theta_hat0 * as.numeric(1/t(theta_hat0) %*% PSD_Sigma %*% theta_hat0)
    }
  }
}

row_ind_set <- NULL; for(i in 1 : (p_kept - 1)){row_ind_set <- c(row_ind_set, rep(i, p_kept - i))}
col_ind_set <- NULL; for(i in 2 : p_kept){col_ind_set <- c(col_ind_set, i : p_kept)}
if((batch_edges - 1) * batch_size + 1 <= length(row_ind_set)){
  if(batch_edges * batch_size < length(row_ind_set)){
    row_ind <- row_ind_set[((batch_edges - 1) * batch_size + 1) : (batch_edges * batch_size)]; 
    col_ind <- col_ind_set[((batch_edges - 1) * batch_size + 1) : (batch_edges * batch_size)];
  }else{
    row_ind <- row_ind_set[((batch_edges - 1) * batch_size + 1) : length(row_ind_set)]; 
    col_ind <- col_ind_set[((batch_edges - 1) * batch_size + 1) : length(row_ind_set)];
  }
  
  #graph_type <- "chain"; measurement_scenario <- 4; n_total <- 20000; seed <- 1
  
  if(nblasso_step1_results$success){
    nblasso_var_est <- find_variances_nblasso_updated(nblasso_step1_results, PSD_Sigma, theta_hat_updated, Obs, row_ind, col_ind)
    results <- data.frame(rows = nblasso_step1_results$kept_nodes[row_ind], cols =  nblasso_step1_results$kept_nodes[col_ind], var_est = nblasso_var_est, method = "nblasso")
    file_name <- paste(csv_path, sprintf("raw_thrs_var_est_expr%d_seed%d_batch%d.csv", exprID, seed, batch_edges), sep = "")
    write.table(results, file = file_name, sep = " ", row.names=FALSE)
  }
}




