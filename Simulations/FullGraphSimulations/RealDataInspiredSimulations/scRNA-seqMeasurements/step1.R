## real data - inspired simulations with real observational patterns: estimation and inference for the whole graph under different missing patterns; compare with baseline estimation and inference methods
## get observational pattern, graph type, seed from the command line
## step 1: for our inference methods, only perform the debiasing step, leave the variance estimation and FMER, FDR for step 2

.libPaths("/home/liliz/R_packages")
wd <- "/home/liliz/Graph_DiffN/simulations"
RData_path <- "scRNA_seq_simulation/step1/RData/"
csv_path <- "scRNA_seq_simulation/F1_results.csv";
Obs_path <- "scRNA_seq_simulation/scRNAseq_Obs.RData"
setwd(wd)
source("fitting_functions_missing.R")
source("fitting_functions_glasso.R")
source("sim2_functions.R")
source("sim_functions.R")
source("baseline_fitting.R")
source("sim2/step1/Inference_methods_step1.R")

library(igraph)
library(huge)

args <- commandArgs(trailingOnly=TRUE);
measurement_setting <- as.numeric(args[[1]])
graph_type <- args[[2]]
seed <- as.numeric(args[[3]])


load(Obs_path)
Obs <- Obs_list[[measurement_setting]]
n_total <- nrow(Obs); p <- ncol(Obs)
N <- t(Obs)%*%Obs

model_seed <- 2021;
edge_val_range <- c(0.6,0.8)
eigmin <- 0.25

#generate model
if(graph_type == "star"){
  model_param <- precision_generation(p,edge_val_range,eigmin,model_seed,graph_type, nstar = 10)
}else if(graph_type == "ER"){
  model_param <- precision_generation(p,edge_val_range,eigmin,model_seed,graph_type, prob = 5 / p)
}else if(graph_type == "small-world"){
  model_param <- precision_generation(p,edge_val_range,eigmin,model_seed,graph_type, degree = 2)
}else if(graph_type == "scale-free"){
  model_param <- precision_generation(p,edge_val_range,eigmin,model_seed,graph_type, new_edge = 1)
  Theta <- model_param[[1]][p:1,p:1]
  Adj_mat <- model_param[[2]][p:1,p:1]
}else{
  model_param <- precision_generation(p,edge_val_range,eigmin,model_seed,graph_type)
}
if(graph_type != "scale-free"){
  Theta <- model_param[[1]]
  Adj_mat <- model_param[[2]]
}
  
diag(Adj_mat) <- 0
Sigma <- solve(Theta)

#generate data
eigen_out <- eigen(Sigma)

set.seed(seed)
print(sprintf("running experiment with measurement setting %d, %s graph, seed =%d", measurement_setting, graph_type, seed))
Z <- matrix(rnorm(p * n_total, 0, 1), n_total, p)%*%eigen_out[[2]]%*%diag(sqrt(eigen_out[[1]]))%*%t(eigen_out[[2]])
X <- Z * Obs; 

#find p-value matrix
thrs_c <- 5; signif_pval <- 0.05
start_t <- Sys.time()
nblasso_step1_results <- nblasso_step1(X, Obs, n_total, p, thrs_c, signif_pval, onetuning = TRUE, symmetry = FALSE)
print(Sys.time() - start_t)
if(nblasso_step1_results$success){
  print("Finished our nblasso testing method step 1")
}else{
  print("Failed our nblasso testing method")
}
save(Obs, nblasso_step1_results, file = paste(RData_path, sprintf("NDS_%d%sexpr%d.RData", measurement_setting, graph_type, seed), sep = ""))


##other naive estimation baselines: graphical lasso; neighborhood lasso; CLIME
#save confusion and F1; better figures;
#baseline_methods <- c("nblasso_and", "nblasso_or", "glasso", "CLIME"); graph_est <- list()
kept_nodes <- which(apply(N, 1, sum) > 0)
baseline_methods <- c("nblasso_and", "nblasso_or", "glasso"); graph_est <- list(list())
for(baseline in baseline_methods){
  start_t <- Sys.time()
  print(paste("running with baseline method", baseline))
  #for(zerofill in c(TRUE, FALSE)){
  zerofill <- FALSE
    graph_est[[baseline]][[2-as.numeric(zerofill)]] <- find_graph_stb_tuning(X, Obs, n_total, p, kept_nodes, baseline, zerofill = zerofill)
    print(paste("Finished running baseline method", baseline, sprintf("zerofill = %s", zerofill)))
    print(Sys.time() - start_t)
  #}
}
save(Obs, nblasso_step1_results, graph_est, file = paste(RData_path, sprintf("NDS_%d%sexpr%d.RData", measurement_setting, graph_type, seed), sep = ""))


##baseline inference methods
#graph_test_glasso <- test_glasso(graph_est[["glasso"]]$Sigma_hat, graph_est[["glasso"]]$Theta_hat[test_results$kept_nodes,test_results$kept_nodes], 
#                                 min(N[test_results$kept_nodes,test_results$kept_nodes]), length(test_results$kept_nodes), test_results$kept_nodes, 0.05)
start_t <- Sys.time()
graph_test_glasso <- test_glasso(graph_est[["glasso"]][[2]]$Sigma_hat, matrix(graph_est[["glasso"]][[2]]$Theta_hat, length(kept_nodes), length(kept_nodes)), 
                                 min(N[kept_nodes, kept_nodes]), p, kept_nodes, 0.05)
#debiased_glasso_graph <- p_val_to_Holm_FDR(graph_test_glasso$p_val, p_kept = length(kept_nodes), p = p)

print("Finished baseline method: graphical lasso inference")
print(Sys.time() - start_t)

save(Obs, nblasso_step1_results, graph_est, graph_test_glasso, file = paste(RData_path, sprintf("NDS_%d%sexpr%d.RData", measurement_setting, graph_type, seed), sep = ""))
