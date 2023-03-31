## simulation 2: estimation and inference for the whole graph under different missing patterns; compare with baseline estimation and inference methods
## get graph type, measurement scenarios, total sample sizes, seed from the command line
## step 1: for our inference methods, only perform the debiasing step, leave the variance estimation and FMER, FDR for step 2

.libPaths("/home/liliz/R_packages")
wd <- "/home/liliz/Graph_DiffN/simulations"
RData_path <- "sim2/step1/RData/"
csv_path <- "sim2/F1_results.csv";
setwd(wd)
source("fitting_functions_missing.R")
source("fitting_functions_glasso.R")
source("sim2_functions.R")
source("sim_functions.R")
source("baseline_fitting.R")
source("sim2/step1/Inference_methods_step1.R")

library(huge)

args <- commandArgs(trailingOnly=TRUE);
graph_type <- args[[1]]
measurement_scenario <- as.numeric(args[[2]])
n_total <- as.numeric(args[[3]])
seed <- as.numeric(args[[4]])

#graph_type <- "chain"; measurement_scenario <- 4; n_total <- 20000; seed <- 1

p <- 200; 
model_seed <- 2021;
edge_val_range <- c(0.6,0.8)
eigmin <- 0.25

#generate model
if(graph_type == "star"){
  model_param <- precision_generation(p,edge_val_range,eigmin,model_seed,graph_type, nstar = 10)
}else if(graph_type == "ER"){
  model_param <- precision_generation(p,edge_val_range,eigmin,model_seed,graph_type, prob = 3 / (p - 1))
}else{
  model_param <- precision_generation(p,edge_val_range,eigmin,model_seed,graph_type)
}

Theta <- model_param[[1]]
Adj_mat <- model_param[[2]]
diag(Adj_mat) <- 0
Sigma <- solve(Theta)

#measurement scenario
Obs_seed <- 1000
Obs <- generate_Obs(p, n_total, Adj_mat, measurement_scenario, Obs_seed)
N <- t(Obs)%*%Obs

#generate data
eigen_out <- eigen(Sigma)

set.seed(seed)
print(sprintf("running experiment with %s graph, measurement scenario %d, n_total = %d, seed =%d", graph_type, measurement_scenario, n_total, seed))
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
save(Obs, nblasso_step1_results, file = paste(RData_path, sprintf("NDS_%s%dn%dexpr%d.RData", graph_type, measurement_scenario, n_total, seed), sep = ""))


##other naive estimation baselines: graphical lasso; neighborhood lasso; CLIME
#save confusion and F1; better figures;
#baseline_methods <- c("nblasso_and", "nblasso_or", "glasso", "CLIME"); graph_est <- list()
baseline_methods <- c("nblasso_and", "nblasso_or", "glasso"); graph_est <- list(list())
for(baseline in baseline_methods){
  start_t <- Sys.time()
  print(paste("running with baseline method", baseline))
  #for(zerofill in c(TRUE, FALSE)){
  zerofill <- FALSE
    graph_est[[baseline]][[2-as.numeric(zerofill)]] <- find_graph_stb_tuning(X, Obs, n_total, p, 1 : p, baseline, zerofill = zerofill)
    record_graph_new(graph_est[[baseline]][[2-as.numeric(zerofill)]]$graph_est, method = paste(sprintf("baseline_zerofill%s", zerofill), baseline, sep = "_"), graph_type, 
                     measurement_scenario, n_total, seed, Adj_mat, csv_path)
    print(paste("Finished running baseline method", baseline, sprintf("zerofill = %s", zerofill)))
    print(Sys.time() - start_t)
  #}
}
save(Obs, nblasso_step1_results, graph_est, file = paste(RData_path, sprintf("NDS_%s%dn%dexpr%d.RData", graph_type, measurement_scenario, n_total, seed), sep = ""))


##baseline inference methods
#graph_test_glasso <- test_glasso(graph_est[["glasso"]]$Sigma_hat, graph_est[["glasso"]]$Theta_hat[test_results$kept_nodes,test_results$kept_nodes], 
#                                 min(N[test_results$kept_nodes,test_results$kept_nodes]), length(test_results$kept_nodes), test_results$kept_nodes, 0.05)
start_t <- Sys.time()
graph_test_glasso <- test_glasso(graph_est[["glasso"]][[2]]$Sigma_hat, graph_est[["glasso"]][[2]]$Theta_hat, 
                                 min(N[N > 5 * sqrt(log(p))]), p, 1 : p, 0.05)
record_graph_new(graph_test_glasso$signif_graph_holm, method = "baseline_debiased_glasso", graph_type, 
                 measurement_scenario, n_total, seed, Adj_mat, csv_path)
print("Finished baseline method: graphical lasso inference")
print(Sys.time() - start_t)

save(Obs, nblasso_step1_results, graph_est, graph_test_glasso, file = paste(RData_path, sprintf("NDS_%s%dn%dexpr%d.RData", graph_type, measurement_scenario, n_total, seed), sep = ""))
