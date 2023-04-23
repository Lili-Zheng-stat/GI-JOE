## simulation 2: estimation and inference for the whole graph under different missing patterns; compare with baseline estimation and inference methods
## get graph type, measurement scenarios, total sample sizes, seed from the command line

source("DataGeneration/DataGenerationFunctions.R")
source("DataGeneration/ExperimentFunctions.R")
source("Methods/Edgewise_GI-JOE.R")
source("Methods/baseline_graph_selection.R")
source("Methods/FullGraph_GI-JOE.R")
library(huge)
csv_path <- "Results/sim2/F1_results.csv";

args <- commandArgs(trailingOnly=TRUE);
graph_type <- args[[1]]
measurement_scenario <- as.numeric(args[[2]])
n_total <- as.numeric(args[[3]])
seed <- as.numeric(args[[4]])

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

## GI-JOE full graph testing
thrs_c <- 5; signif_pval <- 0.05
start_t <- Sys.time()
## step1: first perform neighborhood regression and debiasing for all node pairs
nblasso_step1_results <- GIJOE_step1(X, Obs, n_total, p, thrs_c, signif_pval, onetuning = TRUE, symmetry = FALSE)
print(Sys.time() - start_t)
if(nblasso_step1_results$success){
  print("Finished our nblasso testing method step 1")
}else{
  print("Failed our nblasso testing method")
}
## step2 of GI-JOE full graph testing: variance estimation
p_kept <- length(nblasso_step1_results$kept_nodes)
row_ind <- NULL; for(i in 1 : (p_kept - 1)){row_ind <- c(row_ind, rep(i, p_kept - i))}
col_ind <- NULL; for(i in 2 : p_kept){col_ind <- c(col_ind, i : p_kept)}
nblasso_var_est <- find_variances_GIJOE(nblasso_step1_results, nblasso_step1_results$test_out$PSD_Sigma, nblasso_step1_results$test_out$theta_hat_varest, Obs, row_ind, col_ind)
var_est_results <- data.frame(rows = nblasso_step1_results$kept_nodes[row_ind], cols =  nblasso_step1_results$kept_nodes[col_ind], var_est = nblasso_var_est)

## step3 of GI-JOE full graph testing: compute edge-wise test statistic and apply Holm's correction or FDR control procedure on p-values
record_graphs_GIJOE_step3(nblasso_step1_results, var_est_results, p, signif_pval = 0.05, screenSigma = FALSE, graph_type, 
                            measurement_scenario, n_total, seed, Adj_mat, csv_path)

## run baseline graph selection methods
baseline_methods <- c("Nlasso(AND)", "Nlasso(OR)", "Glasso", "CLIME"); 
kept_nodes <- which(apply(N, 1, sum) > 0)
for(baseline in baseline_methods){
  start_t <- Sys.time()
  print(paste("running with baseline method", baseline))
  graph_est[[baseline]] <- find_graph_stb_tuning(X, Obs, n_total, p, kept_nodes, baseline, zerofill = FALSE)
  record_graph(graph_est[[baseline]]$graph_est, method = baseline, graph_type, 
                   measurement_scenario, n_total, seed, Adj_mat, csv_path)
  print(paste("Finished running baseline method", baseline))
  print(Sys.time() - start_t)
}

## run debiased graphical lasso with minimum pairwise sample size
graph_test_glasso <- test_glasso(graph_est[["Glasso"]]$Sigma_hat, matrix(graph_est[["Glasso"]]$Theta_hat, length(kept_nodes), length(kept_nodes)), 
                                 min(N[kept_nodes, kept_nodes]), p, kept_nodes, 0.05)
debiased_glasso_graph <- p_val_to_Holm_FDR(graph_test_glasso$p_val, p_kept = length(kept_nodes), p = p)
record_graph(debiased_glasso_graph$signif_graph_holm, method = "DBGlasso_Holm", graph_type, 
                 measurement_scenario, n_total, seed, Adj_mat, csv_path)
record_graph(debiased_glasso_graph$signif_graph_FDR, method = "DBGlasso_FDR", graph_type, 
                 measurement_scenario, n_total, seed, Adj_mat, csv_path)

