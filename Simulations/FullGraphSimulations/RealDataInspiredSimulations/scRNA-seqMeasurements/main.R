## real data - inspired simulations with real observational patterns: estimation and inference for the whole graph under different missing patterns; compare with baseline estimation and inference methods
## get observational pattern, graph type, seed from the command line

source("DataGeneration/DataGenerationFunctions.R")
source("DataGeneration/ExperimentFunctions.R")
source("Methods/Edgewise_GI-JOE.R")
source("Methods/baseline_graph_selection.R")
source("Methods/FullGraph_GI-JOE.R")
library(huge)
csv_path <- "Results/realdata_inspired_simulations/scRNA-seqMeasurements/F1_results.csv";
Obs_path <- "Results/realdata_inspired_simulations/scRNA-seqMeasurements/scRNAseq_Obs.RData"

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
  model_param <- precision_generation(p,edge_val_range,eigmin,model_seed,graph_type, prob = 3 / (p - 1))
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
nblasso_step1_results <- GIJOE_step1(X, Obs, n_total, p, thrs_c, signif_pval, onetuning = TRUE, symmetry = FALSE)
print(Sys.time() - start_t)
if(nblasso_step1_results$success){
  print("Finished our GI-JOE testing method step 1")
}else{
  print("Failed our GI-JOE testing method")
}
## step2 of GI-JOE full graph testing: variance estimation
p_kept <- length(nblasso_step1_results$kept_nodes)
row_ind <- NULL; for(i in 1 : (p_kept - 1)){row_ind <- c(row_ind, rep(i, p_kept - i))}
col_ind <- NULL; for(i in 2 : p_kept){col_ind <- c(col_ind, i : p_kept)}
nblasso_var_est <- find_variances_GIJOE(nblasso_step1_results, nblasso_step1_results$PSD_Sigma, nblasso_step1_results$theta_hat_varest, Obs, row_ind, col_ind)
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
graph_test_glasso <- test_glasso(graph_est[["glasso"]]$Sigma_hat, matrix(graph_est[["glasso"]]$Theta_hat, length(kept_nodes), length(kept_nodes)), 
                                 min(N[kept_nodes, kept_nodes]), p, kept_nodes, 0.05)
debiased_glasso_graph <- p_val_to_Holm_FDR(graph_test_glasso$p_val, p_kept = length(kept_nodes), p = p)
record_graph(debiased_glasso_graph$signif_graph_holm, method = "DBGlasso_Holm", graph_type, 
             measurement_scenario, n_total, seed, Adj_mat, csv_path)
record_graph(debiased_glasso_graph$signif_graph_FDR, method = "DBGlasso_FDR", graph_type, 
             measurement_scenario, n_total, seed, Adj_mat, csv_path)
