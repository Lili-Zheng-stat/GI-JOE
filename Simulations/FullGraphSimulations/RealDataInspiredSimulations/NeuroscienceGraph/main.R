## real data - inspired simulations: estimation and inference for the whole graph under different missing patterns; compare with baseline estimation and inference methods
## get graph type, measurement scenarios, total sample sizes, seed from the command line

source("DataGeneration/DataGenerationFunctions.R")
source("DataGeneration/ExperimentFunctions.R")
source("Methods/Edgewise_GI-JOE.R")
source("Methods/baseline_graph_selection.R")
source("Methods/FullGraph_GI-JOE.R")
library(huge)
csv_path <- "Results/realdata_inspired_simulations/NeuroscienceGraph/F1_results.csv";
preprocessing_path <- "Results/NeuroscienceGraph/preprocessed.RData"

args <- commandArgs(trailingOnly=TRUE);
measurement_scenario <- as.numeric(args[[1]])
n_total <- as.numeric(args[[2]])
seed <- as.numeric(args[[3]])

load(preprocessing_path)
p <- 227; 
graph_type <- "ABA"

#generate model
Theta <- (Theta_full + t(Theta_full)) / 2
Theta[abs(Theta)<0.05 & abs(t(Theta))<0.05] <- 0; Adj_mat <- Theta!=0; diag(Adj_mat) <- 0
Sigma <- solve(Theta)


#measurement scenario
Obs_seed <- 1000
Obs <- generate_Obs(p, n_total, Adj_mat, measurement_scenario, Obs_seed)
N <- t(Obs)%*%Obs

#generate data
eigen_out <- eigen(Sigma)

set.seed(seed)
print(sprintf("running experiment with measurement scenario %d, n_total = %d, seed =%d", measurement_scenario, n_total, seed))
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
