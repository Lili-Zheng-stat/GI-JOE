##step 1 of ABA experiment, raw traces, spontaneous period, testing against a threshold
source("DataGeneration/DataGenerationFunctions.R")
source("DataGeneration/ExperimentFunctions.R")
source("Methods/Edgewise_GI-JOE.R")
source("Methods/baseline_graph_selection.R")
source("Methods/FullGraph_GI-JOE.R")
library(huge)


preprocessing_path <- "NeuroscienceDataAnalysis/preprocessed.RData"
load(preprocessing_path)
n_total <- nrow(X_full); p <- ncol(X_full)
neuron_location <- read.csv("NeuroscienceDataAnalysis/neuron_location.csv", header = TRUE, row.names=1)
node_order <- order(neuron_location[,1], decreasing = TRUE)

##missing patterns
Obs_seed <- 1000
#neurons are independently measured with probabilities sqrt(0.1), sqrt(0.5), or sqrt(0.9)
Obs <- generate_Obs_AllenBrainAtlas(p, n_total, Adj_mat, Obs_seed, node_order)
N <- t(Obs)%*%Obs
X <- X_full * Obs

#find p-value matrix
start_t <- Sys.time()
seed <- 2059
set.seed(seed)
nblasso_step1_results <- GIJOE_step1(X, Obs, n_total, p, thrs_c = 5, signif_pval = 0.05, onetuning = TRUE, symmetry = FALSE)
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
nblasso_var_est <- find_variances_GIJOE(nblasso_step1_results, nblasso_step1_results$PSD_Sigma, nblasso_step1_results$theta_hat_varest, Obs, row_ind, col_ind)
var_est_results <- data.frame(rows = nblasso_step1_results$kept_nodes[row_ind], cols =  nblasso_step1_results$kept_nodes[col_ind], var_est = nblasso_var_est)

## step3 of GI-JOE full graph testing: compute edge-wise test statistic and apply Holm's correction or FDR control procedure on p-values
testing_thrs <- 0.12
if(nblasso_step1_results$success){
  GIJOE_results <- GIJOE_thrs_step3(nblasso_step1_results, var_est_results, p, signif_pval = 0.05, testing_thrs = testing_thrs)
}

## debiased graphical lasso
kept_nodes <- which(apply(N, 1, sum) != 0)
baseline <- "Glasso"
start_t <- Sys.time()
print(paste("running with baseline method", baseline))
graph_est <- find_graph_stb_tuning(X, Obs, n_total, p, kept_nodes, baseline, zerofill = zerofill)
print(paste("Finished running baseline method", baseline, sprintf("zerofill = %s", zerofill)))
print(Sys.time() - start_t)

graph_test_glasso <- test_glasso(graph_est$Sigma_hat, matrix(graph_est$Theta_hat, length(kept_nodes), length(kept_nodes)), 
                                 min(N[kept_nodes, kept_nodes]), p, kept_nodes, 0.05, 
                                 testing_thrs = testing_thrs * mean(diag(matrix(graph_est$Theta_hat,length(kept_nodes), length(kept_nodes)))))
debiased_glasso_graph <- p_val_to_Holm_FDR(graph_test_glasso$p_val, p_kept = length(kept_nodes), p = p)

save(graph_test_full, debiased_glasso_graph, GIJOE_results, file = "/Applications/Files/GitHub/GI-JOE/NeuroscienceDataAnalysis/AnalyzedResults.RData")

