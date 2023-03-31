## step 3: summarize debiased and variance estimation results
wd <- "/Applications/Files/research/graph-quilting/graph-learning-varying-sample-size/simulations/inference_simulations"
setwd(wd)
mylibpath <- "/Users/lilybean/Downloads/software/R_packages"
.libPaths(mylibpath)
RData_path <- "ABA_experiments/RData/"
varest_path <- "ABA_experiments/varest_results_updated/"
setwd(wd)
source("fitting_functions_missing.R")
source("fitting_functions_glasso.R")
source("sim2_functions.R")
source("sim_functions.R")
source("baseline_fitting.R")
source("sim2/step3/step3_functions.R")
source("sim2/step2/variance_estimate_methods.R")

##determine exprID and seed and preprocessing seed
preprocessing_seed <- 2059
exprID <- 4
seed <- 2059


#load(preprocessing_path)
load(paste(RData_path, sprintf("raw_thrs_expr%d_seed%d.RData", exprID, seed), sep = ""))
p <- ncol(Obs)

screenSigma <- FALSE
varest_file_name <- paste(varest_path, sprintf("raw_thrs_var_est_expr%d_seed%d_batch1.csv", exprID, seed), sep = "")
var_est_results <- read.table(varest_file_name, sep = " ", header = TRUE)
p_kept <- length(nblasso_step1_results$kept_nodes); 
for(batch in 2 : ceiling(p_kept * (p_kept - 1) / 2 / 200)){
  varest_file_name <- paste(varest_path, sprintf("raw_thrs_var_est_expr%d_seed%d_batch%d.csv", exprID, seed, batch), sep = "")
  if(file.exists(varest_file_name)){
    var_est_results_add <- read.table(varest_file_name, sep = " ", header = TRUE)
    var_est_results <- rbind(var_est_results, var_est_results_add)
  }else{
    print(sprintf("Variance estimate for batch %d not found!", batch))
    break;
  }
}

#summarize and visualize results with different threshold level
#setting testing threshold
testing_thrs <- 0.12
if(nblasso_step1_results$success){
  nblasso_inference_results <- find_graphs_nblasso_step3(nblasso_step1_results, var_est_results, p, signif_pval = 0.05, testing_thrs = testing_thrs)
}
N <- t(Obs)%*%Obs; kept_nodes <- which(apply(N, 1, sum) != 0)
graph_test_glasso <- test_glasso(graph_est[["glasso"]][[2]]$Sigma_hat, matrix(graph_est[["glasso"]][[2]]$Theta_hat, length(kept_nodes), length(kept_nodes)), 
                                 min(N[kept_nodes, kept_nodes]), p, kept_nodes, 0.05, testing_thrs = testing_thrs * mean(diag(matrix(graph_est[["glasso"]][[2]]$Theta_hat,length(kept_nodes), length(kept_nodes)))))
debiased_glasso_graph <- p_val_to_Holm_FDR(graph_test_glasso$p_val, p_kept = length(kept_nodes), p = p)



load(sprintf("ABA_experiments/thrs_raw_preprocessed_seed%d.RData", preprocessing_seed))
#redo the inference with new threshold
graph_test_full <- find_graphs_nblasso_step3(nblasso_step1_results_full, var_est_full, p, signif_pval = 0.05, testing_thrs = testing_thrs)


#visualize graphs using adjacency matrix
neuron_location <- read.csv("/Applications/Files/research/graph-quilting/RealData/AllenBrainAtlasData/neuron_location.csv", header = TRUE, row.names=1)
if(exprID == 2 || exprID == 3){
  node_order <- order(neuron_location[,2])
}else if(exprID == 4){
  node_order <- order(neuron_location[,1], decreasing = TRUE)
}
#node_set <- node_order[151:227]
#node_set <- node_order[76:150]
#node_set <- node_order[1:75]

#F1 results for the whole graph
node_set <- 1:p
a <- find_F1(nblasso_inference_results$signif_graph_FDR[node_set,node_set], graph_test_full$signif_graph_FDR[node_set,node_set])
a <- rbind(a, find_F1(nblasso_inference_results$signif_graph_holm[node_set,node_set], graph_test_full$signif_graph_FDR[node_set,node_set]))
a <- rbind(a, find_F1(debiased_glasso_graph$signif_graph_FDR[node_set,node_set], graph_test_full$signif_graph_FDR[node_set,node_set]))
a <- rbind(a, find_F1(debiased_glasso_graph$signif_graph_holm[node_set,node_set], graph_test_full$signif_graph_FDR[node_set,node_set]))
results <- data.frame(Method =c("GI-JOE(FDR)","GI-JOE(Holm)", "DB-Glasso(FDR)","DB-Glasso(Holm)"), FP = a[,'FP'], FN = a[, 'FN'], TP = a[, 'TP'], TN = a[, 'TN'], F1 = a[, 'F1'])
results$TPR <- results$TP / (results$FN + results$TP)
results$TNR <- results$TN / (results$TN + results$FP)
results$TDR <- results$TP / (results$FP + results$TP)
xtable(results[,c("Method", "TPR","TNR","TDR", "F1")], digits=3)

#F1 results for each node set
node_label <- c(rep("low", floor(p / 3)), rep("median", floor(p / 3)), rep("high", p - floor(p / 3) * 2))
for(node_set_label in c("high", "median", "low")){
  node_set <- node_order[which(node_label == node_set_label)]
  a <- find_F1(nblasso_inference_results$signif_graph_FDR[node_set,node_set], graph_test_full$signif_graph_FDR[node_set,node_set])
  a <- rbind(a, find_F1(nblasso_inference_results$signif_graph_holm[node_set,node_set], graph_test_full$signif_graph_FDR[node_set,node_set]))
  a <- rbind(a, find_F1(debiased_glasso_graph$signif_graph_FDR[node_set,node_set], graph_test_full$signif_graph_FDR[node_set,node_set]))
  a <- rbind(a, find_F1(debiased_glasso_graph$signif_graph_holm[node_set,node_set], graph_test_full$signif_graph_FDR[node_set,node_set]))
  results <- data.frame(Method =c("GI-JOE(FDR)","GI-JOE(Holm)", "DB-Glasso(FDR)","DB-Glasso(Holm)"), FP = a[,'FP'], FN = a[, 'FN'], TP = a[, 'TP'], TN = a[, 'TN'], F1 = a[, 'F1'])
  results$TPR <- results$TP / (results$FN + results$TP)
  results$TNR <- results$TN / (results$TN + results$FP)
  results$TDR <- results$TP / (results$FP + results$TP)
  if(node_set_label == "high"){
    full_results <- results[,c("Method", "TPR","TNR","TDR", "F1")]
  }else{
    full_results <- cbind(full_results, results[,c("TPR","TNR","TDR", "F1")])
  }
}
full_results
xtable(full_results, digits=3)



plot_graph(nblasso_inference_results$signif_graph_FDR[node_set,node_set], graph_test_full$signif_graph_FDR[node_set,node_set])
plot_graph(debiased_glasso_graph$signif_graph_FDR[node_set,node_set], graph_test_full$signif_graph_FDR[node_set,node_set])

plot_graph(graph_test_full$signif_graph_FDR[node_order,node_order], graph_test_full$signif_graph_FDR[node_order,node_order])
plot_graph(nblasso_inference_results$signif_graph_FDR[node_order,node_order], nblasso_inference_results$signif_graph_FDR[node_order,node_order])
plot_graph(debiased_glasso_graph$signif_graph_FDR[node_order,node_order], debiased_glasso_graph$signif_graph_FDR[node_order,node_order])

ggsave(plot_graph(nblasso_inference_results$signif_graph_FDR[node_order,node_order], graph_test_full$signif_graph_FDR[node_order,node_order]),
       filename = sprintf("ABA_experiments/fig_updated/adjacency_GIJOE_raw_thrs%.2f_expr%d_seed%d.png", testing_thrs, exprID, seed))
ggsave(plot_graph(debiased_glasso_graph$signif_graph_FDR[node_order,node_order], graph_test_full$signif_graph_FDR[node_order,node_order]),
       filename = sprintf("ABA_experiments/fig_updated/adjacency_debiasedGlasso_raw_thrs%.2f_expr%d_seed%d.png", testing_thrs, exprID, seed))


# neuronal tuning info
angular_tuning <- read.table("/Applications/Files/research/graph-quilting/RealData/AllenBrainAtlasData/neuron_orient_tuning.csv", sep = ",", header = TRUE)
freq_tuning <- read.table("/Applications/Files/research/graph-quilting/RealData/AllenBrainAtlasData/neuron_freq_tuning.csv", sep = ",", header = TRUE)
angular_tuning[,2] <- angular_tuning[,2] * 45


comp_angular <- matrix(rep(as.vector(angular_tuning[,2]), p), p, p, byrow = TRUE) == matrix(rep(as.vector(angular_tuning[,2]), p), p, p, byrow = FALSE)
comp_freq <- matrix(rep(as.vector(freq_tuning[,2]), p), p, p, byrow = TRUE) == matrix(rep(as.vector(freq_tuning[,2]), p), p, p, byrow = FALSE)

# neuron locations for plots
library(tidyverse)
neuron_location$measurement_level <- factor(c(rep("low", floor(p / 3)), rep("moderate", floor(p / 3)), rep("high", p - 2 * floor(p / 3))), 
                                            levels = c("low", "moderate", "high"))
neuron_location$angular_tuning <- as.factor(angular_tuning[,2])
neuron_location$freq_tuning <- as.factor(freq_tuning[,2])
#given a node set, find F1 scores, plot graphs
ABA_summarize_graphs_F1 <- function(fitted_graph, true_graph, node_set, neuronal_tuning){
  F1_results <- find_F1(fitted_graph[node_set, node_set], true_graph[node_set, node_set])
  TPR <- F1_results['TP'] / (F1_results['TP'] + F1_results['FN'])
  TNR <- F1_results['TN'] / (F1_results['TN'] + F1_results['FP'])
  precision <- F1_results['TP'] / (F1_results['TP'] + F1_results['FP'])
  edge_ratio_same_tuning <- sum(fitted_graph[node_set, node_set]!=0 & neuronal_tuning[node_set, node_set]!=0)/sum(fitted_graph[node_set, node_set]!=0)
  same_tuning_ratio_edge <- sum(fitted_graph[node_set, node_set]!=0 & neuronal_tuning[node_set, node_set]!=0)/sum(neuronal_tuning[node_set, node_set]!=0)
  return(c(F1_results, TPR = TPR, TNR = TNR, precision = precision, precision_tuning = edge_ratio_same_tuning, recall_tuning = same_tuning_ratio_edge))
}



ABA_visualize_graphs <- function(fitted_graph, node_set, neuronal_tuning, node_order = NULL){
  # Plot points in space
  if(neuronal_tuning == "angular"){
    base_plot <- ggplot() +
      geom_point(aes(x = V1, y = V2, col = angular_tuning), data = neuron_location[node_set,], size = 1) +
      theme_bw()
  }else if(neuronal_tuning == "frequency"){
    base_plot <- ggplot() +
      geom_point(aes(x = V1, y = V2, col = freq_tuning), data = neuron_location[node_set,], size = 1) +
      theme_bw()
  }else if(!is.null(node_order)){
    p <- 227; neuron_location$neuron_type <- rep(0, p)
    neuron_location$neuron_type[node_order[1 : floor(p / 3)]] <- "low"
    neuron_location$neuron_type[node_order[(1 + floor(p / 3)) : (2 * floor(p / 3))]] <- "median"
    neuron_location$neuron_type[node_order[(1 + 2 * floor(p / 3)) : p]] <- "high"
    base_plot <- ggplot() +
      geom_point(aes(x = V1, y = V2, col = neuron_type), data = neuron_location[node_set,], size = 1) +
      theme_bw()
  }
  
  #print(base_plot)
  fitted_edges <- as.data.frame(which(fitted_graph[node_set, node_set] != 0, arr.ind = TRUE))
  fitted_edges$x1 <- neuron_location[node_set[fitted_edges[, 1]], "V1"]
  fitted_edges$y1 <- neuron_location[node_set[fitted_edges[, 1]], "V2"]
  fitted_edges$x2 <- neuron_location[node_set[fitted_edges[, 2]], "V1"]
  fitted_edges$y2 <- neuron_location[node_set[fitted_edges[, 2]], "V2"]
  if(neuronal_tuning == "angular"){
    p <- base_plot + 
      geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = fitted_edges,
                   alpha = 0.7, size = 0.2) +
      theme_bw() + 
      labs(x = "neuron location (x)", y = "neuron location (y)", color = "Angular tuning")
  }else if(neuronal_tuning == "frequency"){
    p <- base_plot + 
      geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = fitted_edges,
                   alpha = 0.7, size = 0.2) +
      theme_bw() + 
      labs(x = "neuron location (x)", y = "neuron location (y)", color = "Frequency tuning")
  }else{
    p <- base_plot + 
      geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = fitted_edges,
                   alpha = 0.7, size = 0.2) +
      theme_bw() + 
      labs(x = "neuron location (x)", y = "neuron location (y)", color = "Sample size")
  }
  return(p)
}

ggsave(ABA_visualize_graphs(graph_test_full$signif_graph_FDR, 1:p, "sample_size", node_order),
       filename = sprintf("ABA_experiments/fig_updated/true_graph_raw_thrs%.2f_expr%d_seed%d.png", testing_thrs, exprID, seed), 
       width = 6.5, height = 5)
ggsave(ABA_visualize_graphs(nblasso_inference_results$signif_graph_FDR, 1:p, "sample_size", node_order),
       filename = sprintf("ABA_experiments/fig_updated/GIJOE_raw_thrs%.2f_expr%d_seed%d.png", testing_thrs, exprID, seed),
       width = 6.5, height = 5)
ggsave(ABA_visualize_graphs(debiased_glasso_graph$signif_graph_FDR, 1:p, "sample_size", node_order),
       filename = sprintf("ABA_experiments/fig_updated/debiasedGlasso_raw_thrs%.2f_expr%d_seed%d.png", testing_thrs, exprID, seed),
       width = 6.5, height = 5)

ABA_summarize_graphs_F1(nblasso_inference_results$signif_graph_FDR, graph_test_full$signif_graph_FDR, node_order[151:227], comp_angular)
ABA_summarize_graphs_F1(nblasso_inference_results$signif_graph_holm, graph_test_full$signif_graph_FDR, node_order[151:227], comp_angular)
ABA_summarize_graphs_F1(debiased_glasso_graph$signif_graph_FDR, graph_test_full$signif_graph_FDR, node_order[151:227], comp_angular)
ABA_summarize_graphs_F1(debiased_glasso_graph$signif_graph_holm, graph_test_full$signif_graph_FDR, node_order[151:227], comp_angular)

ABA_summarize_graphs_F1(nblasso_inference_results$signif_graph_FDR, graph_test_full$signif_graph_FDR, 1:p, comp_angular)
ABA_summarize_graphs_F1(nblasso_inference_results$signif_graph_holm, graph_test_full$signif_graph_FDR, 1:p, comp_angular)
ABA_summarize_graphs_F1(debiased_glasso_graph$signif_graph_FDR, graph_test_full$signif_graph_FDR, 1:p, comp_angular)
ABA_summarize_graphs_F1(debiased_glasso_graph$signif_graph_holm, graph_test_full$signif_graph_FDR, 1:p, comp_angular)

