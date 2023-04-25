## summarize real data experiment results
load("NeuroscienceDataAnalysis/AnalyzedResults.RData")
node_order <- order(neuron_location[,1], decreasing = TRUE)

#F1 results for the whole graph
node_set <- 1:p
a <- find_F1(GIJOE_results$signif_graph_FDR[node_set,node_set], graph_test_full$signif_graph_FDR[node_set,node_set])
a <- rbind(a, find_F1(GIJOE_results$signif_graph_holm[node_set,node_set], graph_test_full$signif_graph_FDR[node_set,node_set]))
a <- rbind(a, find_F1(debiased_glasso_graph$signif_graph_FDR[node_set,node_set], graph_test_full$signif_graph_FDR[node_set,node_set]))
a <- rbind(a, find_F1(debiased_glasso_graph$signif_graph_holm[node_set,node_set], graph_test_full$signif_graph_FDR[node_set,node_set]))
results <- data.frame(Method =c("GI-JOE(FDR)","GI-JOE(Holm)", "DB-Glasso(FDR)","DB-Glasso(Holm)"), FP = a[,'FP'], FN = a[, 'FN'], TP = a[, 'TP'], TN = a[, 'TN'], F1 = a[, 'F1'])
results$TPR <- results$TP / (results$FN + results$TP)
results$TNR <- results$TN / (results$TN + results$FP)
results$TDR <- results$TP / (results$FP + results$TP)
results

#F1 results for each node set
node_label <- c(rep("low", floor(p / 3)), rep("median", floor(p / 3)), rep("high", p - floor(p / 3) * 2))
for(node_set_label in c("high", "median", "low")){
  node_set <- node_order[which(node_label == node_set_label)]
  a <- find_F1(GIJOE_results$signif_graph_FDR[node_set,node_set], graph_test_full$signif_graph_FDR[node_set,node_set])
  a <- rbind(a, find_F1(GIJOE_results$signif_graph_holm[node_set,node_set], graph_test_full$signif_graph_FDR[node_set,node_set]))
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


# plot neuron graph 
library(tidyverse)
neuron_location$measurement_level <- factor(c(rep("low", floor(p / 3)), rep("moderate", floor(p / 3)), rep("high", p - 2 * floor(p / 3))), 
                                            levels = c("low", "moderate", "high"))
ABA_visualize_graphs <- function(fitted_graph, node_set, node_order){
  # Plot points in space
  p <- 227; neuron_location$neuron_type <- rep(0, p)
  neuron_location$neuron_type[node_order[1 : floor(p / 3)]] <- "low"
  neuron_location$neuron_type[node_order[(1 + floor(p / 3)) : (2 * floor(p / 3))]] <- "median"
  neuron_location$neuron_type[node_order[(1 + 2 * floor(p / 3)) : p]] <- "high"
  base_plot <- ggplot() +
    geom_point(aes(x = V1, y = V2, col = neuron_type), data = neuron_location[node_set,], size = 1) +
    theme_bw()
  
  fitted_edges <- as.data.frame(which(fitted_graph[node_set, node_set] != 0, arr.ind = TRUE))
  fitted_edges$x1 <- neuron_location[node_set[fitted_edges[, 1]], "V1"]
  fitted_edges$y1 <- neuron_location[node_set[fitted_edges[, 1]], "V2"]
  fitted_edges$x2 <- neuron_location[node_set[fitted_edges[, 2]], "V1"]
  fitted_edges$y2 <- neuron_location[node_set[fitted_edges[, 2]], "V2"]
  p <- base_plot + 
    geom_segment(aes(x = x1, y = y1, xend = x2, yend = y2), data = fitted_edges,
                 alpha = 0.7, size = 0.2) +
    theme_bw() + 
    labs(x = "neuron location (x)", y = "neuron location (y)", color = "Sample size")
  return(p)
}

ABA_visualize_graphs(graph_test_full$signif_graph_FDR, 1:p, node_order)
ABA_visualize_graphs(GIJOE_results$signif_graph_FDR, 1:p,  node_order)
ABA_visualize_graphs(debiased_glasso_graph$signif_graph_FDR, 1:p, node_order)

