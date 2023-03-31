##step 1 of ABA experiment, raw traces, spontaneous period, testing against a threshold
wd <- "/Applications/Files/research/graph-quilting/graph-learning-varying-sample-size/simulations/inference_simulations"
setwd(wd)
mylibpath <- "/Users/lilybean/Downloads/software/R_packages"
.libPaths(mylibpath)
RData_path <- "ABA_experiments/RData/"
source("fitting_functions_missing.R")
source("fitting_functions_glasso.R")
source("sim2_functions.R")
source("sim_functions.R")
source("baseline_fitting.R")
source("sim2/step1/Inference_methods_step1.R")
source("ABA_experiments/ABA_functions")
library(huge)

preprocessing_seed <- 2059
preprocessing_path <- sprintf("ABA_experiments/thrs_raw_preprocessed_seed%d.RData", preprocessing_seed)
seed <- 2059

load(preprocessing_path)
n_total <- nrow(X_full); p <- ncol(X_full)
neuron_location <- read.csv("/Applications/Files/research/graph-quilting/RealData/AllenBrainAtlasData/neuron_location.csv", header = TRUE, row.names=1)
node_order <- order(neuron_location[,1], decreasing = TRUE)

##missing patterns
exprID <- 4
Obs_seed <- 1000
#probabilities 0.1, 0.5, 0.9
measurement_scenario <- 3
Obs <- generate_Obs_ABA(p, n_total, Adj_mat, measurement_scenario, Obs_seed, node_order)
N <- t(Obs)%*%Obs
hist(N)

g <- ggplot(data.frame(node1 = rep(1 : p, each = p), node2= rep(1 : p, p), 
                       sample_sizes = as.vector(N[node_order, node_order])), aes(x = node1, y = node2))+
  geom_raster(aes(fill = sample_sizes)) + scale_fill_gradient(high = "dodgerblue4", low = "ghostwhite")+
  theme_bw() + theme(axis.text.x=element_text(size=9, angle=0, vjust=0.3),
                     axis.text.y=element_text(size=9),
                     plot.title=element_text(size=11)) + 
  guides(fill = guide_legend(title="Sample sizes")) + scale_y_reverse()

X <- X_full * Obs


#find p-value matrix
thrs_c <- 5; signif_pval <- 0.05
start_t <- Sys.time()
set.seed(seed)
nblasso_step1_results <- nblasso_step1(X, Obs, n_total, p, thrs_c, signif_pval, onetuning = TRUE, symmetry = FALSE)
print(Sys.time() - start_t)
if(nblasso_step1_results$success){
  print("Finished our nblasso testing method step 1")
}else{
  print("Failed our nblasso testing method")
}
save(Obs, nblasso_step1_results, file = paste(RData_path, sprintf("raw_thrs_expr%d_seed%d.RData", exprID, seed), sep = ""))


##other naive estimation baselines: graphical lasso; neighborhood lasso; CLIME
#save confusion and F1; better figures;
#baseline_methods <- c("nblasso_and", "nblasso_or", "glasso", "CLIME"); graph_est <- list()
kept_nodes <- which(apply(N, 1, sum) != 0)
baseline_methods <- c("nblasso_and", "nblasso_or", "glasso"); graph_est <- list(list())
for(baseline in baseline_methods){
  start_t <- Sys.time()
  print(paste("running with baseline method", baseline))
  zerofill <- FALSE
  graph_est[[baseline]][[2-as.numeric(zerofill)]] <- find_graph_stb_tuning(X, Obs, n_total, p, kept_nodes, baseline, zerofill = zerofill)
  print(paste("Finished running baseline method", baseline, sprintf("zerofill = %s", zerofill)))
  print(Sys.time() - start_t)
}
save(Obs, nblasso_step1_results, graph_est, file = paste(RData_path, sprintf("raw_thrs_expr%d_seed%d.RData", exprID, seed), sep = ""))


##baseline inference methods
#graph_test_glasso <- test_glasso(graph_est[["glasso"]]$Sigma_hat, graph_est[["glasso"]]$Theta_hat[test_results$kept_nodes,test_results$kept_nodes], 
#                                 min(N[test_results$kept_nodes,test_results$kept_nodes]), length(test_results$kept_nodes), test_results$kept_nodes, 0.05)
start_t <- Sys.time()
graph_test_glasso <- test_glasso(graph_est[["glasso"]][[2]]$Sigma_hat, matrix(graph_est[["glasso"]][[2]]$Theta_hat, length(kept_nodes), length(kept_nodes)), 
                                 min(N[kept_nodes, kept_nodes]), p, kept_nodes, 0.05)
print("Finished baseline method: graphical lasso inference")
print(Sys.time() - start_t)

save(Obs, nblasso_step1_results, graph_est, graph_test_glasso, file = paste(RData_path, sprintf("raw_thrs_expr%d_seed%d.RData", exprID, seed), sep = ""))

