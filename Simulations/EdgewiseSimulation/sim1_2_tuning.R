##find tuning constant c for edge-wise simulations, validating power
## graph type, alpha, n2 and signal strength are determined by command line arguments, 


csv_path <- "Results/sim1_2/simulation1_2_tuning.csv"
source("DataGeneration/DataGenerationFunctions.R")
source("DataGeneration/ExperimentFunctions.R")
source("Methods/Edgewise_GI-JOE.R")


args <- commandArgs(trailingOnly=TRUE);
graph_type <- args[[1]]
alpha <- as.numeric(args[[2]])
n2 <- as.numeric(args[[3]])
signal_strength_ind <- as.numeric(args[[4]])

signal_strength_vec <- c(0, 0.2, 0.4, 0.6, 0.8, 1.2, 2, 3) * sqrt(125 / n2)
signal_strength <- signal_strength_vec[signal_strength_ind]
p <- 200; seed <- 1
target_nodes <- c(2, 3); 
#alpha is the ratio between n2 and n1
#generate the same graph for all replicates
model_seed <- 2021;
edge_val_range <- c(0.6,0.8)
eigmin <- 0.25

if(graph_type == "ER"){
  model_param <- precision_generation(p,edge_val_range,eigmin,model_seed,graph_type, 3/(p - 1), target_nodes = target_nodes, signal = signal_strength)
}else if(graph_type == "star"){
  model_param <- precision_generation(p,edge_val_range,eigmin,model_seed,graph_type, 3, target_nodes = target_nodes, signal = signal_strength)
}else{
  model_param <- precision_generation(p,edge_val_range,eigmin,model_seed,graph_type, target_nodes = target_nodes, signal = signal_strength)
}
Theta <- model_param[[1]];
Adj_mat <- model_param[[2]]
diag(Adj_mat) <- 0
Sigma <- solve(Theta)

n1 <- floor(n2 / alpha)
print(sprintf("running experiment with %s graph, signal_ind = %d, n2 = %d, alpha = %.2f, seed =%d", graph_type, signal_strength_ind, n2, alpha, seed))
start_t <- Sys.time()
test_tmp <- test_pairwise_experiment(p, target_nodes, n1, n2, seed, Theta = Theta, Sigma = Sigma, Adj_mat = Adj_mat)
comp_time <- Sys.time() - start_t
print(comp_time)
tuning_results <- data.frame(graph = graph_type, signal_ind = signal_strength_ind, n2 = n2, alpha = alpha, tuning_c1 = test_tmp$tuning_c1, tuning_c2 = test_tmp$tuning_c2)
if(file.exists(csv_path)){
  write.table(tuning_results, file = csv_path, sep = " ", row.names=FALSE, append = TRUE, col.names = FALSE)
}else{
  write.table(tuning_results, file = csv_path, sep = " ", row.names=FALSE)
}



