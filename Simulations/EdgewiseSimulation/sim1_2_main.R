## simulation 1.1 (H_A)
## graph type, alpha, n2, signal_strength are determined by command line arguments, 
## loop over seed

csv_path <- "Results/sim1_2/simulation1_2_results.csv"
precision_csv_path <- "Results/precision_CI_results.csv"
tuning_path <- "Results/sim1_2/simulation1_2_tuning.csv"
source("DataGeneration/DataGenerationFunctions.R")
source("DataGeneration/ExperimentFunctions.R")
source("Methods/Edgewise_GI-JOE.R")

args <- commandArgs(trailingOnly=TRUE);
graph_type <- args[[1]]
alpha <- as.numeric(args[[2]])
n2 <- as.numeric(args[[3]])
signal_strength_ind <- as.numeric(args[[4]])

signal_strength_vec <- c(0, 0.2, 0.4, 0.6, 0.8, 1.2, 2, 3) *sqrt(125 / n2)
signal_strength <- signal_strength_vec[signal_strength_ind]
p <- 200; seed_vec <- 1 : 200
target_nodes <- c(2, 3); n1 <- floor(n2 / alpha)
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
truesignal <- Theta[target_nodes[1], target_nodes[2]] / Theta[target_nodes[1], target_nodes[1]]


tuning_results <- read.table(file = tuning_path, sep = " ", header = TRUE)
tuning_c <- tuning_results[tuning_results$graph == graph_type & tuning_results$signal_ind == signal_strength_ind & 
                             tuning_results$n2 == n2 & tuning_results$alpha == alpha, c("tuning_c1", "tuning_c2")]
if(nrow(tuning_c) > 1){
  tuning_c <- as.numeric(tuning_c[1, ])
}else{
  tuning_c <- as.numeric(tuning_c)
}

#run experiments for seed in seed_vec
for(kk in 1:length(seed_vec)){
  seed <- seed_vec[kk]
  print(sprintf("running experiment with %s graph, signal_ind = %d, n1 = %d, alpha = %.2f, seed =%d", graph_type, signal_strength_ind, n1, alpha, seed))
  test_tmp <- test_pairwise_experiment(p, target_nodes, n1, n2, seed, tuning_c = tuning_c, Theta = Theta, Sigma = Sigma, Adj_mat = Adj_mat, precisionCI = TRUE)
  result <- data.frame(graph = graph_type, signal_ind = signal_strength_ind, truesignal = truesignal, n2 = n2, alpha = alpha, seed = seed, test_stat = test_tmp$test, var_est = test_tmp$var_est,
                       test_knownvar = test_tmp$test_knownvar, var_true = test_tmp$true_var)
  if(file.exists(csv_path)){
    write.table(result, file = csv_path, sep = " ", row.names=FALSE, append = TRUE, col.names = FALSE)
  }else{
    write.table(result, file = csv_path, sep = " ", row.names=FALSE)
  }
  coverage = (test_tmp$CI_Theta[1] <= Theta[target_nodes[1], target_nodes[2]] & test_tmp$CI_Theta[2] >= Theta[target_nodes[1], target_nodes[2]])
  result_precision <- data.frame(graph = graph_type, p = p, n1 = n1, alpha = alpha, n2 = n2, signal_ind = signal_strength_ind, seed = seed, coverage = coverage)
  if(file.exists(precision_csv_path)){
    write.table(result_precision, file = precision_csv_path, sep = " ", row.names=FALSE, append = TRUE, col.names = FALSE)
  }else{
    write.table(result_precision, file = precision_csv_path, sep = " ", row.names=FALSE)
  }
}
