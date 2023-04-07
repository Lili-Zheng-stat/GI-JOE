## simulation 1.1 (H0)
## graph type, alpha, and p, n1 are determined by command line arguments, 
## loop over seed
csv_path <- "Results/sim1_1/simulation1_1_results.csv"
precision_csv_path <- "Results/precision_CI_results.csv"
tuning_path <- "Results/sim1_1/simulation1_1_tuning.csv"
source("DataGeneration/DataGenerationFunctions.R")
source("DataGeneration/ExperimentFunctions.R")
source("Methods/Edgewise_GI-JOE.R")


args <- commandArgs(trailingOnly=TRUE);
graph_type <- args[[1]]
alpha <- as.numeric(args[[2]])
p <- as.numeric(args[[3]])
n1 <- as.numeric(args[[4]])
target_nodes <- c(2,4); n2 <- floor(n1 * alpha); 
#alpha is the ratio between n2 and n1

#generate the same graph for all replicates
model_seed <- 2021;
edge_val_range <- c(0.6,0.8)
eigmin <- 0.25
seed_vec <- 1 : 200


if(graph_type == "ER"){
  model_param <- precision_generation(p,edge_val_range,eigmin,model_seed,graph_type, prob = 3/(p - 1))
}else if(graph_type == "star"){
  model_param <- precision_generation(p,edge_val_range,eigmin,model_seed,graph_type, nstar = 3)
}else{
  model_param <- precision_generation(p,edge_val_range,eigmin,model_seed,graph_type)
}
Theta <- model_param[[1]]
Adj_mat <- model_param[[2]]
diag(Adj_mat) <- 0
Sigma <- solve(Theta)

#read previously selected tuning constant
tuning_results <- read.table(file = tuning_path, sep = " ", header = TRUE)
tuning_c <- tuning_results[tuning_results$graph == graph_type & tuning_results$p == p & 
                             tuning_results$n1 == n1 & tuning_results$alpha == alpha, c("tuning_c1", "tuning_c2")]
if(nrow(tuning_c) > 1){
  tuning_c <- as.numeric(tuning_c[1, ])
}else{
  tuning_c <- as.numeric(tuning_c)
}

#run experiments for seed in seed_vec
for(kk in 1:length(seed_vec)){
  seed <- seed_vec[kk]
  print(sprintf("running experiment with %s graph, p = %d, n1 = %d, alpha = %.2f, seed =%d", graph_type, p, n1, alpha, seed))
  test_tmp <- test_pairwise_experiment(p, target_nodes, n1, n2, seed, tuning_c = tuning_c, Theta = Theta, Sigma = Sigma, Adj_mat = Adj_mat, precisionCI = TRUE)
  result <- data.frame(graph = graph_type, p = p, n1 = n1, alpha = alpha, seed = seed, test_stat = test_tmp$test, var_est = test_tmp$var_est,
                       test_knownvar = test_tmp$test_knownvar)
  if(file.exists(csv_path)){
    write.table(result, file = csv_path, sep = " ", row.names=FALSE, append = TRUE, col.names = FALSE)
  }else{
    write.table(result, file = csv_path, sep = " ", row.names=FALSE)
  }
  coverage = (test_tmp$CI_Theta[1] <= Theta[target_nodes[1], target_nodes[2]] & test_tmp$CI_Theta[2] >= Theta[target_nodes[1], target_nodes[2]])
  result_precision <- data.frame(graph = graph_type, p = p, n1 = n1, alpha = alpha, n2 = n2, signal_ind = 1, seed = seed, coverage = coverage)
  if(file.exists(precision_csv_path)){
    write.table(result_precision, file = precision_csv_path, sep = " ", row.names=FALSE, append = TRUE, col.names = FALSE)
  }else{
    write.table(result_precision, file = precision_csv_path, sep = " ", row.names=FALSE)
  }
}
