## simulation 2: estimation and inference for the whole graph under different missing patterns; compare with baseline estimation and inference methods
## get graph type, measurement scenarios, total sample sizes, seed, batch of edges from the command line
## step 2: for our inference methods, finding the variance estimation for the batch edges 

.libPaths("/home/liliz/R_packages")
wd <- "/home/liliz/Graph_DiffN/simulations"
RData_step1_path <- "sim2/step1/RData/"
csv_path <- "sim2/step2/varest_results/"
setwd(wd)
source("fitting_functions_missing.R")
source("fitting_functions_glasso.R")
source("sim2_functions.R")
source("sim_functions.R")
source("baseline_fitting.R")
source("sim2/step2/variance_estimate_methods.R")


args <- commandArgs(trailingOnly=TRUE);
graph_type <- args[[1]]
measurement_scenario <- as.numeric(args[[2]])
n_total <- as.numeric(args[[3]])
seed <- as.numeric(args[[4]])
batch_edges <- as.numeric(args[[5]])

p <- 200; batch_size <- p;
#model_seed <- 2021;
#edge_val_range <- c(0.6,0.8)
#eigmin <- 0.25

#generate model
#if(graph_type == "star"){
#  model_param <- precision_generation(p,edge_val_range,eigmin,model_seed,graph_type, nstar = 10)
#}else if(graph_type == "ER"){
#  model_param <- precision_generation(p,edge_val_range,eigmin,model_seed,graph_type, prob = 5 / p)
#}else{
#  model_param <- precision_generation(p,edge_val_range,eigmin,model_seed,graph_type)
#}

#Theta <- model_param[[1]]
#Adj_mat <- model_param[[2]]
#diag(Adj_mat) <- 0
#Sigma <- solve(Theta)

#measurement scenario
#Obs_seed <- 1000
#Obs <- generate_Obs(p, n_total, Adj_mat, measurement_scenario, Obs_seed)

load(file = paste(RData_step1_path, sprintf("NDS_%s%dn%dexpr%d.RData", graph_type, measurement_scenario, n_total, seed), sep = ""))
p_kept <- length(nblasso_step1_results$kept_nodes)
row_ind_set <- NULL; for(i in 1 : (p_kept - 1)){row_ind_set <- c(row_ind_set, rep(i, p_kept - i))}
col_ind_set <- NULL; for(i in 2 : p_kept){col_ind_set <- c(col_ind_set, i : p_kept)}
if((batch_edges - 1) * batch_size + 1 <= length(row_ind_set)){
  if(batch_edges * batch_size < length(row_ind_set)){
    row_ind <- row_ind_set[((batch_edges - 1) * batch_size + 1) : (batch_edges * batch_size)]; 
    col_ind <- col_ind_set[((batch_edges - 1) * batch_size + 1) : (batch_edges * batch_size)];
  }else{
    row_ind <- row_ind_set[((batch_edges - 1) * batch_size + 1) : length(row_ind_set)]; 
    col_ind <- col_ind_set[((batch_edges - 1) * batch_size + 1) : length(row_ind_set)];
  }

#graph_type <- "chain"; measurement_scenario <- 4; n_total <- 20000; seed <- 1

  if(nblasso_step1_results$success){
    nblasso_var_est <- find_variances_nblasso(nblasso_step1_results, Obs, row_ind, col_ind)
    results <- data.frame(rows = nblasso_step1_results$kept_nodes[row_ind], cols =  nblasso_step1_results$kept_nodes[col_ind], var_est = nblasso_var_est, method = "nblasso")
    file_name <- paste(csv_path, sprintf("var_est_%s%dn%dexpr%dbatch%d.csv", graph_type, measurement_scenario, n_total, seed, batch_edges), sep = "")
    write.table(results, file = file_name, sep = " ", row.names=FALSE)
  }
}

