## real data - inspired simulations: estimation and inference for the whole graph under different missing patterns; compare with baseline estimation and inference methods
## get graph type, measurement scenarios, total sample sizes, seed, batch of edges from the command line
## step 2: for our inference methods, finding the variance estimation for the batch edges 

.libPaths("/home/liliz/R_packages")
wd <- "/home/liliz/Graph_DiffN/simulations"
RData_step1_path <- "ABA_simulation/step1/RData/"
csv_path <- "ABA_simulation/step2/varest_results/";
setwd(wd)
source("fitting_functions_missing.R")
source("fitting_functions_glasso.R")
source("sim2_functions.R")
source("sim_functions.R")
source("baseline_fitting.R")
source("sim2/step2/variance_estimate_methods.R")


args <- commandArgs(trailingOnly=TRUE);
measurement_scenario <- as.numeric(args[[1]])
n_total <- as.numeric(args[[2]])
seed <- as.numeric(args[[3]])
batch_edges <- as.numeric(args[[4]])

p <- 227; batch_size <- p; graph_type <- "ABA"

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

