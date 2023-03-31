## real data - inspired simulations (real measurement patterns): estimation and inference for the whole graph under different missing patterns; compare with baseline estimation and inference methods
## get graph type, measurement scenarios, total sample sizes, seed, batch of edges from the command line
## step 3: for our inference methods, find the p-values and tested graphs

wd <- "/Applications/Files/research/graph-quilting/graph-learning-varying-sample-size/simulations/inference_simulations"
setwd(wd)
mylibpath <- "/Users/lilybean/Downloads/software/R_packages"
.libPaths(mylibpath)
library(igraph)

RData_step1_path <- "realdata_inspired_simulations/scRNA-seq/step1/RData/"
varest_path <- "realdata_inspired_simulations/scRNA-seq/step2/varest_results_updated/"
csv_path <- "realdata_inspired_simulations/scRNA-seq/Updated_F1_results.csv";
Obs_path <- "realdata_inspired_simulations/scRNA-seq/scRNAseq_Obs.RData"
setwd(wd)
source("fitting_functions_missing.R")
source("fitting_functions_glasso.R")
source("sim2_functions.R")
source("sim_functions.R")
source("baseline_fitting.R")
source("sim2/step3/step3_functions.R")


ms_sc_vec <- rep(1 : 3, each = 4);
graph_type_vec <- rep(c("small-world", "scale-free", "chain", "star"), 3)

load(Obs_path)
p <- 200
model_seed <- 2021;
edge_val_range <- c(0.6,0.8)
eigmin <- 0.25

summarize_F1_3steps_func <- function(graph_type, measurement_setting, seed){
  screenSigma <- FALSE
  start_t <- Sys.time()
  data_exists <- FALSE
  
  #generate model
  if(graph_type == "star"){
    model_param <- precision_generation(p,edge_val_range,eigmin,model_seed,graph_type, nstar = 10)
  }else if(graph_type == "ER"){
    model_param <- precision_generation(p,edge_val_range,eigmin,model_seed,graph_type, prob = 3 / (p-1))
  }else if(graph_type == "small-world"){
    model_param <- precision_generation(p,edge_val_range,eigmin,model_seed,graph_type, degree = 2)
  }else if(graph_type == "scale-free"){
    model_param <- precision_generation(p,edge_val_range,eigmin,model_seed,graph_type, new_edge = 1)
    Theta <- model_param[[1]][p:1,p:1]
    Adj_mat <- model_param[[2]][p:1,p:1]
  }else{
    model_param <- precision_generation(p,edge_val_range,eigmin,model_seed,graph_type)
  }
  if(graph_type != "scale-free"){
    Theta <- model_param[[1]]
    Adj_mat <- model_param[[2]]
  }
  
  diag(Adj_mat) <- 0
  
  RDatafilename <- paste(RData_step1_path, sprintf("NDS_%d%sexpr%d.RData", measurement_setting, graph_type, seed), sep = "")
  if(file.exists(RDatafilename)){
    load(RDatafilename, verbose = TRUE)
    N <- t(Obs)%*%Obs
    varest_file_name <- paste(varest_path, sprintf("var_est_%d%sexpr%dbatch1.csv", measurement_setting, graph_type, seed), sep = "")
    if(file.exists(varest_file_name)){
      data_exists <- TRUE
      var_est_results <- read.table(varest_file_name, sep = " ", header = TRUE)
      p_kept <- length(nblasso_step1_results$kept_nodes); 
      for(batch in 2 : ceiling(p_kept * (p_kept - 1) / 2 / p)){
        varest_file_name <- paste(varest_path, sprintf("var_est_%d%sexpr%dbatch%d.csv", measurement_setting, graph_type, seed, batch), sep = "")
        if(file.exists(varest_file_name)){
          var_est_results_add <- read.table(varest_file_name, sep = " ", header = TRUE)
          var_est_results <- rbind(var_est_results, var_est_results_add)
        }else{
          data_exists <- FALSE
          break;
        }
      }
    }
  }
  if(data_exists){
    if(nblasso_step1_results$success){
      #our_graphs <- find_graphs_nblasso_step3(nblasso_step1_results, var_est_results, p, signif_pval = 0.05)
      record_graphs_nblasso_step3(nblasso_step1_results, var_est_results, p, signif_pval = 0.05, screenSigma, graph_type, 
                                  measurement_setting, nrow(Obs), seed, Adj_mat, csv_path)
    }
    print(Sys.time()-start_t)
    print(sprintf("Job completed for %s!", RDatafilename))
  }
  return(data_exists)
}



for(i in 2:12){
  measurement_setting <- ms_sc_vec[i]; graph_type <- graph_type_vec[i]
  seed <- 1
  data_exists <- summarize_F1_3steps_func(graph_type, measurement_setting, seed)
  while(data_exists){
    seed <- seed + 1
    data_exists <- summarize_F1_3steps_func(graph_type, measurement_setting, seed)
  }
}
