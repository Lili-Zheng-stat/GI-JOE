## Updated simulation 2 with new variance estimate: estimation and inference for the whole graph under different missing patterns; compare with baseline estimation and inference methods
## get graph type, measurement scenarios, total sample sizes, seed, batch of edges from the command line
## step 3: for our inference methods, find the p-values and tested graphs

wd <- "/Applications/Files/research/graph-quilting/graph-learning-varying-sample-size/simulations/inference_simulations"
setwd(wd)
mylibpath <- "/Users/lilybean/Downloads/software/R_packages"
.libPaths(mylibpath)

RData_step1_path <- "sim2/step1/RData/"
varest_path <- "sim2/step2/varest_results_updated/"
csv_path <- "sim2/updated_F1_results.csv";
setwd(wd)
source("fitting_functions_missing.R")
source("fitting_functions_glasso.R")
source("sim2_functions.R")
source("sim_functions.R")
source("baseline_fitting.R")
source("sim2/step3/step3_functions.R")
library(huge)

ms_sc_vec <- rep(1 : 3, each = 6);
n_vec <- rep(c(600,800,1500,3000, 20000, 30000), each = 3)
graph_type_vec <- rep(c("chain", "star", "ER"), 6)

summarize_F1_3steps_func <- function(graph_type, measurement_scenario, n_total, seed){
  p <- 200; screenSigma <- FALSE
  model_seed <- 2021;
  edge_val_range <- c(0.6,0.8)
  eigmin <- 0.25
  
  start_t <- Sys.time()
  #generate model
  if(graph_type == "star"){
    model_param <- precision_generation(p,edge_val_range,eigmin,model_seed,graph_type, nstar = 10)
  }else if(graph_type == "ER"){
    model_param <- precision_generation(p,edge_val_range,eigmin,model_seed,graph_type, prob = 5 / p)
  }else{
    model_param <- precision_generation(p,edge_val_range,eigmin,model_seed,graph_type)
  }
  
  Theta <- model_param[[1]]
  Adj_mat <- model_param[[2]]
  diag(Adj_mat) <- 0
  data_exists <- FALSE
  RDatafilename <- paste(RData_step1_path, sprintf("NDS_%s%dn%dexpr%d.RData", graph_type, measurement_scenario, n_total, seed), sep = "")
  if(file.exists(RDatafilename)){
    load(RDatafilename)
    varest_file_name <- paste(varest_path, sprintf("var_est_%s%dn%dexpr%dbatch1.csv", graph_type, measurement_scenario, n_total, seed), sep = "")
    if(file.exists(varest_file_name)){
      data_exists <- TRUE
      var_est_results <- read.table(varest_file_name, sep = " ", header = TRUE)
      p_kept <- length(nblasso_step1_results$kept_nodes); 
      for(batch in 2 : ceiling(p_kept * (p_kept - 1) / 2 / p)){
        varest_file_name <- paste(varest_path, sprintf("var_est_%s%dn%dexpr%dbatch%d.csv", graph_type, measurement_scenario, n_total, seed, batch), sep = "")
        if(file.exists(varest_file_name)){
          var_est_results_add <- read.table(varest_file_name, sep = " ", header = TRUE)
          var_est_results <- rbind(var_est_results, var_est_results_add)
        }else{
          data_exists <- FALSE
          break;
        }
      }
    }else{
      varest_file_name <- paste(varest_path, sprintf("var_est_%s%dn%dexpr%d.csv", graph_type, measurement_scenario, n_total, seed), sep = "")
      if(file.exists(varest_file_name)){
        data_exists <- TRUE
        var_est_results <- read.table(varest_file_name, sep = " ", header = TRUE)
      }
    }
    
    if(data_exists){
      if(nblasso_step1_results$success){
        record_graphs_nblasso_step3(nblasso_step1_results, var_est_results, p, signif_pval = 0.05, screenSigma, graph_type, 
                                    measurement_scenario, n_total, seed, Adj_mat, csv_path)
      }
      print(Sys.time()-start_t)
      print(sprintf("Job completed for %s!", RDatafilename))
    }
  }
  return(data_exists)
}

for(i in 1:18){
  measurement_scenario <- ms_sc_vec[i]; n_total <- n_vec[i]; graph_type <- graph_type_vec[i]
  seed <- 1
  data_exists <- summarize_F1_3steps_func(graph_type, measurement_scenario, n_total, seed)
  while(data_exists){
    seed <- seed + 1
    data_exists <- summarize_F1_3steps_func(graph_type, measurement_scenario, n_total, seed)
  }
}


results <- read.table(file = csv_path, sep = " ", header =  TRUE)
unique(results$est_method)
results$F1[results$est_method=="nblasso_signif_FDR"]
results[results$est_method=="nblasso_signif_FDR",][results$F1[results$est_method=="nblasso_signif_FDR"] <= 0.7,]

