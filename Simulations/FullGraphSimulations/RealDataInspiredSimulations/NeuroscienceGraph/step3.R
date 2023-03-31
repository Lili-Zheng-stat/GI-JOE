## real data - inspired simulations: estimation and inference for the whole graph under different missing patterns; compare with baseline estimation and inference methods
## get graph type, measurement scenarios, total sample sizes, seed, batch of edges from the command line
## step 3: for our inference methods, find the p-values and tested graphs

wd <- "/Applications/Files/research/graph-quilting/graph-learning-varying-sample-size/simulations/inference_simulations"
setwd(wd)
mylibpath <- "/Users/lilybean/Downloads/software/R_packages"
.libPaths(mylibpath)

RData_step1_path <- "realdata_inspired_simulations/neuroscience/step1/RData/"
varest_path <- "realdata_inspired_simulations/neuroscience/step2/varest_results_updated/"
csv_path <- "realdata_inspired_simulations/neuroscience/updated_F1_results.csv";
preprocessing_path <- "ABA_experiments/preprocessed.RData"
setwd(wd)
source("fitting_functions_missing.R")
source("fitting_functions_glasso.R")
source("sim2_functions.R")
source("sim_functions.R")
source("baseline_fitting.R")
source("sim2/step3/step3_functions.R")


ms_sc_vec <- rep(1 : 3, each = 2);
n_vec <- c(8000,12000,5000,8000, 80000, 120000)
graph_type_vec <- rep("ABA", 6)

load(preprocessing_path)
p <- 227; 
Theta <- Theta_full;Theta[abs(Theta)<0.05] <- 0; Adj_mat <- Theta!=0; diag(Adj_mat) <- 0
Sigma <- solve(Theta)


summarize_F1_3steps_func <- function(graph_type, measurement_scenario, n_total, seed){
  screenSigma <- FALSE
  start_t <- Sys.time()
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
  return(data_exists)
}



for(i in 1:6){
  measurement_scenario <- ms_sc_vec[i]; n_total <- n_vec[i]; graph_type <- graph_type_vec[i]
  seed <- 1
  data_exists <- summarize_F1_3steps_func(graph_type, measurement_scenario, n_total, seed)
  while(data_exists){
    seed <- seed + 1
    data_exists <- summarize_F1_3steps_func(graph_type, measurement_scenario, n_total, seed)
  }
}

