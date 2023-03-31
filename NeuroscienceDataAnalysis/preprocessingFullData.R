
#estimation and testing for a certain threshold with the full raw data, raw traces, spontaneous period, testing against a threshold
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
source("sim2/step2/variance_estimate_methods.R")
source("sim2/step3/step3_functions.R")
library(huge)

##read data and preprocessing
data_path <- "/Applications/Files/research/graph-quilting/RealData/AllenBrainAtlasData/raw_fl_traces.csv"
raw_traces <- read.table(data_path, sep = ",")
raw_traces <- data.matrix(raw_traces[66938 : 75868, 2 : 228,])
nrow(raw_traces)
ncol(raw_traces)

##centralize and normalize the trace data
mean_traces <- matrix(rep(apply(raw_traces, 2, mean), nrow(raw_traces)), nrow(raw_traces), 
                      ncol(raw_traces), byrow = TRUE)
centralized_traces <- raw_traces - mean_traces
sd_traces <- matrix(rep(apply(raw_traces, 2, sd), nrow(raw_traces)), nrow(raw_traces), 
                    ncol(raw_traces), byrow = TRUE)
X_full <- centralized_traces / sd_traces
Corr_full <- t(X_full) %*% X_full / (nrow(X_full) - 1)

#fit a "ground truth" graph without low-rank adjustment for latent factors
seed <- 2059
set.seed(seed)
out_full <- huge(X_full, method = "glasso", nlambda = 20)
out_full$lambda
out_full_selection <- huge.select(out_full, criterion="stars")
out_full_selection$opt.lambda
out_full_selection$opt.sparsity
Theta_full <- out_full_selection$opt.icov
sum(Theta_full!=0)
Adj_mat_full <- Theta_full!=0; diag(Adj_mat_full) <- 0
p <- dim(X_full)[2]; n <- dim(X_full)[1]
sum(Adj_mat_full)/p

# run our inference method (variance estimation needs to be recoded so that it is not that slow!)
Obs <- matrix(rep(1, n * p), n, p)
thrs_c <- 5; signif_pval <- 0.05
start_t <- Sys.time()
nblasso_step1_results_full <- nblasso_step1(X_full, Obs, n, p, thrs_c, signif_pval, onetuning = TRUE, symmetry = FALSE)
print(Sys.time() - start_t)

# variance estimation for each pair; loop over node a
var_est_full <- matrix(rep(0, p^2), p, p)
if(nblasso_step1_results_full$success){
  for(a in 1 : (p - 1)){
    var_est_full[a, (a + 1) : p] <- find_variances_nblasso_sameN(nblasso_step1_results_full, n, a)
  }
}

graph_test_full <- find_graphs_nblasso_step3(nblasso_step1_results_full, var_est_full, p, signif_pval, testing_thrs = 0.1)
save(X_full, Adj_mat_full, Theta_full, graph_test_full, nblasso_step1_results_full, var_est_full, 
     file = sprintf("ABA_experiments/thrs_raw_preprocessed_seed%d.RData", seed))


