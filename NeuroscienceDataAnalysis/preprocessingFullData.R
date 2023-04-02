
#estimation and testing for a certain threshold with the full raw data, raw traces, spontaneous period, testing against a threshold
source("DataGeneration/DataGenerationFunctions.R")
source("DataGeneration/ExperimentFunctions.R")
source("Methods/Edgewise_GI-JOE.R")
source("Methods/baseline_graph_selection.R")
source("Methods/FullGraph_GI-JOE.R")
library(huge)

##read data and preprocessing
data_path <- "NeuroscienceDataAnalysis/raw_fl_traces.csv"
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
#set a seed for reproducibility concern (the stability selection in huge is random)
set.seed(seed)
out_full <- huge(X_full, method = "glasso", nlambda = 20)
out_full_selection <- huge.select(out_full, criterion="stars")
Theta_full <- out_full_selection$opt.icov
Adj_mat_full <- Theta_full!=0; diag(Adj_mat_full) <- 0
p <- dim(X_full)[2]; n <- dim(X_full)[1]

# run our inference method on the full data set 
Obs <- matrix(rep(1, n * p), n, p)
thrs_c <- 5; signif_pval <- 0.05
start_t <- Sys.time()
nblasso_step1_results_full <- GIJOE_step1(X_full, Obs, n, p, thrs_c, signif_pval, onetuning = TRUE, symmetry = FALSE)
print(Sys.time() - start_t)

# variance estimation for each pair; loop over node a
var_est_full <- matrix(rep(0, p^2), p, p)
if(nblasso_step1_results_full$success){
  for(a in 1 : (p - 1)){
    var_est_full[a, (a + 1) : p] <- find_variances_GIJOE_sameN(nblasso_step1_results_full, n, a)
  }
}
graph_test_full <- GIJOE_thrs_step3(nblasso_step1_results_full, var_est_full, p, signif_pval, testing_thrs = 0.12)
save(X_full, Adj_mat_full, Theta_full, graph_test_full, nblasso_step1_results_full, var_est_full, 
     file = "NeuroscienceDataAnalysis/preprocessed.RData")


