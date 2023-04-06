##summarize real data inspired simulations results (real graph)
library(ggplot2)
library(dplyr)
library(kableExtra)
source("DataGeneration/ExperimentFunctions.R")

F1_results <- read.table(file = "Results/realdata_inspired_simulations/NeuroscienceGraph/F1_results.csv", sep = " ", header =  TRUE)
F1_results$F1 <- 2*F1_results$TP/(2*F1_results$TP+F1_results$FN+F1_results$FP)
F1_results$FDR <- F1_results$FP / (F1_results$FP + F1_results$TP)
F1_results$precision <- F1_results$TP / (F1_results$FP + F1_results$TP)
F1_results$recall <- F1_results$TP / (F1_results$FN + F1_results$TP)

F1_results_summary <- as.data.frame(F1_results%>%group_by(est_method, ms_sc, graph_type, n, seed)%>%
                                      summarise(precision = mean(precision), recall = mean(recall), F1 = mean(F1), 
                                                TP = mean(TP), FP = mean(FP), TN = mean(TN), FN = mean(FN)))

Est_methods <- c("Nlasso_JOE(AND)", "Nlasso_JOE(OR)", "Nlasso(AND)", "Nlasso(OR)", "Glasso", "CLIME")
Test_methods <- c("DBGlasso_Holm", "DBGlasso_FDR", "GIJOE_Holm", "GIJOE_FDR")
F1_results_summary$method_type <- rep(0, nrow(F1_results_summary))
for(i in 1 : nrow(F1_results_summary)){
  if(F1_results_summary$est_method[i] %in% Est_methods){
    F1_results_summary$method_type[i] <- "Estimation"
  }else if(F1_results_summary$est_method[i] %in% Test_methods){
    F1_results_summary$method_type[i] <- "Inference"
  }else{
    F1_results_summary$method_type[i] <- "Others"
  }
}
##summarize estimation results
F1_results_est <- F1_results_summary[F1_results_summary$method_type == "Estimation", ]

F1_results_test <- F1_results_summary[F1_results_summary$method_type == "Inference", ]

## table of TP, FP, FN, TN (mean and sd)
summarized_results_est <- as.data.frame(F1_results_est%>%group_by(est_method, ms_sc, graph_type, n)%>%
                                          summarise(TP_mean = mean(TP), TP_sd = sd(TP), FP_mean = mean(FP), FP_sd = sd(FP), 
                                                    TN_mean = mean(TN), TN_sd = sd(TN), FN_mean = mean(FN), FN_sd = sd(FN), 
                                                    TDR_mean = mean(precision), TDR_sd = sd(precision), 
                                                    TPR_mean = mean(recall), TPR_sd = sd(recall), 
                                                    TNR_mean = mean(TN / (TN + FP)), TNR_sd = sd(TN / (TN + FP)),
                                                    F1_mean = mean(F1), F1_sd = sd(F1)))
summarized_results_test <- as.data.frame(F1_results_test%>%group_by(est_method, ms_sc, graph_type, n)%>%
                                           summarise(TP_mean = mean(TP), TP_sd = sd(TP), FP_mean = mean(FP), FP_sd = sd(FP), 
                                                     TN_mean = mean(TN), TN_sd = sd(TN), FN_mean = mean(FN), FN_sd = sd(FN), 
                                                     TDR_mean = mean(precision), TDR_sd = sd(precision), 
                                                     TPR_mean = mean(recall), TPR_sd = sd(recall), 
                                                     TNR_mean = mean(TN / (TN + FP)), TNR_sd = sd(TN / (TN + FP)),
                                                     F1_mean = mean(F1), F1_sd = sd(F1)))
#recall is also true positive rate (TPR = TP / (TP + FN))
#FDP (1-precision)
#TNR = TN / (TN + FP)
graph_type <- "ABA"
n_mat <- matrix(c(8000, 12000, 5000, 8000, 80000, 120000),3,2,byrow = TRUE)

# choose a measurement scenario to produce table
ms_sc <- 1 # 2, 3

ii <- 1;
for(n in n_mat[ms_sc,]){
  if(ii == 1){
    table_tmp <- summarized_results_est[summarized_results_est$ms_sc == ms_sc & summarized_results_est$n == n & summarized_results_est$graph_type == graph_type, 
                                            c("est_method", "TPR_mean", "TPR_sd", "TNR_mean", "TNR_sd",
                                              "TDR_mean", "TDR_sd", "F1_mean", "F1_sd")] 
    for(j in 1:4){
      table_tmp[, (ncol(table_tmp) + 1)] <- paste(round(table_tmp[, 2*j], 3), " (", round(table_tmp[, (2*j+1)], 3), ")", sep = "")
    }
    colnames(table_tmp)[10 : 13] <- c("TPR", "TNR", "TDR", "F1")
    table_results <- table_tmp[order(table_tmp$est_method), c(1, 10 : 13)]
  }else{
    table_tmp <- summarized_results_est[summarized_results_est$ms_sc == ms_sc & summarized_results_est$n == n & summarized_results_est$graph_type == graph_type, 
                                          c("est_method", "TPR_mean", "TPR_sd", "TNR_mean", "TNR_sd",
                                            "TDR_mean", "TDR_sd", "F1_mean", "F1_sd")] 
    for(j in 1:4){
      table_tmp[, (ncol(table_tmp) + 1)] <- paste(round(table_tmp[, 2*j], 3), " (", round(table_tmp[, (2*j+1)], 3), ")", sep = "")
    }
    results_tmp <- table_tmp[order(table_tmp$est_method), c(1, 10 : 13)]
    table_results <- cbind(table_results, results_tmp[,2:ncol(results_tmp)])
    colnames(table_results)[(ncol(table_results) - 3) : ncol(table_results)] <- c("TPR", "TNR", "TDR", "F1")
  }
  ii <- ii + 1;
}
row.names(table_results) <- NULL
kbl(table_results) %>%
  kable_paper("striped", full_width = F) %>%
  add_header_above(c(" " = 1, "smaller n" = 4, "larger n" = 4))


ii <- 1;
for(n in n_mat[ms_sc,]){
  if(ii == 1){
    table_tmp <- summarized_results_test[summarized_results_test$ms_sc == ms_sc & summarized_results_test$n == n & summarized_results_test$graph_type == graph_type, 
                                             c("est_method", "TPR_mean", "TPR_sd", "TNR_mean", "TNR_sd",
                                               "TDR_mean", "TDR_sd", "F1_mean", "F1_sd")] 
    for(j in 1:4){
      table_tmp[, (ncol(table_tmp) + 1)] <- paste(round(table_tmp[, 2*j], 3), " (", round(table_tmp[, (2*j+1)], 3), ")", sep = "")
    }
    colnames(table_tmp)[10 : 13] <- c("TPR", "TNR", "TDR", "F1")
    table_results <- table_tmp[order(table_tmp$est_method), c(1, 10 : 13)]
  }else{
    table_tmp <- summarized_results_test[summarized_results_test$ms_sc == ms_sc & summarized_results_test$n == n & summarized_results_test$graph_type == graph_type, 
                                           c("est_method", "TPR_mean", "TPR_sd", "TNR_mean", "TNR_sd",
                                             "TDR_mean", "TDR_sd", "F1_mean", "F1_sd")] 
    for(j in 1:4){
      table_tmp[, (ncol(table_tmp) + 1)] <- paste(round(table_tmp[, 2*j], 3), " (", round(table_tmp[, (2*j+1)], 3), ")", sep = "")
    }
    results_tmp <- table_tmp[order(table_tmp$est_method), c(1, 10 : 13)]
    table_results <- cbind(table_results, results_tmp[,2:ncol(results_tmp)])
    colnames(table_results)[(ncol(table_results) - 3) : ncol(table_results)] <- c("TPR", "TNR", "TDR", "F1")
  }
  ii <- ii + 1;
}
row.names(table_results) <- NULL
kbl(table_results) %>%
  kable_paper("striped", full_width = F) %>%
  add_header_above(c(" " = 1, "smaller n" = 4, "larger n" = 4))




