##summarize simulation 2 results 
library(ggplot2)
library(dplyr)
library(xtable)
source("DataGeneration/ExperimentFunctions.R")

F1_results <- read.table(file = "Results/sim2/F1_results.csv", sep = " ", header =  TRUE)
F1_results$F1 <- 2*F1_results$TP/(2*F1_results$TP+F1_results$FN+F1_results$FP)
F1_results$FDR <- F1_results$FP / (F1_results$FP + F1_results$TP)
#precision is TDR
F1_results$precision <- F1_results$TP / (F1_results$FP + F1_results$TP)
#recall is TPR
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

n_mat <- matrix(c(600,800,1500,3000, 20000, 30000),3,2,byrow = TRUE)
ii <- 1;
for(graph_type in c("chain", "star", "ER")){
  for(ms_sc in c(1,3)){
    for(n in n_mat[ms_sc,]){
      if(ii == 1){
        table_results <- summarized_results_est[summarized_results_est$ms_sc == ms_sc & summarized_results_est$n == n & summarized_results_est$graph_type == graph_type, 
                                                c("est_method", "F1_mean", "F1_sd")] 
        table_results <- table_results[order(table_results$est_method),]
      }else{
        results_tmp <- summarized_results_est[summarized_results_est$ms_sc == ms_sc & summarized_results_est$n == n & summarized_results_est$graph_type == graph_type, 
                                              c("est_method", "F1_mean", "F1_sd")] 
        results_tmp <- results_tmp[order(results_tmp$est_method),]
        table_results <- cbind(table_results, results_tmp[,2:3])
      }
      ii <- ii + 1;
    }
  }
}
table_results

ii <- 1;
for(graph_type in c("chain", "star", "ER")){
  for(ms_sc in 1:3){
    for(n in n_mat[ms_sc,]){
      if(ii == 1){
        table_results <- summarized_results_test[summarized_results_test$ms_sc == ms_sc & summarized_results_test$n == n & summarized_results_test$graph_type == graph_type, 
                                                 c("est_method", "F1_mean", "F1_sd")] 
        table_results <- table_results[order(table_results$est_method),]
      }else{
        results_tmp <- summarized_results_test[summarized_results_test$ms_sc == ms_sc & summarized_results_test$n == n & summarized_results_test$graph_type == graph_type, 
                                               c("est_method", "F1_mean", "F1_sd")] 
        results_tmp <- results_tmp[order(results_tmp$est_method),]
        table_results <- cbind(table_results, results_tmp[,2:3])
      }
      ii <- ii + 1;
    }
  }
}
table_results

