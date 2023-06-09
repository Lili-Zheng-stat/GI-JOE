##summarize real data inspired simulations results (real measurement patterns)
library(ggplot2)
library(dplyr)
library(kableExtra)
source("DataGeneration/ExperimentFunctions.R")

F1_results <- read.table(file = "Results/realdata_inspired_simulations/scRNA-seqMeasurements/F1_results.csv", sep = " ", header =  TRUE)
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
ii <- 1;
for(ms_sc in c(1,2)){
  for(graph_type in c("scale-free", "small-world")){
    if(ii == 1){
      table_tmp <- summarized_results_est[summarized_results_est$ms_sc == ms_sc & summarized_results_est$graph_type == graph_type, 
                                              c("est_method", "F1_mean", "F1_sd")] 
      table_tmp[,4] <- paste(round(table_tmp[,2], 3), " (", round(table_tmp[,3], 3), ")", sep = "")
      colnames(table_tmp)[4] <- graph_type
      table_results <- table_tmp[order(table_tmp$est_method),c(1,4)]
    }else{
      table_tmp <- summarized_results_est[summarized_results_est$ms_sc == ms_sc & summarized_results_est$graph_type == graph_type, 
                                            c("est_method", "F1_mean", "F1_sd")] 
      table_tmp[,4] <- paste(round(table_tmp[,2], 3), " (", round(table_tmp[,3], 3), ")", sep = "")
      results_tmp <- table_tmp[order(table_tmp$est_method),c(1,4)]
      table_results <- cbind(table_results, results_tmp[,2])
      colnames(table_results)[ncol(table_results)] <- graph_type
    }
    ii <- ii + 1;
  }
}
row.names(table_results) <- NULL
kbl(table_results) %>%
  kable_paper("striped", full_width = F) %>%
  add_header_above(c(" " = 1, "chu measurement" = 2, "darmanis measurement" = 2))

ii <- 1;
for(ms_sc in c(1,2)){
  for(graph_type in c("scale-free", "small-world")){
    if(ii == 1){
      table_tmp <- summarized_results_test[summarized_results_test$ms_sc == ms_sc & summarized_results_test$graph_type == graph_type, 
                                               c("est_method", "F1_mean", "F1_sd")] 
      table_tmp[,4] <- paste(round(table_tmp[,2], 3), " (", round(table_tmp[,3], 3), ")", sep = "")
      colnames(table_tmp)[4] <- graph_type
      table_results <- table_tmp[order(table_tmp$est_method),c(1,4)]
    }else{
      table_tmp <- summarized_results_test[summarized_results_test$ms_sc == ms_sc & summarized_results_test$graph_type == graph_type, 
                                             c("est_method", "F1_mean", "F1_sd")] 
      table_tmp[,4] <- paste(round(table_tmp[,2], 3), " (", round(table_tmp[,3], 3), ")", sep = "")
      results_tmp <- table_tmp[order(table_tmp$est_method),c(1,4)]
      table_results <- cbind(table_results, results_tmp[,2])
      colnames(table_results)[ncol(table_results)] <- graph_type
    }
    ii <- ii + 1;
  }
}
row.names(table_results) <- NULL
kbl(table_results) %>%
  kable_paper("striped", full_width = F) %>%
  add_header_above(c(" " = 1, "chu measurement" = 2, "darmanis measurement" = 2))



