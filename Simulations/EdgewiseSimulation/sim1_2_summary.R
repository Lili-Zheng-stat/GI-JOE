##summarize simulation 1.2 (edge-wise power validation)
library(ggplot2)
library(dplyr)
library(latex2exp)
##read data for simulation 1.1
source("DataGeneration/ExperimentFunctions.R")

sim1_2_results <- read.table(file = "Results/sim1_2/simulation1_2_results.csv", sep = " ", header = TRUE)
sim1_2_results$snr <- sim1_2_results$truesignal/sqrt(sim1_2_results$var_true)
sim1_2_results$p_val <- (1-pnorm(abs(sim1_2_results$test_stat)))*2
sim1_2_results$reject <- (sim1_2_results$p_val < 0.05)

sim1_2_summary <- sim1_2_results %>% group_by(graph, signal_ind, alpha, n2)
#power vs. signal_strength; facets depend on n2, alpha; graph_type determines lines
sim1_2_power <- sim1_2_summary %>% summarise(rejection_rate = mean(reject), snr = mean(snr), nrep = n())
sim1_2_power$rejection_ci <- sqrt(sim1_2_power$rejection_rate * (1-sim1_2_power$rejection_rate)/sim1_2_power$nrep)*qnorm(0.975)
n.labs <- c("n2 = 125", "n2 = 250", "n2 = 500"); names(n.labs) <- c("125", "250", "500");
alpha.labs <- c("n2/n1 = 1", "n2/n1 = 1.2", "n2/n1=1.5"); 
names(alpha.labs) <- c("1", "1.2", "1.5")
g <- ggplot(sim1_2_power, aes(x = snr, y = rejection_rate, col = as.factor(n2))) + 
  geom_errorbar(aes(ymin = pmax(rejection_rate - rejection_ci, 0), ymax = pmin(rejection_rate + rejection_ci,1)), width=.5, position=position_dodge(0.05)) +
  labs(x = TeX("$\\frac{\\Theta^*_{a,b}}{\\Theta^*_{a,a}\\sigma_n(a,b)}$"), y = "power") + ylim(0, 1) +
  geom_line() + facet_grid(alpha~graph, labeller = labeller(alpha = alpha.labs, n2 = n.labs)) + scale_color_discrete(name = "n2", labels = c("125", "250", "500")) +
  theme(axis.text.x=element_text(size=10, angle=0, vjust=0.3),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=12),axis.title.y=element_text(size=12),
        plot.title=element_text(size=11), legend.title=element_text(size=15), legend.text=element_text(size=15))
g




