##summarize simulation 1.1 (edge-wise type I error validation)
library(ggplot2)
library(dplyr)
library(latex2exp)
source("DataGeneration/ExperimentFunctions.R")

##read results for simulation 1.1
sim1_1_results <- read.table(file = "Results/sim1_1/simulation1_1_results.csv", sep = " ", header = TRUE)
sim1_1_results$p_val <- (1-pnorm(abs(sim1_1_results$test_stat)))*2
sim1_1_results$reject <- (sim1_1_results$p_val < 0.05)

sim1_1_summary <- sim1_1_results %>% group_by(graph, p, alpha, n1)
#type 1 error (1-coverage) vs. n1; facets depend on p, alpha; graph_type determines lines
sim1_1_type1 <- sim1_1_summary %>% summarise(rejection_rate = mean(reject), nrep = n())
sim1_1_type1$rejection_ci <- sqrt(sim1_1_type1$rejection_rate * (1-sim1_1_type1$rejection_rate)/sim1_1_type1$nrep)*qnorm(0.975)
p.labs <- c("p=50", "p=100", "p=200"); names(p.labs) <- c("50", "100", "200");
alpha.labs <- c("n2/n1=1", "n2/n1=1.2", "n2/n1=1.5")
names(alpha.labs) <- c("1", "1.2", "1.5")
g <- ggplot(sim1_1_type1, aes(x = n1, y = rejection_rate, col = graph)) + 
  geom_errorbar(aes(ymin = pmax(rejection_rate - rejection_ci, 0), 
                    ymax = pmin(rejection_rate + rejection_ci,1)), width=50) +
  labs(x = "n1", y = "type I error") + ylim(0, 0.25) + xlim(0, 1100) +
  geom_line() + facet_grid(alpha~p, labeller = labeller(alpha = alpha.labs, p = p.labs)) + 
  scale_color_discrete(name = "graph", labels = c("Chain", "Erdős–Rényi", "Star")) +
  theme(axis.text.x=element_text(size=10, angle=0, vjust=0.3),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=12),axis.title.y=element_text(size=12),
        plot.title=element_text(size=11), legend.title=element_text(size=15), legend.text=element_text(size=15))
g <- g + geom_hline(yintercept=0.05, linetype="dashed", color = "black")
g



