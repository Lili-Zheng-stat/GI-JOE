##summarize precision inference results 
library(ggplot2)
library(dplyr)
library(latex2exp)
source("DataGeneration/ExperimentFunctions.R")
csv_path <- "Results/precision_CI_results.csv"

simTheta_results <- read.table(file = csv_path, sep = " ", header = TRUE)
simTheta_summary <- simTheta_results %>% group_by(graph, signal_ind, n1)


#power vs. signal_strength; facets depend on n2, alpha; graph_type determines lines
simTheta_summary <- simTheta_summary %>% summarise(coverage = mean(coverage), nrep = n())

simTheta_summary$coverage_ci <- sqrt(simTheta_summary$coverage * (1-simTheta_summary$coverage)/simTheta_summary$nrep)*qnorm(0.975)

signal.labs <- c("signal = 0", "signal = 2.2", "signal = 4.5"); names(signal.labs) <- c("1", "2", "3");
g <- ggplot(simTheta_summary, aes(x = n1, y = coverage, col = graph)) + 
  geom_errorbar(aes(ymin = pmax(coverage - coverage_ci, 0), 
                    ymax = pmin(coverage + coverage_ci,1)), width=50) +
  labs(x = "n1", y = "Coverage rate") + ylim(0.7, 1) + xlim(0, 1100) +
  geom_line() + facet_wrap(~signal_ind, labeller = labeller(signal_ind = signal.labs)) + 
  scale_color_discrete(name = "graph", labels = c("Chain", "Erdős–Rényi", "Star")) +
  theme(axis.text.x=element_text(size=10, angle=0, vjust=0.3),
        axis.text.y=element_text(size=12),
        axis.title.x=element_text(size=12),axis.title.y=element_text(size=12),
        plot.title=element_text(size=11), legend.title=element_text(size=15), legend.text=element_text(size=15))
g <- g + geom_hline(yintercept=0.95, linetype="dashed", color = "black")
g
