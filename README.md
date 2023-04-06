# GI-JOE: Graph Inference when Joint Observations are Erose
This repository contains the R code for running GI-JOE method and reproducing the empirical results in the paper ``Graphical Model Inference with Erosely
Measured Data". To reproduce the results, the working directory of all code should be set as the home directory of this repository.

## Main Functions
The folder *Methods* contains the R code for running the GI-JOE method and running some baseline methods under the erose measurement settings.
- *Edgewise_GI-JOE.R* contains the functions for performing edge-wise testing, including code for stability selection
- *fullGraph_GI-JOE.R* contains the functions for performing whole graph testing, including the Bonferroni correction, Holm's correction, and the GI-JOE (FDR) method. Data-driven stability tuning code for the whole graph is also included.
- *baseline_graph_selection.R* contains the functions for running some baseline methods in the erose measurement setting, including the neighborhood lasso, graphical lasso, CLIME, debiased graphical lasso based on the projected positive semidefinite covariance estimate. The stability tuning code applicable for the erose measurements is also included.

## Simulations
The folder *Simulations* contains the R code for running and summarizing the simulations included in the paper. 
- *EdgewiseSimulation/sim1_1_main.R* runs the simulation for validating the Type I error of edge-wise testing. The graph structure, dimension, sample sizes can be input from the command line. 
- *EdgewiseSimulation/sim1_1_main.R* runs the simulation for validating the power of edge-wise testing.
- Both *sim1_1_main.R* and *sim1_1_main.R* read previously saved tuning result (*Results/sim1_1/simulation1_1_tuning.csv*, *Results/sim1_2/simulation1_2_tuning.csv*), which can be reproduced by running *sim1_1_tuning.R* and *sim1_2_tuning.R*
- After running *sim1_1_main.R* and *sim1_2_main.R* for the settings described in the manuscript, one would obtain *Results/sim1_1/simulation1_1_results.csv*, and *Results/sim1_2/simulation1_2_results.csv*, which can be summarized by *sim1_1_summary.R* and *sim1_2_summary.R*, to reproduce the figures in the manuscript.

- *FullGraphSimulations/sim2_main.R* runs the comparative studies with synthetic data and synthetic set-ups. The results are saved in *Results/sim2/F1_results.csv* and can be summarized by *sim2_summary.R*.
- *FullGraphSimulations/RealDataInspiredSimulations* include the codes for real data-inspired simulations, either with real neuroscience graphs or with real scRNA-seq measurement patterns. Both can be reproduced by running *main.R* and then *summary.R* in the corresponding subfolder.

- These simulation code used functions in the *Methods* and *DataGeneration* folder. 

## Real Data Experiments
The main code for real data experiments are included in *NeuroscienceDataAnalysis*. 
- *preprocessingFullData.R* preprocesses the original data file *raw_fl_traces.csv*, and save the result in *preprocessed.RData*
- *main.R* takes the preprocessed data and run our GI-JOE method and debiased graphical lasso. Results are saved in *AnalyzedResults.RData*
- *summary.R* summarizes the graph learning results and produce the figures and tables in the manuscript.
