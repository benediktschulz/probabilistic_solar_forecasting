# Aggregating distribution forecasts from deep ensembles

This repository provides R-code accompanying the paper

> Schulz, B. and Lerch, S. (2022). 
> Aggregating distribution forecasts from deep ensembles.
> Preprint available at https://arxiv.org/abs/2204.02291.

In particular, code for the implementation of the networks, aggregation methods, simulation study, case study and evaluation is available.

## Summary

We conduct a systematic comparison of aggregation methods for deep ensembles of distributional forecasts. Three network variants (Bernstein Quantile Network, Distributional Regression Network, Histogram Estimation Network) resulting in three different forecast types and two general aggregation approaches (Linear Pool, Vincentization) are considered. The methods are compared in a simulation study and a case study with wind gust predictions.


## Code

| File | Description |
| ---- | ----------- | 
| `fn_basic` | Helper functions for network functions. |
| `fn_eval` | Functions for the evaluation of probabilistic forecasts. |
| `fn_nn_cs` | Functions for the network variants in the case study. |
| `fn_nn_ss` | Functions for the network variants in the simulation study. |
| `figures` | Generation of the figures in the paper. |
| ---- | ----------- | 
| `cs_1_ensemble` | Case study: Deep ensemble generation. |
| `cs_2_aggregation` | Case study: Deep ensemble aggregation. |
| `cs_3_scores` | Case study: Evaluation of deep ensemble and aggregated forecasts. |
| ---- | ----------- | 
| `ss_0_data` | Simulation study: Data generation. |
| `ss_1_ensemble` | Simulation study: Deep ensemble generation. |
| `ss_2_aggregation` | Simulation study: Deep ensemble aggregation. |
| `ss_3_scores` | Simulation study: Evaluation of deep ensemble and aggregated forecasts. |
| ---- | ----------- | 
| `data/` | Directory for the evaluation data. |
| `plots/` | Directory for the plots generated by the `figure`-file. |

Note that the functions `fn_nn_cs` and `fn_nn_cs` mainly differ in the forecast distributions applied in the simulation and case study. The functions of the case study are tailored to the wind gust data, thus we recommend using the file of the simulation study to investigate the network variants, unless one is interested in a strictly positive predictive distribution or the station embedding applied in the case study.

## Data

### Simulation study

The simulated data, the deep ensembles and the aggregated forecasts result in files that are too large to be stored in this repository, thus we supply only the data on the evaluation of the forecasts. The files corresponding to the simulation study can be used to replicate the data, apply the network variants and aggregate the deep ensemble forecasts.

### Case study

The data was supplied by the German weather service (Deutscher Wetterdienst, DWD) and is not publicly available. For more information, we refer to the repository of the original study at https://github.com/benediktschulz/paper_pp_wind_gusts.

