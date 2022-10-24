# Probabilistic Solar Forecasting

This repository provides R-code accompanying the paper

> Gneiting, T., Lerch, S. and Schulz, B. (2022). 
> Probabilistic Solar Forecasting: Benchmarks, Post-processing, Verification.
> TBD-TBD

In particular, code for the implementation of the postprocessing methods, case study and evaluation is available.

## Summary

Here, we provide code for the replication of the results in the above-mentioned paper. The setting and data was adopted from Yang et al. (2022), who theirselves provide data and code for the replication of their work. In the case study, probabilistic one day-ahead forecasts of Global Horizontal Irradiance (GHI) are issued for each hour of the following day at seven locations in the continental US. Five predictor variables are given, of which we highlight two. One is a deterministic GHI forecast from a numerical weather prediction model is (ECMWF HRES), the other is the REST2 clear-sky irradiance. For further information on the data and setting, we refer to Yang et al. (2022).

We use four methods for postprocessing of solar radiation forecasts, namely the Analogue Ensemble (AnEn), Isotonic Distributional Regression (IDR), the Distributional Regression Network (DRN) and the Bernstein Quantile Network (BQN). While AnEn is replicated from Yang et al. (2022), the other three methods are based on the code of Schulz and Lerch (2022), who use these methods for probabilistic wind gust forecasting. For the evaluation of the methods, we again rely on code from Schulz and Lerch (2022) and for the generation of quantile reliability diagrams also on code from Gneiting et al. (2023).


## Code

This repository includes three directories. 'code' includes the R-code that can be used for replication, 'data' the corresponding data and 'figures' the figures in the paper.

The following table lists the scripts provided in the 'code'-directory:

| File | Description |
| ---- | ----------- |  
| `functions_basic` | Basic functions (based on code of Schulz and Lerch, 2022). |
| `functions_data` | Functions for handling of data (based on code of Schulz and Lerch, 2022). |
| `functions_pp` | Implementation of the postprocessing methods (based on code of Schulz and Lerch, 2022). |
| `functions_eval` | Functions for evaluation of the forecasts (based on code of Schulz and Lerch, 2022). |
| `functions_quantile_rd` | Functions for quantile reliability diagrams (based on code of Gneiting et al., 2023). |
| ---- | ----------- | 
| `data_preprocessing` | Preprocessing of the data analogous to Yang et al. (2022). |
| ---- | ----------- | 
| `pp_AnEn` | Implementation of the AnEn-method taken from Yang et al. (2022). |
| `pp_idr` | Postprocessing via IDR. |
| `pp_drn` | Postprocessing via DRN. |
| `pp_bqn` | Postprocessing via BQN. |
| ---- | ----------- | 
| `evaluation_scores` | Summary of the evaluation measures for the postprocessing methods. |
| `evaluation_quantile_data` | Calculate data needed for quantile reliability diagrams (based on code of Gneiting et al., 2023). |
| `figures_paper` | Generation of figures from the paper. |
| ---- | ----------- |


## Comments on the implementation of the postprocessing methods

# DRN

For DRN, we do not forecast GHI directly but instead the bias of the ECMWF HRES forecast (forecast minus observation), which we model with a normal distribution. The GHI forecast is then obtained by subtracting the location parameter of the forecast distribution from the ECMWF HRES forecast. However, this does not ensure a positive forecast. Therefore, the forecast distribution is truncated in zero. The loss function used is the CRPS, the configuration is described later.

We also provide an implementation of DRN that directly models GHI based on a truncated normal or logistic distribution.

# BQN

As for DRN, we do not forecast GHI directly but instead the bias of the ECMWF HRES forecast. Here, we then have to substract the quantile function from the ECMWF HRES forecast, which also includes a change of orientation of the quantile forecast. After wards, we censor the forecasts in zero to obtain non-negativity. The loss function used to estimate the network parameters is the mean pinball loss over 99 equidistant quantiles on the unit interval, i.e. at steps of 1% at the levels 1%,..., 99%.

The PIT values of the BQN forecasts cannot be calculated directly, thus we rely on a generalization of PIT values, the uPIT (as in Schulz and Lerch, 2022). The evaluation of the BQN forecasts is based on a set of 99 equidistant quantiles for those measures that cannot be calculated exactly (CRPS, Brier score, CORP reliability diagrams).

We also provide an implementation of BQN that directly models GHI.

# IDR

As AnEN, IDR is estimated separately for eachs station.

As mentioned in the paper, we here compare two variants of IDR. CSD-IDR was used in the paper and is based on the clear-sky model REST2, while GHI-IDR is based on the ECMWF HRES forecast. 

| CRPS | BON | DRA | FPK | GWN | PSU | SXF | TBL |
| ---- | --- | --- |  
| CSD-IDR | 82.9 | 37.9 | 62.4 | 85.6 | 84.6 | 74.8 | 74.2 | 
| GHI-IDR | 58.4 | 32.8 | 50.6 | 62.7 | 62.7 | 57.8 | 61.2 | 
| ---- | --- | --- |   
| Brier Score (250 W/m^2) | BON | DRA | FPK | GWN | PSU | SXF | TBL |
| ---- | --- | --- |  
| CSD-IDR | 0.141 | 0.040 | 0.107 | 0.133 | 0.142 | 0.124 | 0.104 | 
| GHI-IDR | 0.097 | 0.033 | 0.084 | 0.091 | 0.092 | 0.090 | 0.083 | 
| ---- | --- | --- |   
| PI Length | BON | DRA | FPK | GWN | PSU | SXF | TBL |
| ---- | --- | --- |  
| CSD-IDR | 387.8 | 233.0 | 324.6 | 414.0 | 401.6 | 386.6 | 398.4 | 
| GHI-IDR | 290.7 | 184.1 | 264.9 | 308.2 | 301.2 | 295.8 | 330.6 | 
| ---- | --- | --- |   
| Pinball Loss (8.33%) | BON | DRA | FPK | GWN | PSU | SXF | TBL |
| ---- | --- | --- |  
| CSD-IDR | 23.9 | 21.0 | 22.1 | 24.7 | 22.1 | 24.1 | 26.1 | 
| GHI-IDR | 18.7 | 15.0 | 16.8 | 19.8 | 17.1 | 19.9 | 21.3 | 
| ---- | --- | --- |   
| Pinball Loss (91.67%) | BON | DRA | FPK | GWN | PSU | SXF | TBL |
| ---- | --- | --- |  
| CSD-IDR | 13.5 | 6.1 | 10.0 | 13.9 | 15.6 | 11.8 | 11.9 | 
| GHI-IDR | 17.7 | 10.6 | 15.8 | 19.7 | 21.4 | 17.1 | 16.1 | 
| ---- | --- | --- |   

## Network hyperparameters

The following configurations have been used for the network methods. Both network methods use station embedding for a global estimation (in contrast to AnEn and IDR that are estimated locally) and have been implemented using the Keras and Tensorflow-Interfaces to R. For further details, we refer to Schulz and Lerch (2022).

| Hyperparmaeter | DRN | BQN |
| ---- | --- | --- |  
| Learning rate | 5 x 10^-4 | 5 x 10^-4 |
| Epochs | 150 | 150 |
| Patience | 10 | 10 |
| Batch size | 32 | 32 |
| Embedding dimension | 5 | 5 |
| Hidden layers | 3 | 3 |
| Nodes per layer | (96, 48, 24) | (96, 48, 24) |
| Activation | Softplus | Softplus |
| Output nodes | 2 | 13 |
| Output activation | (Linear, Sofplus) | (Linear, 12xSoftplus) |
| Size of network ensemble | 10 | 10 |
| ---- | --- | --- |  

## References

- Gneiting, T., Wolffram, D., Resin, J., Kraus, K., Bracher, J., Dimitriadis, T., Hagenmeyer, V., Jordan, A.I., Lerch, S., Phipps, K., Schienle, M., 2023. Model diagnostics and forecast evaluation for quantiles. Annual Review of Statistics and Its Application 10, TBD–TBD.
- Schulz, B., Lerch, S., 2022. Machine learning methods for postprocessing ensemble forecasts of wind gusts: A systematic comparison. Monthly Weather Review 150, 235–257.
- Yang, D., Wang, W., Hong, T., 2022. A historical weather forecast dataset from the European Centre for Medium-Range Weather Forecasts (ECMWF) for energy forecasting. Solar Energy 232, 263–274.
