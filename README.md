# Probabilistic Solar Forecasting

This repository provides R code and data accompanying the paper

> Gneiting, T., Lerch, S. and Schulz, B. (2022). 
> Probabilistic Solar Forecasting: Benchmarks, Post-processing, Verification.
> TBD-TBD

In particular, code for the implementation of the postprocessing methods, application to the case study and forecast evaluation is available.

## Summary

The setting of the case study and the dataset were adopted from Yang et al. (2022), who themselves provide data and code for the replication of their work. In the case study, probabilistic one day-ahead forecasts of Global Horizontal Irradiance (GHI) are issued for each hour of the following day at seven locations in the continental US. Five predictor variables are available, of which we highlight two. One is a deterministic GHI forecast from a numerical weather prediction model (ECMWF HRES), the other is the REST2 clear-sky irradiance. For further information on the data and setting, we refer to Yang et al. (2022). The data and code are provided the supplemental material of their work that is accessible via the DOI noted in the references.

We use four methods for postprocessing of solar radiation forecasts, namely the Analogue Ensemble (AnEn), Isotonic Distributional Regression (IDR), the Distributional Regression Network (DRN) and the Bernstein Quantile Network (BQN). While AnEn is replicated from Yang et al. (2022), the other three methods are based on the code of Schulz and Lerch (2022), who adapt and apply these methods to probabilistic wind gust forecasting. For the evaluation of the probabilistic forecasts, we again build on code from Schulz and Lerch (2022) and for the generation of quantile reliability diagrams on code from Gneiting et al. (2023).


## Data

We are grateful to Dazhi Yang for generous advice on the handling of the benchmark data from Yang et al. (2022) as well as the permission to redistribute the data.

The original data provided in `/data/original_data/` is described in Appendix B of Yang et al. (2022) and also available at the supplemental material at https://doi.org/10.1016/j.solener.2021.12.011 under a 'CC BY 4.0' license. The data sets `data_bon`,..., `data_tbl` as well as `data_total` in the `/data/` directory have been generated via the `data_preprocessing` script for preprocessing of the original data.


## Code

This repository includes three directories. `/code/` includes the R code to replicate our results, `/data/` the corresponding data and `/figures/` the figures in the paper.

The following table lists the scripts provided in the `/code/`-directory:

| File | Description |
| ---- | ----------- |  
| `functions_basic` | Basic functions (based on code from Schulz and Lerch, 2022). |
| `functions_pp` | Implementation of the postprocessing methods (based on code from Schulz and Lerch, 2022). |
| `functions_eval` | Functions for evaluation of the forecasts (based on code from Schulz and Lerch, 2022). |
| `functions_quantile_rd` | Functions for quantile reliability diagrams (based on code from Gneiting et al., 2023). |
| `data_preprocessing` | Preprocessing of the data (analogous to Yang et al., 2022). |
| `pp_AnEn` | Implementation of the AnEn-method (taken from Yang et al., 2022). |
| `pp_idr` | Postprocessing via IDR (based on code from Schulz and Lerch, 2022). |
| `pp_drn` | Postprocessing via DRN (based on code from Schulz and Lerch, 2022). |
| `pp_bqn` | Postprocessing via BQN (based on code from Schulz and Lerch, 2022). |
| `evaluation_scores` | Summary of the evaluation measures for the postprocessing methods. |
| `evaluation_quantile_data` | Calculations required for quantile reliability diagrams (based on code from Gneiting et al., 2023). |
| `figures_paper` | Generation of figures from the paper. |


## Licenses

The original data from Yang et al. (2022) and the additional data in the `/data/`-directory, as well as the figures in the `/figures/`-directory are distributed under a 'CC BY 4.0' license. The code in the `/code/`-directory is distributed under a 'MIT' license. Please refer to the 'LICENSE_data/figures/code' files in the respective folders for the corresponding license texts.


## Instructions & computational requirements

We ran the analysis in R (version 3.6.1) using the following packages:

  - IDR: `isodistrreg` (0.1.0.9000)
  - CORP reliability diagrams: `reliabilitydiag` (0.2.0)
  - (Proper) Scoring Rules: `scoringRules` (1.0.0)
  - `dplyr` (1.0.10)
  - `ggplot2` (3.3.6)
  - `isotone` (1.1-0)
  - `lubridate` (1.8.0)
  - `RANN` (2.6.1)
  - `tidyverse` (1.3.2)
  - `xtable` (1.8-4)
  
DRN and BQN are run with the R-Interface to Keras (v2.4.3; R-package version 2.3.0.0) and the TensorFlow backend (v2.3.0; R-package version 2.2.0), which is based on Python (v3.6). The implementation of DRN makes use of the tf-probability extension (0.11.1).
  
## Comments on the implementation of the postprocessing methods

# DRN

For DRN, we do not forecast GHI directly but instead use the bias of the ECMWF HRES forecast (forecast minus observation) as a target variable, which we model with a normal distribution. The GHI forecast is then obtained by subtracting the location parameter of the forecast distribution from the ECMWF HRES forecast. However, this does not ensure a positive forecast. Therefore, the forecast distribution is left-truncated at zero. The loss function used is the CRPS, the configuration is described below.

We also provide an alternative implementation of DRN that directly models GHI based on a truncated normal or logistic distribution (not used in the paper).

# BQN

As for DRN, we do not forecast GHI directly but instead the bias of the ECMWF HRES forecast. Here, we then substract the quantile function from the ECMWF HRES forecast, which also includes a change of orientation of the quantile forecast. Afterwards, we left-censor the forecasts at zero to ensure non-negativity. The loss function used to estimate the network parameters is the mean pinball loss over 99 equidistant quantiles on the unit interval, i.e. at steps of 1% at the levels 1%,..., 99%.

The PIT values of the BQN forecasts cannot be calculated directly, we thus rely on a generalization of PIT values, the uPIT (as in Schulz and Lerch, 2022). The evaluation of the BQN forecasts is based on a set of 99 equidistant quantiles for those measures that cannot be calculated exactly (CRPS, Brier score, CORP reliability diagrams).

We also provide an alternative implementation of BQN that directly models GHI (not used in the paper).

# IDR

Analogous to AnEN, IDR is applied separately to each station.

As mentioned in the paper, we here compare two variants of IDR. CSD-IDR was used in the paper and is based on the clear-sky model REST2, while GHI-IDR is based on the ECMWF HRES forecast. 

| CRPS | BON | DRA | FPK | GWN | PSU | SXF | TBL |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |  
| CSD-IDR | 82.9 | 37.9 | 62.4 | 85.6 | 84.6 | 74.8 | 74.2 | 
| GHI-IDR | 58.4 | 32.8 | 50.6 | 62.7 | 62.7 | 57.8 | 61.2 | 
| **Brier Score (250 W/m^2)** | **BON** | **DRA** | **FPK** | **GWN** | **PSU** | **SXF** | **TBL** |
| CSD-IDR | 0.141 | 0.040 | 0.107 | 0.133 | 0.142 | 0.124 | 0.104 | 
| GHI-IDR | 0.097 | 0.033 | 0.084 | 0.091 | 0.092 | 0.090 | 0.083 | 
| **PI Length** | **BON** | **DRA** | **FPK** | **GWN** | **PSU** | **SXF** | **TBL** |
| CSD-IDR | 387.8 | 233.0 | 324.6 | 414.0 | 401.6 | 386.6 | 398.4 | 
| GHI-IDR | 290.7 | 184.1 | 264.9 | 308.2 | 301.2 | 295.8 | 330.6 | 
| **Pinball Loss (8.33%)** | **BON** | **DRA** | **FPK** | **GWN** | **PSU** | **SXF** | **TBL** |
| CSD-IDR | 23.9 | 21.0 | 22.1 | 24.7 | 22.1 | 24.1 | 26.1 | 
| GHI-IDR | 18.7 | 15.0 | 16.8 | 19.8 | 17.1 | 19.9 | 21.3 | 
| **Pinball Loss (91.67%)** | **BON** | **DRA** | **FPK** | **GWN** | **PSU** | **SXF** | **TBL** |
| CSD-IDR | 13.5 | 6.1 | 10.0 | 13.9 | 15.6 | 11.8 | 11.9 | 
| GHI-IDR | 17.7 | 10.6 | 15.8 | 19.7 | 21.4 | 17.1 | 16.1 | 

## Network hyperparameters

The following configurations have been used for the network methods. Both network methods use station embedding for a global estimation (in contrast to AnEn and IDR that are estimated locally) and have been implemented using the Keras and Tensorflow-Interfaces to R. For further details, we refer to Schulz and Lerch (2022).

| Hyperparmaeter | DRN | BQN |
| --- | :---: | :---: |  
| Learning rate | $$5 \cdot 10^{-4}$$ | $$5 \cdot 10^{-4}$$ |
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


## References

- Gneiting, T., Wolffram, D., Resin, J., Kraus, K., Bracher, J., Dimitriadis, T., Hagenmeyer, V., Jordan, A.I., Lerch, S., Phipps, K., Schienle, M., 2023. Model diagnostics and forecast evaluation for quantiles. Annual Review of Statistics and Its Application 10, TBD–TBD. Replication material available at https://github.com/dwolffram/replication-ARSIA2023.
- Schulz, B., Lerch, S., 2022. Machine learning methods for postprocessing ensemble forecasts of wind gusts: A systematic comparison. Monthly Weather Review 150, 235–257. https://doi.org/10.1175/MWR-D-21-0150.1.
- Yang, D., Wang, W., Hong, T., 2022. A historical weather forecast dataset from the European Centre for Medium-Range Weather Forecasts (ECMWF) for energy forecasting. Solar Energy 232, 263–274. https://doi.org/10.1016/j.solener.2021.12.011.

