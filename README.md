## Playeraging
This will be repository for data and code from Schuckers, Lopez and Macdonald 
_Estimation of Player Aging Curves Using Regression and Imputation_
Annals of Operations Research 325(1): 681-699 (2023)

Here is the link to the arXiv version of this paper: https://arxiv.org/abs/2110.14017

### Files

*NewAgingModelsTesting.R* is a file of functions for generation of simulated player aging data and for fitting the models in this paper.

*Fit_Models_NHL_Data.R* fits the models in our paper to data for NHL players found in the file "forwards_zppg_missing_long_2022.csv" and outputs the results in 
the file "Forwards_model_out_071022.csv".

*forwards_zppg_missing_long_2022.csv* is a file with data on players from the National Hockey League in a long format tabular format with z-scores

*Forwards_model_out_071022.csv* is the results of fitting the models described in our paper.

*testsim_1_1000_0.8_4.R* is an example file that generates a specific set of simuated data, fits models to those data and reports the estimated mean aging curve from those models.  This code outputs to files, a file "outsingle<_>.csv" which is the data for all players for a single simulation and "outtest_<_>.csv" which is the 
Primary function here is _playeragingsim_.  This function also produces _outsingle*.csv_

*NHL_forwards_models_out_071022.csv* is a example file that has an instance of the estimates of the 'aging curve' for each of the model discussed in the paper. 

*forwards_nhl_aging2020_long.csv* is a data file where each row is an observed season for a player.

*out_sim_all_073121.csv* is a realization of 1000 simulated datasets 

*outsingle062821_0_600_0.4_4q.csv* is a single realization of the simulated data 







*
