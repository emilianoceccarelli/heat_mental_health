# Code and anonymized data for article https://doi.org/10.1111/risa.70267

This repository contains the R code and anonymized data related to the article:

**Heat-Related Mental Health Hospitalizations in Italy: A Global Sensitivity Analysis Approach to Evaluate Generalized Additive Model Assumptions**

## Data

The file `db_anon_MJJAS.rds` contains a randomized version of the original dataset.  
Due to privacy restrictions, the original data cannot be publicly shared. To make the code executable, the outcome variable has been randomly permuted while preserving the overall structure of the dataset.

As a consequence, this dataset allows users to reproduce the computational workflow, but not the numerical results reported in the article.

## GAM model

The file `01_GAM.R` contains the R code used to estimate the GAM model described in the article.

The script fits a negative binomial GAM including the distributed lag non-linear model crossbasis for daily maximum temperature, temporal covariates, population offset, and random effects for geographical macro-area.

## Global sensitivity analysis

The file `02_GSA.R` contains the code used to perform the global sensitivity analysis (GSA).

The GSA evaluates the sensitivity of the model outputs to alternative modelling choices, including:

- the number and placement of knots for the temperature-response function;
- the number and placement of knots for the lag-response function;
- the inclusion or exclusion of September from the analysis period.

The GSA script is structured to be run on a high-performance computing (HPC) system. In particular, the code uses parallel computing through the `doParallel` and `foreach` packages to distribute model fitting across multiple cores.

## Important note

Because the dataset is randomized for privacy reasons, the estimated effects produced by the code are going to differ from those shown in the article.
