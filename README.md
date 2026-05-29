# Code and anonymized data for article https://doi.org/10.1111/risa.70267

This repository contains the R code and anonymized data related to the article:

**Heat-Related Mental Health Hospitalizations in Italy:A Global Sensitivity Analysis Approach to Evaluate Generalized Additive Model Assumptions**

## Data

The file `db_anon_MJJAS.rds` contains a randomized version of the original dataset.  
Due to privacy restrictions, the original data cannot be publicly shared. To make the code executable, the outcome variable has been randomly permuted while preserving the overall structure of the dataset.

As a consequence, this dataset allows users to reproduce the computational workflow, but not the numerical results reported in the article.

## Model code

The file `01_GAM.R` contains the R code used to estimate the GAM model described in the article.

## Important note

Because the dataset is randomized for privacy reasons, the estimated effects, relative risks, confidence intervals, tables, and figures produced by the code may differ from those shown in the article.
