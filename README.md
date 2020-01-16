This repository contains code from the paper "Parametric g-formula implementations for causal survival analyses", co-authored by Young, Robins, Hernan and Wen. Folders “static_treatment” and “dynamic_treatment” and “random_treatment” contain codes to generate data for the static, dynamic and random treatment interventions  described in the main paper. In each of these folders, the files to reproduce the results found in the main paper (and in the supplementary materials) include

1. datagen.R: contains the function to generate data sets
2. deterministic_datagen_wide_true.R: code to produce the true parameter estimates 
3. HB_ICE.R: codes to produce the parameter estimates and standard errors from the pooled hazard-extended ICE estimator
4. HB_ICE_strat.R: codes to produce the parameter estimates and standard errors from the stratified hazard-extended ICE estimator 
5. sim_analysis.R: code to analyze the simulated data sets
