#!/bin/bash

#############################################
## Bayesian semiparametric Item Response Theory models using NIMBLE 
## Sally Paganin
## January 2021
#############################################
## All paths are relative to the paper_reproducible_code folder within the github repo

## first create a folder for saving all of the models and figures

mkdir -p figures

mkdir -p output
mkdir -p output/prior_samples                 
mkdir -p output/posterior_samples              
mkdir -p output/posterior_samples_elaborated  
mkdir -p output/mcmc_time

#############################################
## 1) RUN ALL MODELS 
#############################################
## UNIMODAL SIMULATION 
#############################################

## Parametric models

FILES=(models/parametric/*.R)

Rscript 1_runModels.R  \
--model=${FILES[$SLURM_ARRAY_TASK_ID]} \
--data=data/simulation_unimodal.rds \
--niter=50000 \
--nburnin=5000 \
--nthin=10 \
--mode=default

Rscript 1_runModels.R  \
--model=${FILES[$SLURM_ARRAY_TASK_ID]} \
--data=data/simulation_unimodal.rds \
--niter=50000 \
--nburnin=5000 \
--nthin=10 \
--mode=centered


## Semiparametric models

FILES=(models/bnp/*.R)

Rscript 1_runModels.R  \
--model=${FILES[$SLURM_ARRAY_TASK_ID]} \
--data=data/simulation_unimodal.rds \
--niter=50000 \
--nburnin=5000 \
--nthin=10 \
--mode=default

Rscript 1_runModels.R  \
--model=${FILES[$SLURM_ARRAY_TASK_ID]} \
--data=data/simulation_unimodal.rds \
--niter=50000 \
--nburnin=5000 \
--nthin=10 \
--mode=centered


#############################################
## BIMODAL SIMULATION 
#############################################
## Parametric models

FILES=(models/parametric/*.R)

Rscript 1_runModels.R  \
--model=${FILES[$SLURM_ARRAY_TASK_ID]} \
--data=data/simulation_bimodal.rds \
--niter=50000 \
--nburnin=5000 \
--nthin=10 \
--mode=default

Rscript 1_runModels.R  \
--model=${FILES[$SLURM_ARRAY_TASK_ID]} \
--data=data/simulation_bimodal.rds \
--niter=50000 \
--nburnin=5000 \
--nthin=10 \
--mode=centered


## Semiparametric models

FILES=(models/bnp/*.R)

Rscript 1_runModels.R  \
--model=${FILES[$SLURM_ARRAY_TASK_ID]} \
--data=data/simulation_bimodal.rds \
--niter=50000 \
--nburnin=5000 \
--nthin=10 \
--mode=default

Rscript 1_runModels.R  \
--model=${FILES[$SLURM_ARRAY_TASK_ID]} \
--data=data/simulation_bimodal.rds \
--niter=50000 \
--nburnin=5000 \
--nthin=10 \
--mode=centered

#############################################
## Stan models
#############################################

Rscript 2_runStanModel.R  \
--data=data/simulation_unimodal.rds \
--niter=4000 \
--nburnin=4000 

Rscript 2_runStanModel.R  \
--data=data/simulation_bimodal.rds \
--niter=4000 \
--nburnin=4000 

#############################################
## 2) extract results
#############################################
