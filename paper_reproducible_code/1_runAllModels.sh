#!/bin/bash

#############################################
## Bayesian semiparametric Item Response Theory models using NIMBLE 
## Sally Paganin
## January 2021

mkdir -p figures
mkdir -p figures/dataForFigures

mkdir -p output
mkdir -p output/prior_samples                 
mkdir -p output/posterior_samples              
mkdir -p output/posterior_samples_elaborated  
mkdir -p output/mcmc_time


#############################################
## 0) Simulate from prior predictive 
#############################################

## The script simulate from prior predictive for all models 
## results can be used to adjust model hyperpriors 
## as suggested in Sec 5

Rscript sec5_priorMatching.R

#############################################
## 1) RUN ALL MODELS 
#############################################
## Run all models for simulated data 
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
## 2) extract MCMC samples and postprocess them
#############################################

## Extract results for parametric models 

for filename in output/posterior_samples/simulation_unimodal/parametric/*.rds; do
	echo $filename
	Rscript 3_extractResults.R --resFileName=$filename
done

for filename in output/posterior_samples/simulation_bimodal/parametric/*.rds; do
	echo $filename
	Rscript 3_extractResults.R --resFileName=$filename
done
############################################################
## Extract results for semiparametric models 
 
for filename in output/posterior_samples/simulation_unimodal/bnp/*.rds; do
	echo $filename
	Rscript 3_extractResults.R --resFileName=$filename
done

for filename in output/posterior_samples/simulation_bimodal/bnp/*.rds; do
	echo $filename
	Rscript 3_extractResults.R --resFileName=$filename
done

#############################################
## 3) Simulate from DP measure 
#############################################
## Using more efficient MCMC strategies for each 
## simulated scenario

Rscript 4_simulateFromDPmeasure.R \
--dataName="simulation_unimodal" \
--modelName="bnp_SI_unconstrained"

Rscript simulateFromDPmeasure.R \
--dataName="simulation_bimodal" \
--modelName="bnp_IRT_unconstrained"

# Rscript 4_simulateFromDPmeasure.R \
# --dataName="data_health" \
# --modelName="bnp_IRT_unconstrained"

# Rscript 4_simulateFromDPmeasure.R \
# --dataName="data_timss" \
# --modelName="bnp_IRT_unconstrained"

