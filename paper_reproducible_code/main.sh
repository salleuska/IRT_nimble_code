#!/bin/bash

#############################################
## Bayesian semiparametric Item Response Theory models using NIMBLE 
## Sally Paganin
## July 2022
#############################################

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

## The script simulate from the prior predictive for all models.
## Results can be used to select model hyperpriors 
## as suggested in Section 4

Rscript sec4_priorMatching.R
Rscript sec4_makePlots.R

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
--nthin=1 \
--mode=default

Rscript 1_runModels.R  \
--model=${FILES[$SLURM_ARRAY_TASK_ID]} \
--data=data/simulation_unimodal.rds \
--niter=50000 \
--nburnin=5000 \
--nthin=1 \
--mode=centered


## Semiparametric models

FILES=(models/bnp/*.R)

Rscript 1_runModels.R  \
--model=${FILES[$SLURM_ARRAY_TASK_ID]} \
--data=data/simulation_unimodal.rds \
--niter=50000 \
--nburnin=5000 \
--nthin=1 \
--mode=default

Rscript 1_runModels.R  \
--model=${FILES[$SLURM_ARRAY_TASK_ID]} \
--data=data/simulation_unimodal.rds \
--niter=50000 \
--nburnin=5000 \
--nthin=1 \
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
--nthin=1 \
--mode=default

Rscript 1_runModels.R  \
--model=${FILES[$SLURM_ARRAY_TASK_ID]} \
--data=data/simulation_bimodal.rds \
--niter=50000 \
--nburnin=5000 \
--nthin=1 \
--mode=centered

## Semiparametric models

FILES=(models/bnp/*.R)

Rscript 1_runModels.R  \
--model=${FILES[$SLURM_ARRAY_TASK_ID]} \
--data=data/simulation_bimodal.rds \
--niter=50000 \
--nburnin=5000 \
--nthin=1 \
--mode=default

Rscript 1_runModels.R  \
--model=${FILES[$SLURM_ARRAY_TASK_ID]} \
--data=data/simulation_bimodal.rds \
--niter=50000 \
--nburnin=5000 \
--nthin=1 \
--mode=centered

#############################################
## TIMSS DATA 
#############################################
## Parametric models

FILES=(models/parametric_long/*.R)

Rscript 1_runModels.R  \
--model=${FILES[$SLURM_ARRAY_TASK_ID]} \
--data=data/data_timss.rds \
--niter=50000 \
--nburnin=5000 \
--nthin=1 \
--mode=default

Rscript 1_runModels.R  \
--model=${FILES[$SLURM_ARRAY_TASK_ID]} \
--data=data/data_timss.rds \
--niter=50000 \
--nburnin=5000 \
--nthin=1 \
--mode=centered

## Semiparametric models

FILES=(models/bnp_long/*.R)

Rscript 1_runModels.R  \
--model=${FILES[$SLURM_ARRAY_TASK_ID]} \
--data=data/data_timss.rds \
--niter=50000 \
--nburnin=5000 \
--nthin=1 \
--mode=default

Rscript 1_runModels.R  \
--model=${FILES[$SLURM_ARRAY_TASK_ID]} \
--data=data/data_timss.rds \
--niter=50000 \
--nburnin=5000 \
--nthin=1 \
--mode=centered

#############################################
## HEALTH DATA 
#############################################
## Parametric models

FILES=(models/parametric/*.R)

Rscript 1_runModels.R  \
--model=${FILES[$SLURM_ARRAY_TASK_ID]} \
--data=data/data_health.rds \
--niter=50000 \
--nburnin=5000 \
--nthin=1 \
--mode=default

Rscript 1_runModels.R  \
--model=${FILES[$SLURM_ARRAY_TASK_ID]} \
--data=data/data_health.rds \
--niter=50000 \
--nburnin=5000 \
--nthin=1 \
--mode=centered

## Semiparametric models

FILES=(models/bnp/*.R)

Rscript 1_runModels.R  \
--model=${FILES[$SLURM_ARRAY_TASK_ID]} \
--data=data/data_health.rds \
--niter=50000 \
--nburnin=5000 \
--nthin=1 \
--mode=default

Rscript 1_runModels.R  \
--model=${FILES[$SLURM_ARRAY_TASK_ID]} \
--data=data/data_health.rds \
--niter=50000 \
--nburnin=5000 \
--nthin=1 \
--mode=centered

#############################################
## Stan models
#############################################

Rscript 1_runStanModel.R  \
--data=data/simulation_unimodal.rds \
--nsamples=10000 \
--nwarmup=5000 

Rscript 1_runStanModel.R  \
--data=data/simulation_bimodal.rds \
--nsamples=10000 \
--nwarmup=5000 

Rscript 1_runStanModel.R  \
--data=data/data_health.rds \
--nsamples=10000 \
--nwarmup=5000 

Rscript 1_runStanModel.R  \
--data=data/data_timss.rds \
--nsamples=10000 \
--nwarmup=5000 

#############################################
## 2) extract MCMC samples and postprocess them
#############################################

## Extract results for parametric models 

for filename in output/posterior_samples/simulation_unimodal/parametric/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename
done

for filename in output/posterior_samples/simulation_bimodal/parametric/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename
done

for filename in output/posterior_samples/data_timss/parametric/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename
done

for filename in output/posterior_samples/data_health/parametric/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename
done

############################################################
## Extract results for semiparametric models 

for filename in output/posterior_samples/simulation_unimodal/bnp/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename
done

for filename in output/posterior_samples/simulation_bimodal/bnp/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename
done

for filename in output/posterior_samples/data_timss/bnp/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename
done


for filename in output/posterior_samples/data_health/bnp/*.rds; do
	echo $filename
	Rscript 2_extractResults.R --resFileName=$filename
done

#############################################
## 3) Simulate from DP measure 
#############################################
## Using more efficient MCMC strategies for each 
## simulated scenario

Rscript 3_simulateFromDPmeasure.R \
--dataName="simulation_unimodal" \
--modelName="bnp_IRT_unconstrained"

Rscript 3_simulateFromDPmeasure.R \
--dataName="simulation_bimodal" \
--modelName="bnp_IRT_unconstrained"

Rscript 3_simulateFromDPmeasure.R \
--dataName="simulation_multimodal" \
--modelName="bnp_IRT_unconstrained"

Rscript 3_simulateFromDPmeasure.R \
--dataName="data_health" \
--modelName="bnp_IRT_unconstrained"

Rscript 3_simulateFromDPmeasure.R \
--dataName="data_timss" \
--modelName="bnp_IRT_unconstrained"

## 3PL model for TIMSS data
Rscript 3_simulateFromDPmeasure.R \
--dataName="data_timss" \
--modelName="bnp3PL_IRT_unconstrained"

#############################################
## 4) Make plots
############################################

## Extract some quantities to reproduce the plots
Rscript 4_computeQuantitesForFigures.R

## Plots Section 6 
Rscript sec6_makePlots.R

## Plots Section 7
Rscript sec7_dataResults.R
Rscript sec7_simulationsResults.R


