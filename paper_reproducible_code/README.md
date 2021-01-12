------------
## Code to reproduce analysis in  "Bayesian semiparametric Item Response Theory models using NIMBLE" 
------------

### Folder organization

`paper_reproducible_code` contains all the code to reproduc analysis in the paper. It is organized as

```bash
├── R_functions                       # R functions
├── data                              # simulated data with R scripts for simulation
├── models                            # nimble code for all IRT models + stan code for 2PL model
```
Two other folders will be created when running the `main.sh` scripts. An `output` folder containing all outputs generated using the models 

```bash
├── output                            # output from prior/posterior simulations
│   ├── prior_samples                 # samples from prior predictive simulation [sec5_priorMatching.R]
│   ├── posterior_samples             # raw MCMC samples [1_runNimbleModels.R] 
│   ├── posterior_samples_elaborated  # postprocessed MCMC samples [3_extractResults.R]
│   └── mcmc_time					  # time and ESS for different MCMC to to compute efficiencies [3_extractResults.R]
``` 
and a `figures` folder containing all the plots in the paper and data to generate them in the `dataForFigures` folder.





### `R` scripts

Scripts `1` to `5` run all models, extract and postprocess posterior samples, and use them to obtain quantities for inference. 
Scripts startiting with `sec_` reproduce prior analysis and plots relative to the section number. 

```bash
├── 1_runNimbleModels.R               # run models coded in NIMBLE
├── 2_runStanModel.R                  # run models coded in Stan
├── 3_extractResults.R                # postprocess raw posterior samples 
├── 4_simulateFromDPmeasure.R         # simulate from DP posterior 
├── 5_computeQuantitesForFigures.R    # use postprocessed samples to compute quantities for figures
├── sec5_DPpriorNumberOfClusters.R    # calculate prior expectation and variance for n. of clusters of DP prior
├── sec5_makePlots.R                  # make plots for sec. 5
├── sec5_priorMatching.R              # simulate from model priors
├── sec6_makePlots.R                  # make plots for sec. 6
├── sec7_makePlots.R                  # make plots for sec. 7
```


### Usage


<!-- 
`R`  1_runModels.R --model --dirResults --data --niter --nburnin --mode

```bash
 --model=         # path to the model code to to run  
 --dirResults=    # directory to results  
 --data=          # directory to data   
 --niter=  	      # number of iterations  
 --nburnin=       # number of burnin iteration  
 --nthin=  	      # thinning interval for random effects  (will be thin2 in nimble) 
 --mode=	      # sampler types (default, centered, default_centered)  
```

 -->
<!-- 

data -- real + synthetic data + scripts for simulations

models -- nimble model code  
  parametric       -- parametric 2PL models  
  parametric_long  -- parametric 2PL models, data in long format  
  bnp              -- semi-parametric 2PL models  
  bnp_long         -- semi-parametric 2PL models, data in long format  

Name format for models

[parametric|bnp]_constraintType_parametrization  



results -- res + posterior samples for each data/model type + markdown reports  
  unimodal_parametric  
  bimodal_parametric  
  health_data  
  timss_data  

util -- R scripts with utilities   
  customSamplers.R -- implemented custom samplers (centered)  

prior_simulations -- simulated probabilities form the prior distribution
  parametric
  bnp
------------

------------
 USAGE
------------

Bash call to run models 1/2/3/4_*.sh

Rscript 1_runModels.R --model --dirResults --data --niter --nburnin --mode

arguments  
 --model=         path to the model code to to run  
 --dirResults=    directory to results  
 --data=          directory to data   
 --niter=  	      number of iterations  
 --nburnin=       number of burnin iteration  
 --nthin=  	      thinning interval for random effects  (will be thin2 in nimble) 
 --mode=	        sampler types (default, centered, default_centered)  
------------------------------------------------------------
EXAMPLES
------------------------------------------------------------
Rscript runModels.R --model=models/parametric/2PL_unconstrained_gamma.R --dirResults=results/data_pisa/parametric/res/ --data=data/data_pisa.rds --niter=20000 --nburnin=0 --mode=default_constrained

Rscript runModels.R --model=models/parametric/2PL_unconstrained_gamma.R --dirResults=results/data_health/parametric/res/ --data=data/data_health.rds --niter=20000 --nburnin=0 --mode=default_constrained -->