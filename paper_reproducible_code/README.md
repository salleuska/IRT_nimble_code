------------
## Code to reproduce analysis in  "Computational strategies and estimation performance with Bayesian semiparametric Item Response Theory model" 
------------

Make sure to have the following `R` packages installled.
```r
install.packages("nimble")
install.packages("rstan")

## packages for output processing and visualizion
install.packages("reshape2")
install.packages("ggplot2")
install.packages("bayestestR")
install.packages("cowplot")

install.packages("sn") ## skew normal for simulations

```

### Folder organization

`paper_reproducible_code` contains all the code to reproduce analysis in the paper. It is organized as

```bash
├── R_functions                       # R functions
├── data                              # simulated data with R scripts for simulation
├── models                            # nimble code for all IRT 2PL and 3PL models + stan code for 2PL model
```
Two other folders will be created when running the `main.sh` scripts. An `output` folder containing all outputs generated using the models 

```bash
├── output                            # output from prior/posterior simulations
│   ├── prior_samples                 # samples from prior predictive simulation [sec5_priorMatching.R]
│   ├── posterior_samples             # raw MCMC samples [1_runNimbleModels.R] 
│   ├── posterior_samples_elaborated  # postprocessed MCMC samples [3_extractResults.R]
│   └── mcmc_time                     # time and ESS for different MCMC to to compute efficiencies [3_extractResults.R]
``` 
and a `figures` folder containing all the plots in the paper and data to generate them in the `dataForFigures` folder.


### `R` scripts

Scripts `1` to `3` run all models, extract and postprocess posterior samples, and use them to obtain quantities for inference. 
Scripts startiting with `sec_` reproduce prior analysis and plots relative to the section number. 

```bash
├── 1_runNimbleModels.R               # run models coded in NIMBLE
├── 1_runStanModel.R                  # run models coded in Stan
├── 2_extractResults.R                # postprocess raw posterior samples 
├── 3_simulateFromDPmeasure.R         # simulate from DP posterior 
├── computeQuantitesForFigures.R      # use postprocessed samples to compute quantities for figures
├── sec5_DPpriorNumberOfClusters.R    # calculate prior expectation and variance for n. of clusters of DP prior
├── sec5_makePlots.R                  # make plots for sec. 5
├── sec5_priorMatching.R              # simulate from model priors
├── sec6_makePlots.R                  # make plots for sec. 6
├── sec7_makePlots.R                  # make plots for sec. 7
```

### Usage

The bash script `main.sh` runs all NIMBLE/Stan models on the simulated data, simulates from the prior, extracts and postprocesses posterior samples, populating the `output` folder. Running all the models can be memory/time consuming. 

Scripts `1` to `3` are parametrized and can be run separately. Details are in the comments at the beginning of the file. For example:

`Rscript 1_runNimbleModels.R --model --dirResults --data --niter --nburnin --mode`

```bash
 --model=           # path to the model code to to run  
 --dirResults=      # directory to results  
 --data=            # directory to data   
 --niter=           # number of iterations  
 --nburnin=         # number of burnin iteration  
 --nthin=           # thinning interval for random effects  (will be thin2 in nimble) 
 --mode=            # sampler types (default, centered, default_centered)  
```

