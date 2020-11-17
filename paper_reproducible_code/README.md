------------
## Code to reproduce analysis in  "Bayesian semiparametric Item Response Theory models using NIMBLE" 
------------

paper_reproducible_code

`
├── data                              # simulated and real world data 
├── figures                           # figures
├── models                            # nimble code for all IRT models
├── output                            # output from prior/posterior simulations
├── R_functions                       # R functions
├── README.md
├── Sec5_DPpriorNumberOfClusters.R
└── Sec5_priorMatching.R
`


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
 --nthin=  	      thinning interval for random effects   
 --mode=	        sampler types (default, centered, default_centered)  
------------------------------------------------------------
EXAMPLES
------------------------------------------------------------
Rscript runModels.R --model=models/parametric/2PL_unconstrained_gamma.R --dirResults=results/data_pisa/parametric/res/ --data=data/data_pisa.rds --niter=20000 --nburnin=0 --mode=default_constrained

Rscript runModels.R --model=models/parametric/2PL_unconstrained_gamma.R --dirResults=results/data_health/parametric/res/ --data=data/data_health.rds --niter=20000 --nburnin=0 --mode=default_constrained -->