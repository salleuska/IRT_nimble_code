##-----------------------------------------#
## Bayesian semiparametric Item Response Theory models using NIMBLE 
## Sally Paganin
## November 2020
##-----------------------------------------#
## This files contain the theme and set colors
## for all the plots in the paper
##-----------------------------------------#
library(ggplot2)

## Global theme
theme_set(theme_bw() +
          theme(axis.text = element_text(size= 11),
                axis.title=element_text(size=11),
                strip.text = element_text(size=11),
                strip.background =element_rect(fill="white"),
                plot.title = element_text(size=12)))

##################################################
##################################################
## Set colors for parametrization
SI_colors  <- c("#419e45",  "#a0b530", "#D2B21F", "#c26a00")
IRT_colors <- c("#0099ff", "#2dbede", "#935eaf")
stan_color <- "#B3385E"


## Set colors for modelling assumption
bnpColor  <- "#E69F00"
paraColor <- "#56B4E9"


labelData <- data.frame(R_label = c("IRT_stan", 
									"IRT_constrainedItem", 
									 "IRT_constrainedAbilities", 
									 "IRT_unconstrained", 
									 "SI_constrainedItem", "SI_constrainedAbilities", 
		 							 "SI_unconstrained",
		 							 "SI_unconstrained_centered"),
		 				 plot_label = c("IRT HMC (stan)", 
		 				 				"IRT constrained item",
		 				 				"IRT constrained abilities",
										"IRT unconstrained", 
										"SI constrained item", "SI constrained abilities", 
										"SI unconstrained",
										"SI unconstrained centered"),
		 				 colors = c(stan_color, IRT_colors, SI_colors))
labelData$colors <- as.character(labelData$colors)
## add colors to dataset