##-----------------------------------------#
## Bayesian semiparametric Item Response Theory models using NIMBLE 
## Sally Paganin
## November 2020
##-----------------------------------------#
source("R_functions/ggplot_settings.R")
##-----------------------------------------#
## Set dimensions
plot_width  <- 21 ## a4 paper
plot_height <- plot_width/5*2
unit   <- "cm" 

dir <- "output/posterior_samples/mcmc_time/"
paraFileName <- "parametric_efficiency.txt"
bnpFileName <- "bnp_efficiency.txt"
##-----------------------------------------#
## Plot efficiency results for simulation
##-----------------------------------------#

unimodal <- read.table(paste0("output/mcmc_time/simulation_unimodal/", paraFileName), header = T)
bimodal  <- read.table(paste0("output/mcmc_time/simulation_bimodal/", paraFileName), header = T)

unimodal$ESS_second <- unimodal$ess_coda/unimodal$runningTime
bimodal$ESS_second  <- bimodal$ess_coda/bimodal$runningTime

bimodal$simulation  <- "Bimodal simulation"
unimodal$simulation <- "Unimodal simulation"

unimodal$labels <- gsub("parametric_", "", unimodal$fileName)
bimodal$labels  <- gsub("parametric_", "", bimodal$fileName)

## match R labels to plot labels
unimodal$labels <- labelData[match(unimodal$labels, labelData$R_label), ]$plot_label
bimodal$labels  <- labelData[match(bimodal$labels, labelData$R_label), ]$plot_label




## data frame for plotting
dfParametricEff <- data.frame(rbind(unimodal[, c("labels", "ESS_second", "simulation")],
						                        bimodal[, c("labels", "ESS_second", "simulation")]))

colnames(dfParametricEff) <- c("Strategy", "ESS", "Simulation")

dfParametricEff$Strategy  <- factor(dfParametricEff$Strategy, levels = labelData$plot_label)
dfParametricEff$Simulation <- factor(dfParametricEff$Simulation, levels = c("Unimodal simulation", "Bimodal simulation") )

##-----------------------------------------#
## Plot Figure 2a
##-----------------------------------------#
ylab <- paste0("min ESS/second (total time)")
p <-  ggplot(dfParametricEff,  aes_string(x = "Strategy", y= "ESS", fill = "Strategy")) +
      geom_bar(position= position_dodge(),stat='identity',colour = "black",
       width = 0.8) +
      facet_wrap(~ Simulation, ncol=2, scales='fixed') +
      ylab("min ESS/second (total time)") + xlab("") + 
      theme(legend.position = "none") +
      coord_flip() +
      scale_fill_manual(values = labelData$colors) +
      scale_x_discrete(limits = rev(levels(dfParametricEff$Strategy))) 

p


ggsave(filename = "figures/fig3a_simulation_efficiencies.png", plot = p,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')

##-----------------------------------------#
## Plot Figure 2b
## Comparison with sampling times
##-----------------------------------------#
bimodal$ESS_second2 <- bimodal$ess_coda/bimodal$samplingTime
unimodal$ESS_second2 <- unimodal$ess_coda/unimodal$samplingTime

dfParametricEffSampling <- dfParametricEff
dfParametricEffSampling$ESS <- c(unimodal$ESS_second2, bimodal$ESS_second2)


ylabel <- 'min ESS/second'
title <- paste0("Minimum effective sample size per second (sampling time)")
p <-  ggplot(dfParametricEffSampling,  aes_string(x = "Strategy", y= "ESS", fill = "Strategy")) +
      geom_bar(position= position_dodge(),stat='identity',colour = "black", width = 0.8) +
      facet_wrap(~ Simulation, ncol=2, scales='fixed') +
      ylab("min ESS/second (sampling time)") + xlab("") + 
      theme(legend.position = "none") +
      coord_flip() +
      scale_x_discrete(limits = rev(levels(dfParametricEffSampling$Strategy))) +
      scale_fill_manual(values = labelData$colors) 
p

ggsave(filename = "figures/fig3b_simulation_efficiencies_sampling.png", plot = p,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')

##-----------------------------------------#
## Efficiency fof bnp vs parametric - simulation
##-----------------------------------------#


bnpUnimodal <- read.table(paste0("output/mcmc_time/simulation_unimodal/", bnpFileName), header = T)
bnpBimodal  <- read.table(paste0("output/mcmc_time/simulation_bimodal/", bnpFileName), header = T)

bnpUnimodal$ESS_second <- bnpUnimodal$ess_coda/bnpUnimodal$runningTime
bnpBimodal$ESS_second <- bnpBimodal$ess_coda/bnpBimodal$runningTime

bnpBimodal$simulation <- "Bimodal simulation"
bnpUnimodal$simulation <- "Unimodal simulation"

bnpUnimodal$labels <- gsub("bnp_", "", bnpUnimodal$fileName)
bnpBimodal$labels <- gsub("bnp_", "", bnpBimodal$fileName)

## match R labels to plot labels
bnpUnimodal$labels <- droplevels(labelData[match(bnpUnimodal$labels, labelData$R_label), ]$plot_label)
bnpBimodal$labels  <- droplevels(labelData[match(bnpBimodal$labels, labelData$R_label), ]$plot_label)

## select parametric models with bnp equivalent and create data frame for plotting
dfParametricBnp <- data.frame(rbind(unimodal[which(unimodal$label %in% levels(bnpUnimodal$labels)), 
							c("labels", "ESS_second", "simulation")], 
						   bimodal[which(bimodal$label %in% levels(bnpBimodal$labels)), 
							c("labels", "ESS_second", "simulation")])) 
dfParametricBnp <- droplevels(dfParametricBnp)
dfParametricBnp$model <- "Parametric"

bnpUnimodal$model <- "Semiparametric"						  
bnpBimodal$model  <- "Semiparametric"						  

dfParametricBnp <- rbind(dfParametricBnp, 
                         bnpUnimodal[, c("labels", "ESS_second", "simulation", "model")],
 				                 bnpBimodal[, c("labels", "ESS_second", "simulation", "model")])


colnames(dfParametricBnp) <- c("Strategy", "ESS", "Simulation", "Model")
dfParametricBnp$Strategy  <- droplevels(factor(dfParametricBnp$Strategy, levels = labelData$plot_label))
## set level order
dfParametricBnp$Simulation <- factor(dfParametricBnp$Simulation, levels = c("Unimodal simulation", "Bimodal simulation"))

colorsParametricBnp <- labelData$colors[match(levels(dfParametricBnp$Strategy), labelData$plot_label)]
##-----------------------------------------#
## Plot Figure 3a
##-----------------------------------------#

ylabel <- 'Effective sample size per second'
title <- paste0("Minimum effective sample size per second (total time)")

p <-  ggplot(dfParametricBnp,  aes_string(x = "Strategy", y= "ESS", fill = "Strategy")) +
      geom_bar(position= position_dodge(),stat='identity',colour = "black",
       width = 0.8) +
      facet_wrap(~ Model + Simulation,ncol=2, scales='fixed') +
      ylab("min ESS/second (total time)") + xlab("") + 
      theme(legend.position = "none") +
      coord_flip() +
      scale_x_discrete(limits = rev(levels(dfParametricBnp$Strategy))) +
      scale_fill_manual(values = colorsParametricBnp) 
p

ggsave(filename = "figures/fig4_simulation_efficiencies_bnp.png", plot = p,
        width = plot_width, height = plot_height/2*3 , dpi = 300, units = unit, device='png')

##-----------------------------------------#
## Real data - efficiency comparisons
##-----------------------------------------#
paraHealth <- read.table(paste0("output/mcmc_time/data_health/", paraFileName), header = T)
paraTimss  <- read.table(paste0("output/mcmc_time/data_timss/", paraFileName), header = T)

paraHealth$ESS_second <- paraHealth$ess_coda/paraHealth$runningTime
paraTimss$ESS_second  <- paraTimss$ess_coda/paraTimss$runningTime

paraTimss$simulation  <- "TIMSS data"
paraHealth$simulation <- "Health data"

paraHealth$labels <- gsub("parametric_", "", paraHealth$fileName)
paraTimss$labels  <- gsub("parametric_", "", paraTimss$fileName)

## match R labels to plot labels
paraHealth$labels <- labelData[match(paraHealth$labels, labelData$R_label), ]$plot_label
paraTimss$labels  <- labelData[match(paraTimss$labels, labelData$R_label), ]$plot_label

paraHealth$model <- "Parametric"             
paraTimss$model  <- "Parametric"             

## data frame for plotting
dfParametricEff <- data.frame(rbind(paraHealth[, c("labels", "ESS_second", "simulation", "model")],
                                    paraTimss[, c("labels", "ESS_second", "simulation","model")]))

colnames(dfParametricEff) <- c("Strategy", "ESS", "Simulation", "Model")

dfParametricEff$Strategy  <- droplevels(factor(dfParametricEff$Strategy, levels = labelData$plot_label))
dfParametricEff <- droplevels(dfParametricEff[-grep("constrained item", dfParametricEff$Strategy), ])

colorsParametricData <- labelData$colors[match(levels(dfParametricEff$Strategy), labelData$plot_label)]
##-----------------------------------------#
## Plot Figure 5
##-----------------------------------------#
ylab <- paste0("min ESS/second (total time)")
p <-  ggplot(dfParametricEff,  aes_string(x = "Strategy", y= "ESS", fill = "Strategy")) +
      geom_bar(position= position_dodge(),stat='identity',colour = "black",
       width = 0.8) +
      facet_wrap(~ Model + Simulation, ncol=2, scales='free_x') +
      ylab("min ESS/second (total time)") + xlab("") + 
      theme(legend.position = "none") +
      coord_flip() +
      scale_fill_manual(values = colorsParametricData) +
      scale_x_discrete(limits = rev(levels(dfParametricEff$Strategy))) 

p

ggsave(filename = "figures/fig5_data_efficiencies.png", plot = p,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')


##-----------------------------------------#
## Compare bnp efficiencies - data
##-----------------------------------------#

bnpHealth <- read.table(paste0("output/mcmc_time/data_health/", bnpFileName), header = T)
bnpTimss  <- read.table(paste0("output/mcmc_time/data_timss/", bnpFileName), header = T)

bnpHealth$ESS_second <- bnpHealth$ess_coda/bnpHealth$runningTime
bnpTimss$ESS_second <- bnpTimss$ess_coda/bnpTimss$runningTime

bnpTimss$simulation  <- "TIMSS data"
bnpHealth$simulation <- "Health data"

bnpHealth$labels <- gsub("bnp_", "", bnpHealth$fileName)
bnpTimss$labels <- gsub("bnp_", "", bnpTimss$fileName)

## match R labels to plot labels
bnpHealth$labels <- labelData[match(bnpHealth$labels, labelData$R_label), ]$plot_label
bnpTimss$labels  <- labelData[match(bnpTimss$labels, labelData$R_label), ]$plot_label

bnpHealth$model <- "Semiparametric"             
bnpTimss$model  <- "Semiparametric"             
## data frame for plotting
dfBnpEff <- data.frame(rbind(bnpHealth[, c("labels", "ESS_second", "simulation", "model")],
                                    bnpTimss[, c("labels", "ESS_second", "simulation","model")]))

colnames(dfBnpEff) <- c("Strategy", "ESS", "Simulation", "Model")

dfBnpEff$Strategy  <- droplevels(factor(dfBnpEff$Strategy, levels = labelData$plot_label))
dfBnpEff <- droplevels(dfBnpEff[-grep("constrained item", dfBnpEff$Strategy), ])

colorsBnp <- labelData$colors[match(levels(dfBnpEff$Strategy), labelData$plot_label)]
##-----------------------------------------#
## Plot Figure 5b
##-----------------------------------------#

ylabel <- 'Effective sample size per second'
title <- paste0("Minimum effective sample size per second (total time)")

p <-  ggplot(dfBnpEff,  aes_string(x = "Strategy", y= "ESS", fill = "Strategy")) +
      geom_bar(position= position_dodge(),stat='identity',colour = "black",
       width = 0.8) +
      facet_wrap(~ Model + Simulation,ncol=2, scales='free_x') +
      ylab("min ESS/second (total time)") + xlab("") + 
      theme(legend.position = "none") +
      coord_flip() +
      scale_x_discrete(limits = rev(levels(dfBnpEff$Strategy))) +
      scale_fill_manual(values = colorsBnp) 
p

ggsave(filename = "figures/fig5b_data_efficiencies_bnp.png", plot = p,
        width = plot_width, height = plot_height/6*4 , dpi = 300, units = unit, device='png')

