##-----------------------------------------#
## Computational strategies and estimation performance with Bayesian semiparametric Item Response Theory model
## Sally Paganin
## last update: August 2022
## R version 4.2.0 (2022-04-22) -- "Vigorous Calisthenics"
## nimble version 0.12.2
##-----------------------------------------#
source("R_functions/ggplot_settings.R") 
##-----------------------------------------#
## Set dimensions
plot_width  <- 21 ## a4 paper
plot_height <- plot_width/3.5
unit   <- "cm" 

dir <- "output/posterior_samples/mcmc_time/"
paraFileName <- "parametric_efficiency.txt"
bnpFileName <- "bnp_efficiency.txt"
##-----------------------------------------#
## Plot efficiency results for simulation
##-----------------------------------------#
unimodal <- read.table(paste0("output/mcmc_time/simulation_unimodal/", paraFileName), header = T)
bimodal  <- read.table(paste0("output/mcmc_time/simulation_bimodal/", paraFileName), header = T)
multimodal  <- read.table(paste0("output/mcmc_time/simulation_multimodal/", paraFileName), header = T)

unimodal$ESS_second <- unimodal$multiEssItemsAbility/unimodal$runningTime
bimodal$ESS_second  <- bimodal$multiEssItemsAbility/bimodal$runningTime
multimodal$ESS_second  <- multimodal$multiEssItemsAbility/multimodal$runningTime

unimodal$simulation     <- "Unimodal simulation"
bimodal$simulation      <- "Bimodal simulation"
multimodal$simulation   <- "Multimodal simulation"

unimodal$labels <- gsub("parametric_", "", unimodal$fileName)
bimodal$labels  <- gsub("parametric_", "", bimodal$fileName)
multimodal$labels <- gsub("parametric_", "", multimodal$fileName)


## match R labels to plot labels
unimodal$labels <- labelData[match(unimodal$labels, labelData$R_label), ]$plot_label
bimodal$labels  <- labelData[match(bimodal$labels, labelData$R_label), ]$plot_label
multimodal$labels  <- labelData[match(multimodal$labels, labelData$R_label), ]$plot_label

## Remove extra strategy (SI constrained abilities centered) 
if(sum(is.na(unimodal$labels)) > 0) unimodal <- droplevels(unimodal[!is.na(unimodal$labels), ])
if(sum(is.na(bimodal$labels)) > 0) bimodal <- droplevels(bimodal[!is.na(bimodal$labels), ])
if(sum(is.na(multimodal$labels)) > 0) multimodal <- droplevels(multimodal[!is.na(multimodal$labels), ])

# levels(unimodal$labels)[grep("stan", levels(unimodal$labels))] <- paste0(levels(unimodal$labels)[grep("stan", levels(unimodal$labels))], "*")
# levels(bimodal$labels)[grep("stan", levels(bimodal$labels))] <- paste0(levels(bimodal$labels)[grep("stan", levels(bimodal$labels))], "*")
# levels(multimodal$labels)[grep("stan", levels(multimodal$labels))] <- paste0(levels(multimodal$labels)[grep("stan", levels(multimodal$labels))], "*")

## data frame for plotting
dfParametricEff <- data.frame(rbind(unimodal[, c("labels", "ESS_second", "simulation")],
					rbind(bimodal[, c("labels", "ESS_second", "simulation")],
                                  multimodal[, c("labels", "ESS_second", "simulation")]) ))

colnames(dfParametricEff) <- c("Strategy", "ESS", "Simulation")

dfParametricEff$Strategy  <- factor(dfParametricEff$Strategy, levels = labelData$plot_label)
dfParametricEff$Simulation <- factor(dfParametricEff$Simulation, levels = c("Unimodal simulation", "Bimodal simulation",  "Multimodal simulation") )

##-----------------------------------------#
## Plot Figure 2a
##-----------------------------------------#
ylab <- paste0("mESS/second (total time)")
p <-  ggplot(dfParametricEff,  aes_string(x = "Strategy", y= "ESS", fill = "Strategy")) +
      geom_bar(position= position_dodge(),stat='identity',colour = "black",
       width = 0.8) +
      facet_wrap(~ Simulation, ncol=3, scales='fixed') +
      ylab("mESS/second (total time)") + xlab("") + 
#      ylim(c(0,6)) + 
      theme(legend.position = "none") +
      coord_flip() +
      scale_fill_manual(values = labelData$colors) +
      scale_x_discrete(limits = rev(levels(dfParametricEff$Strategy))) 

p

ggsave(filename = "figures/fig3a_simulation_efficiencies.png", plot = p,
        width = plot_width, height = plot_height, scale = 1.2, dpi = 300, units = unit, device='png')

##-----------------------------------------#
## Plot Figure 2b
## Comparison with sampling times
##-----------------------------------------#
bimodal$ESS_second2 <- bimodal$multiEssItemsAbility/bimodal$samplingTime
unimodal$ESS_second2 <- unimodal$multiEssItemsAbility/unimodal$samplingTime
multimodal$ESS_second2 <- multimodal$multiEssItemsAbility/multimodal$samplingTime

dfParametricEffSampling <- dfParametricEff
dfParametricEffSampling$ESS <- c(unimodal$ESS_second2, bimodal$ESS_second2,multimodal$ESS_second2)


ylabel <- 'mESS/second'
title <- paste0("Minimum effective sample size per second (sampling time)")
p <-  ggplot(dfParametricEffSampling,  aes_string(x = "Strategy", y= "ESS", fill = "Strategy")) +
      geom_bar(position= position_dodge(),stat='identity',colour = "black", width = 0.8) +
      facet_wrap(~ Simulation, ncol=3, scales='fixed') +
      ylab("mESS/second (sampling time)") + xlab("") + 
      theme(legend.position = "none") +
      coord_flip() +
      scale_x_discrete(limits = rev(levels(dfParametricEffSampling$Strategy))) +
      scale_fill_manual(values = labelData$colors) 
p

ggsave(filename = "figures/fig3b_simulation_efficiencies_sampling.png", plot = p,
        width = plot_width, height = plot_height, scale = 1.2, dpi = 300, units = unit, device='png')

##-----------------------------------------#
## Efficiency fof bnp vs parametric - simulation
##-----------------------------------------#
bnpUnimodal <- read.table(paste0("output/mcmc_time/simulation_unimodal/", bnpFileName), header = T)
bnpBimodal  <- read.table(paste0("output/mcmc_time/simulation_bimodal/", bnpFileName), header = T)
bnpMultimodal  <- read.table(paste0("output/mcmc_time/simulation_multimodal/", bnpFileName), header = T)

bnpUnimodal$ESS_second <- bnpUnimodal$multiEssItemsAbility/bnpUnimodal$runningTime
bnpBimodal$ESS_second <- bnpBimodal$multiEssItemsAbility/bnpBimodal$runningTime
bnpMultimodal$ESS_second <- bnpMultimodal$multiEssItemsAbility/bnpMultimodal$runningTime


bnpBimodal$simulation <- "Bimodal simulation"
bnpUnimodal$simulation <- "Unimodal simulation"
bnpMultimodal$simulation <- "Multimodal simulation"

bnpUnimodal$labels <- gsub("bnp_", "", bnpUnimodal$fileName)
bnpBimodal$labels <- gsub("bnp_", "", bnpBimodal$fileName)
bnpMultimodal$labels <- gsub("bnp_", "", bnpMultimodal$fileName)

## match R labels to plot labels
bnpUnimodal$labels <- droplevels(labelData[match(bnpUnimodal$labels, labelData$R_label), ]$plot_label)
bnpBimodal$labels  <- droplevels(labelData[match(bnpBimodal$labels, labelData$R_label), ]$plot_label)
bnpMultimodal$labels  <- droplevels(labelData[match(bnpMultimodal$labels, labelData$R_label), ]$plot_label)

## select parametric models with bnp equivalent and create data frame for plotting
dfParametricBnp <- data.frame(rbind(unimodal[which(unimodal$label %in% levels(bnpUnimodal$labels)), 
							c("labels", "ESS_second", "simulation")], 
						   bimodal[which(bimodal$label %in% levels(bnpBimodal$labels)), 
							c("labels", "ESS_second", "simulation")])) 

dfParametricBnp <- rbind(dfParametricBnp, multimodal[which(multimodal$label %in% levels(bnpMultimodal$labels)), 
                                          c("labels", "ESS_second", "simulation")])

dfParametricBnp <- droplevels(dfParametricBnp)
dfParametricBnp$model <- "Parametric"

bnpUnimodal$model <- "Semiparametric"						  
bnpBimodal$model  <- "Semiparametric"						  
bnpMultimodal$model  <- "Semiparametric"                                     

dfParametricBnp <- rbind(dfParametricBnp, 
                         bnpUnimodal[, c("labels", "ESS_second", "simulation", "model")],
 				 bnpBimodal[, c("labels", "ESS_second", "simulation", "model")], 
                         bnpMultimodal[, c("labels", "ESS_second", "simulation", "model")])


colnames(dfParametricBnp) <- c("Strategy", "ESS", "Simulation", "Model")
dfParametricBnp$Strategy  <- droplevels(factor(dfParametricBnp$Strategy, levels = labelData$plot_label))
## set level order
dfParametricBnp$Simulation <- factor(dfParametricBnp$Simulation, levels = c("Unimodal simulation", "Bimodal simulation", "Multimodal simulation"))

colorsParametricBnp <- labelData$colors[match(levels(dfParametricBnp$Strategy), labelData$plot_label)]
##-----------------------------------------#
## Plot Figure 3a
##-----------------------------------------#

ylabel <- 'Effective sample size per second'
title <- paste0("Multivariate effective sample size per second (total time)")

p <-  ggplot(dfParametricBnp,  aes_string(x = "Strategy", y= "ESS", fill = "Strategy")) +
      geom_bar(position= position_dodge(),stat='identity',colour = "black",
       width = 0.8) +
      facet_wrap(~ Model + Simulation,ncol=3, scales='fixed') +
      ylab("mESS/second (total time)") + xlab("") + 
      theme(legend.position = "none") +
      coord_flip() +
      scale_x_discrete(limits = rev(levels(dfParametricBnp$Strategy))) +
      scale_fill_manual(values = colorsParametricBnp) 
p

ggsave(filename = "figures/fig4_simulation_efficiencies_bnp.png", plot = p,
        width = plot_width, height = plot_height/2*3, scale = 1.2,  dpi = 300, units = unit, device='png')

##-----------------------------------------#
## Real data - efficiency comparisons
##-----------------------------------------#
paraHealth <- read.table(paste0("output/mcmc_time/data_health/", paraFileName), header = T)
paraTimss  <- read.table(paste0("output/mcmc_time/data_timss/", paraFileName), header = T)
paraTimss3PL  <- read.table(paste0("output/mcmc_time/data_timss/parametric3PL_efficiency.txt"), header = T)

# paraTimss <- droplevels(paraTimss[!is.na(paraTimss$labels), ])

paraHealth$model <- "Parametric 2PL"             
paraTimss$model  <- "Parametric 2PL"             
paraTimss3PL$model <- "Parametric 3PL"


paraHealth$ESS_second <- paraHealth$multiEssItemsAbility/paraHealth$runningTime
paraTimss$ESS_second  <- paraTimss$multiEssItemsAbility/paraTimss$runningTime
paraTimss3PL$ESS_second <- paraTimss3PL$multiEssItemsAbility/paraTimss3PL$runningTime

paraHealth$simulation <- "Health data"
paraTimss$simulation  <- "TIMSS data"
paraTimss3PL$simulation <- "TIMSS data"

paraHealth$labels <- gsub("parametric_", "", paraHealth$fileName)
paraTimss$labels  <- gsub("parametric_", "", paraTimss$fileName)
paraTimss3PL$labels <- gsub("parametric3PL_", "", paraTimss3PL$fileName)

## match R labels to plot labels
paraHealth$labels <- labelData[match(paraHealth$labels, labelData$R_label), ]$plot_label
paraTimss$labels  <- labelData[match(paraTimss$labels, labelData$R_label), ]$plot_label
paraTimss3PL$labels  <- labelData[match(paraTimss3PL$labels, labelData$R_label), ]$plot_label

# paraHealth$model <- "Parametric"             
# paraTimss$model  <- "Parametric"             
paraHealth$xmin <- 0
paraTimss$xmin  <- 0
paraTimss3PL$xmin <- 0

paraHealth$xmax <- 5
paraTimss$xmax  <- 2
paraTimss3PL$xmax <- 1.5

## data frame for plotting
dfParametricEff <- data.frame(rbind(paraHealth[, c("labels", "ESS_second", "simulation", "model", "xmin", "xmax")],
                                    paraTimss[, c("labels", "ESS_second", "simulation","model", "xmin", "xmax")], 
                                    paraTimss3PL[, c("labels", "ESS_second", "simulation","model", "xmin", "xmax")]))

colnames(dfParametricEff) <- c("Strategy", "ESS", "Simulation", "Model", "xmin", "xmax")

dfParametricEff$Strategy  <- droplevels(factor(dfParametricEff$Strategy, levels = labelData$plot_label))
dfParametricEff <- droplevels(dfParametricEff[-grep("constrained item", dfParametricEff$Strategy), ])

## Remove Stan results because of variability in the mESS
dfParametricEff <- droplevels(dfParametricEff[-which(dfParametricEff$Strategy == "IRT HMC (stan)**"), ])


colorsParametricData <- labelData$colors[match(levels(dfParametricEff$Strategy), labelData$plot_label)]
##-----------------------------------------#
## Plot Figure 5
##-----------------------------------------#
ylab <- paste0("mESS/second (total time)")
p <-  ggplot(dfParametricEff,  aes_string(x = "Strategy", y= "ESS", fill = "Strategy")) +
      geom_bar(position= position_dodge(),stat='identity',colour = "black",
       width = 0.8) +
      facet_wrap(~ Simulation + Model, ncol=3, scales='free_x') +
      geom_blank(aes(y = xmin)) +
      geom_blank(aes(y = xmax)) + 
      ylab("mESS/second (total time)") + xlab("") + 
      theme(legend.position = "none") +
      coord_flip() +
      scale_fill_manual(values = colorsParametricData) +
      scale_x_discrete(limits = rev(levels(dfParametricEff$Strategy))) 

p

ggsave(filename = "figures/fig5_data_efficiencies.png", plot = p,
        width = plot_width, height = plot_height/6*5 ,scale = 1.2, dpi = 300, units = unit, device='png')


##-----------------------------------------#
## Compare bnp efficiencies - data
##-----------------------------------------#

bnpHealth <- read.table(paste0("output/mcmc_time/data_health/", bnpFileName), header = T)
bnpTimss  <- read.table(paste0("output/mcmc_time/data_timss/", bnpFileName), header = T)
bnpTimss3PL  <- read.table(paste0("output/mcmc_time/data_timss/bnp3PL_efficiency.txt"), header = T)

bnpHealth$ESS_second <- bnpHealth$multiEssItemsAbility/bnpHealth$runningTime
bnpTimss$ESS_second <- bnpTimss$multiEssItemsAbility/bnpTimss$runningTime
bnpTimss3PL$ESS_second <- bnpTimss3PL$multiEssItemsAbility/bnpTimss3PL$runningTime

bnpHealth$simulation <- "Health data"
bnpTimss$simulation  <- "TIMSS data"
bnpTimss3PL$simulation  <- "TIMSS data"

bnpHealth$labels <- gsub("bnp_", "", bnpHealth$fileName)
bnpTimss$labels <- gsub("bnp_", "", bnpTimss$fileName)
bnpTimss3PL$labels <- gsub("bnp3PL_", "", bnpTimss3PL$fileName)

## match R labels to plot labels
bnpHealth$labels <- labelData[match(bnpHealth$labels, labelData$R_label), ]$plot_label
bnpTimss$labels  <- labelData[match(bnpTimss$labels, labelData$R_label), ]$plot_label
bnpTimss3PL$labels  <- labelData[match(bnpTimss3PL$labels, labelData$R_label), ]$plot_label

bnpHealth$model <- "Semiparametric 2PL"             
bnpTimss$model  <- "Semiparametric 2PL"             
bnpTimss3PL$model  <- "Semiparametric 3PL"        

bnpHealth$xmin <- 0
bnpTimss$xmin  <- 0
bnpTimss3PL$xmin <- 0     

bnpHealth$xmax <- 5
bnpTimss$xmax  <- 2
bnpTimss3PL$xmax <- 1.5

## data frame for plotting
dfBnpEff <- data.frame(rbind(bnpHealth[, c("labels", "ESS_second", "simulation", "model", "xmin", "xmax")],
                             bnpTimss[, c("labels", "ESS_second", "simulation","model", "xmin", "xmax")], 
                            bnpTimss3PL[, c("labels", "ESS_second", "simulation","model", "xmin", "xmax")]))

colnames(dfBnpEff) <- c("Strategy", "ESS", "Simulation", "Model", "xmin", "xmax")

dfBnpEff$Strategy  <- droplevels(factor(dfBnpEff$Strategy, levels = labelData$plot_label))
dfBnpEff <- droplevels(dfBnpEff[-grep("constrained item", dfBnpEff$Strategy), ])

colorsBnp <- labelData$colors[match(levels(dfBnpEff$Strategy), labelData$plot_label)]
##-----------------------------------------#
## Plot Figure 5b
##-----------------------------------------#


ylabel <- 'Effective sample size per second'
title <- paste0("Multivariate effective sample size per second (total time)")

p <-  ggplot(dfBnpEff,  aes_string(x = "Strategy", y= "ESS", fill = "Strategy", group = "Model")) +
      geom_bar(position= position_dodge(),stat='identity',colour = "black",
       width = 0.8) +
      facet_wrap(~ Simulation + Model ,ncol=3, scales='free_x') +
      geom_blank(aes(y = xmin)) +
      geom_blank(aes(y = xmax)) + 
      ylab("mESS/second (total time)") + xlab("") + 
      theme(legend.position = "none") +
      coord_flip() +
      scale_x_discrete(limits = rev(levels(dfBnpEff$Strategy))) +
      scale_fill_manual(values = colorsBnp) 
p

ggsave(filename = "figures/fig5b_data_efficiencies_bnp.png", plot = p,
        width = plot_width, height = plot_height/6*4 ,scale = 1.2, dpi = 300, units = unit, device='png')
