##-----------------------------------------#
## Computational strategies and estimation performance with Bayesian semiparametric Item Response Theory model
## Sally Paganin
## last update: August 2022
## R version 4.2.0 (2022-04-22) -- "Vigorous Calisthenics"
## nimble version 0.12.2
##-----------------------------------------#
source("R_functions/ggplot_settings.R")
library(cowplot)
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
fileList <- list.files("output/mcmc_time", full.names = TRUE)

##-----------------------------------------#
## Unimodal - extraSimulation
##-----------------------------------------#

unimodalFiles <- fileList[grep("unimodal", fileList)]

unimodalList <- list()
for(i in 1:length(unimodalFiles)){
      unimodalList[[i]] <- read.table(paste0(unimodalFiles[i], "/", paraFileName), header = TRUE )
      unimodalList[[i]]$simulation <- strsplit(unimodalFiles[i], "\\/")[[1]][3]
}     


unimodalDf <- as.data.frame(do.call(rbind, unimodalList))

# unique(unimodalDf$fileName)
# unique(labelData$R_label)

# unimodalDf$ESS_second <- unimodalDf$essCodaLogLik/unimodalDf$runningTime
# unimodalDf$ESS_second <- unimodalDf$essCodaLogPostItemsAbility/unimodalDf$runningTime
unimodalDf$ESS_second <- unimodalDf$multiEssItemsAbility/unimodalDf$runningTime

unimodalDf$labels <- gsub("parametric_", "", unimodalDf$fileName)
## match R labels to plot labels
unimodalDf$labels <- labelData[match(unimodalDf$labels, labelData$R_label), ]$plot_label
unimodalDf <- droplevels(unimodalDf[!is.na(unimodalDf$labels), ])

## data frame for plotting
dfParametricEff <- data.frame(unimodalDf[, c("labels", "ESS_second", "simulation")])

colnames(dfParametricEff) <- c("Strategy", "ESS", "Simulation")
dfParametricEff$Strategy  <- factor(dfParametricEff$Strategy, levels = labelData$plot_label)
dfParametricEff$Simulation <- factor(dfParametricEff$Simulation)

levels(dfParametricEff$Simulation) <-  c("N = 2000, I = 15", 
                                        "N = 1000, I = 10", 
                                        "N = 5000, I = 10", 
                                        "N = 1000, I = 30", 
                                        "N = 5000, I = 30")
dfParametricEff$Simulation <- factor(dfParametricEff$Simulatio, levels = levels(dfParametricEff$Simulation)[c(5,3, 1,4,2)])

## Remove Stan results because of variability in the mESS
dfParametricEff <- droplevels(dfParametricEff[-which(dfParametricEff$Strategy == "IRT HMC (stan)"), ])

yLabel <- paste0("mESS/second (total time)")

pUni <- ggplot(dfParametricEff, 
      aes(x = Simulation, y= ESS)) +
  geom_line(aes(group = Strategy, color = Strategy), position = position_dodge(width = 0.2)) +
  geom_point(aes(group = Strategy, color = Strategy), position = position_dodge(width = 0.2)) +
  ylab(yLabel) + xlab("") + 
  theme_bw() + 
  ylim(c(0,110)) +
  scale_color_manual(values = labelData$colors[-1]) +
  theme(legend.position = "bottom") + coord_flip() + 
  ggtitle("Unimodal Simulation")

pUni

ggsave(filename = "figures/SM_unimodalMultiESS.png", plot = pUni,
        width = plot_width/3, height = plot_height , 
        dpi = 300, scale = 1.4, units = unit, device='png')

##-----------------------------------------#
## Bimodal
##-----------------------------------------#

bimodalFiles <- fileList[grep("bimodal", fileList)]

bimodalList <- list()
for(i in 1:length(bimodalFiles)){
      bimodalList[[i]] <- read.table(paste0(bimodalFiles[i], "/", paraFileName), header = TRUE )
      bimodalList[[i]]$simulation <- strsplit(bimodalFiles[i], "\\/")[[1]][3]

}     

bimodalDf <- as.data.frame(do.call(rbind, bimodalList))

bimodalDf$ESS_second <- bimodalDf$multiEssItemsAbility/bimodalDf$runningTime

bimodalDf$labels <- gsub("parametric_", "", bimodalDf$fileName)
## match R labels to plot labels
bimodalDf$labels <- labelData[match(bimodalDf$labels, labelData$R_label), ]$plot_label
bimodalDf <- droplevels(bimodalDf[!is.na(bimodalDf$labels), ])

## data frame for plotting
dfParametricEff <- data.frame(bimodalDf[, c("labels", "ESS_second", "simulation")])

colnames(dfParametricEff) <- c("Strategy", "ESS", "Simulation")
dfParametricEff$Strategy  <- factor(dfParametricEff$Strategy, levels = labelData$plot_label)
dfParametricEff$Simulation <- factor(dfParametricEff$Simulation)

levels(dfParametricEff$Simulation) <-  c("N = 2000, I = 15", 
                                        "N = 1000, I = 10", 
                                        "N = 5000, I = 10", 
                                        "N = 1000, I = 30", 
                                        "N = 5000, I = 30")
dfParametricEff$Simulation <- factor(dfParametricEff$Simulation, levels = levels(dfParametricEff$Simulation)[c(5,3, 1,4,2)])

## Remove Stan results because of variability in the mESS
dfParametricEff <- droplevels(dfParametricEff[-which(dfParametricEff$Strategy == "IRT HMC (stan)"), ])

yLabel <- paste0("mESS/second (total time)")

pBi <- ggplot(dfParametricEff, 
      aes(x = Simulation, y= ESS)) +
  geom_line(aes(group = Strategy, color = Strategy), position = position_dodge(width = 0.2)) +
  geom_point(aes(group = Strategy, color = Strategy), position = position_dodge(width = 0.2)) +
  ylab(yLabel) + xlab("") + 
  theme_bw() + 
  ylim(c(0,110)) +
  scale_color_manual(values = labelData$colors[-1]) +
  theme(legend.position = "bottom") + coord_flip() + 
  ggtitle("Bimodal Simulation")

pBi

ggsave(filename = "figures/SM_bimodalMultiESS.png", plot = pBi,
        width = plot_width/3, height = plot_height , 
        dpi = 300, scale = 1.4, units = unit, device='png')

##-----------------------------------------#
## multimodal
##-----------------------------------------#

multimodalFiles <- fileList[grep("multimodal", fileList)]

multimodalList <- list()
for(i in 1:length(multimodalFiles)){
      multimodalList[[i]] <- read.table(paste0(multimodalFiles[i], "/", paraFileName), header = TRUE )
      multimodalList[[i]]$simulation <- strsplit(multimodalFiles[i], "\\/")[[1]][3]

}     

multimodalDf <- as.data.frame(do.call(rbind, multimodalList))

unique(multimodalDf$fileName)
unique(labelData$R_label)

multimodalDf$ESS_second <- multimodalDf$multiEssItemsAbility/multimodalDf$runningTime

multimodalDf$labels <- gsub("parametric_", "", multimodalDf$fileName)
## match R labels to plot labels
multimodalDf$labels <- labelData[match(multimodalDf$labels, labelData$R_label), ]$plot_label
multimodalDf <- droplevels(multimodalDf[!is.na(multimodalDf$labels), ])


## data frame for plotting
dfParametricEff <- data.frame(multimodalDf[, c("labels", "ESS_second", "simulation")])

colnames(dfParametricEff) <- c("Strategy", "ESS", "Simulation")
dfParametricEff$Strategy  <- factor(dfParametricEff$Strategy, levels = labelData$plot_label)
dfParametricEff$Simulation <- factor(dfParametricEff$Simulation)

levels(dfParametricEff$Simulation) <-  c("N = 2000, I = 15", 
                                        "N = 1000, I = 10", 
                                        "N = 5000, I = 10", 
                                        "N = 1000, I = 30", 
                                        "N = 5000, I = 30")

dfParametricEff$Simulation <- factor(dfParametricEff$Simulation, levels = levels(dfParametricEff$Simulation)[c(5,3, 1,4,2)])

## Remove Stan results because of variability in the mESS
dfParametricEff <- droplevels(dfParametricEff[-which(dfParametricEff$Strategy == "IRT HMC (stan)"), ])


pMulti <- ggplot(dfParametricEff, 
      aes(x = Simulation, y= ESS)) +
  geom_line(aes(group = Strategy, color = Strategy), position = position_dodge(width = 0.2)) +
  geom_point(aes(group = Strategy, color = Strategy), position = position_dodge(width = 0.2)) +
  ylab(yLabel) + xlab("") + 
  theme_bw() + 
  ylim(c(0,110)) +
  scale_color_manual(values = labelData$colors[-1]) +
  theme(legend.position = "bottom") + coord_flip() + 
  ggtitle("Multimodal Simulation")

pMulti

ggsave(filename = "figures/SM_multimodalMultiESS.png", plot = pMulti,
        width = plot_width/3, height = plot_height , 
        dpi = 300, scale = 1.4, units = unit, device='png')
##-----------------------------------------#


allPlotsTMP<- plot_grid(
    plot_grid(
    pUni  + theme(legend.position = "none"),
    pBi  + theme(legend.position = "none"),
    pMulti  + theme(legend.position = "none"),
    nrow = 1, align = "h"
   ),
   get_legend(pUni + theme(legend.position = "bottom")), 
   rel_heights = c(1, .1), nrow=2)

# title <- ggdraw() + draw_label("Parametric 2PL model", 
#       fontface='bold')

#allPlots <-  plot_grid( title, allPlotsTMP, rel_heights = c(.1, 1), nrow=2)

ggsave(filename = "figures/SM_fig1_allScenarioMultiESS.png", 
        plot = allPlotsTMP,
        width = plot_width, height = plot_height , 
        dpi = 300, scale = 1.4, units = unit, device='png')


##-----------------------------------------#
##-----------------------------------------#
## min ESS 
##-----------------------------------------#

unimodalFiles <- fileList[grep("unimodal", fileList)]

unimodalList <- list()
for(i in 1:length(unimodalFiles)){
      unimodalList[[i]] <- read.table(paste0(unimodalFiles[i], "/", paraFileName), header = TRUE )
      unimodalList[[i]]$simulation <- strsplit(unimodalFiles[i], "\\/")[[1]][3]
}     


unimodalDf <- as.data.frame(do.call(rbind, unimodalList))

unimodalDf$ESS_second <- unimodalDf$essCodaItemsAbility/unimodalDf$runningTime

unimodalDf$labels <- gsub("parametric_", "", unimodalDf$fileName)
## match R labels to plot labels
unimodalDf$labels <- labelData[match(unimodalDf$labels, labelData$R_label), ]$plot_label
unimodalDf <- droplevels(unimodalDf[!is.na(unimodalDf$labels), ])

## data frame for plotting
dfParametricEff <- data.frame(unimodalDf[, c("labels", "ESS_second", "simulation")])

colnames(dfParametricEff) <- c("Strategy", "ESS", "Simulation")
dfParametricEff$Strategy  <- factor(dfParametricEff$Strategy, levels = labelData$plot_label)
dfParametricEff$Simulation <- factor(dfParametricEff$Simulation)

levels(dfParametricEff$Simulation) <-  c("N = 2000, I = 15", 
                                        "N = 1000, I = 10", 
                                        "N = 5000, I = 10", 
                                        "N = 1000, I = 30", 
                                        "N = 5000, I = 30")
dfParametricEff$Simulation <- factor(dfParametricEff$Simulatio, levels = levels(dfParametricEff$Simulation)[c(5,3, 1,4,2)])

## Remove Stan results because of variability in the mESS
dfParametricEff <- droplevels(dfParametricEff[-which(dfParametricEff$Strategy == "IRT HMC (stan)"), ])

yLabel <- paste0("minESS/second (total time)")

dfParametricEff$Simulation

pUni <- ggplot(dfParametricEff, 
      aes(x = Simulation, y= ESS)) +
  geom_line(aes(group = Strategy, color = Strategy), position = position_dodge(width = 0.2)) +
  geom_point(aes(group = Strategy, color = Strategy), position = position_dodge(width = 0.2)) +
  ylab(yLabel) + xlab("") + 
  theme_bw() + 
  ylim(c(0,5)) +
  scale_color_manual(values = labelData$colors[-1]) +
  theme(legend.position = "bottom") + coord_flip() + 
  ggtitle("Unimodal Simulation")

pUni

ggsave(filename = "figures/SM_unimodalminESS.png", plot = pUni,
        width = plot_width/3, height = plot_height , 
        dpi = 300, scale = 1.4, units = unit, device='png')

##-----------------------------------------#
## Bimodal
##-----------------------------------------#

bimodalFiles <- fileList[grep("bimodal", fileList)]

bimodalList <- list()
for(i in 1:length(bimodalFiles)){
      bimodalList[[i]] <- read.table(paste0(bimodalFiles[i], "/", paraFileName), header = TRUE )
      bimodalList[[i]]$simulation <- strsplit(bimodalFiles[i], "\\/")[[1]][3]

}     

bimodalDf <- as.data.frame(do.call(rbind, bimodalList))

bimodalDf$ESS_second <- bimodalDf$essCodaItemsAbility/bimodalDf$runningTime

bimodalDf$labels <- gsub("parametric_", "", bimodalDf$fileName)
## match R labels to plot labels
bimodalDf$labels <- labelData[match(bimodalDf$labels, labelData$R_label), ]$plot_label
bimodalDf <- droplevels(bimodalDf[!is.na(bimodalDf$labels), ])

## data frame for plotting
dfParametricEff <- data.frame(bimodalDf[, c("labels", "ESS_second", "simulation")])

colnames(dfParametricEff) <- c("Strategy", "ESS", "Simulation")
dfParametricEff$Strategy  <- factor(dfParametricEff$Strategy, levels = labelData$plot_label)
dfParametricEff$Simulation <- factor(dfParametricEff$Simulation)

levels(dfParametricEff$Simulation) <-  c("N = 2000, I = 15", 
                                        "N = 1000, I = 10", 
                                        "N = 5000, I = 10", 
                                        "N = 1000, I = 30", 
                                        "N = 5000, I = 30")
dfParametricEff$Simulation <- factor(dfParametricEff$Simulation, levels = levels(dfParametricEff$Simulation)[c(5,3, 1,4,2)])

## Remove Stan results because of variability in the mESS
dfParametricEff <- droplevels(dfParametricEff[-which(dfParametricEff$Strategy == "IRT HMC (stan)"), ])

yLabel <- paste0("minESS/second (total time)")

pBi <- ggplot(dfParametricEff, 
      aes(x = Simulation, y= ESS)) +
  geom_line(aes(group = Strategy, color = Strategy), position = position_dodge(width = 0.2)) +
  geom_point(aes(group = Strategy, color = Strategy), position = position_dodge(width = 0.2)) +
  ylab(yLabel) + xlab("") + 
  theme_bw() + 
  ylim(c(0,5)) +
  scale_color_manual(values = labelData$colors[-1]) +
  theme(legend.position = "bottom") + coord_flip() + 
  ggtitle("Bimodal Simulation")

pBi

ggsave(filename = "figures/SM_bimodalminESS.png", plot = pBi,
        width = plot_width/3, height = plot_height , 
        dpi = 300, scale = 1.4, units = unit, device='png')


##-----------------------------------------#
## multimodal
##-----------------------------------------#

multimodalFiles <- fileList[grep("multimodal", fileList)]

multimodalList <- list()
for(i in 1:length(multimodalFiles)){
      multimodalList[[i]] <- read.table(paste0(multimodalFiles[i], "/", paraFileName), header = TRUE )
      multimodalList[[i]]$simulation <- strsplit(multimodalFiles[i], "\\/")[[1]][3]

}     

multimodalDf <- as.data.frame(do.call(rbind, multimodalList))

unique(multimodalDf$fileName)
unique(labelData$R_label)

multimodalDf$ESS_second <- multimodalDf$essCodaItemsAbility/multimodalDf$runningTime

multimodalDf$labels <- gsub("parametric_", "", multimodalDf$fileName)
## match R labels to plot labels
multimodalDf$labels <- labelData[match(multimodalDf$labels, labelData$R_label), ]$plot_label
multimodalDf <- droplevels(multimodalDf[!is.na(multimodalDf$labels), ])


## data frame for plotting
dfParametricEff <- data.frame(multimodalDf[, c("labels", "ESS_second", "simulation")])


colnames(dfParametricEff) <- c("Strategy", "ESS", "Simulation")
dfParametricEff$Strategy  <- factor(dfParametricEff$Strategy, levels = labelData$plot_label)
dfParametricEff$Simulation <- factor(dfParametricEff$Simulation)
## Remove Stan results because of variability in the mESS
dfParametricEff <- droplevels(dfParametricEff[-which(dfParametricEff$Strategy == "IRT HMC (stan)"), ])

levels(dfParametricEff$Simulation) <-  c("N = 2000, I = 15", 
                                        "N = 1000, I = 10", 
                                        "N = 5000, I = 10", 
                                        "N = 1000, I = 30", 
                                        "N = 5000, I = 30")

dfParametricEff$Simulation <- factor(dfParametricEff$Simulation, levels = levels(dfParametricEff$Simulation)[c(5,3, 1,4,2)])

pMulti <- ggplot(dfParametricEff, 
      aes(x = Simulation, y= ESS)) +
  geom_line(aes(group = Strategy, color = Strategy), position = position_dodge(width = 0.2)) +
  geom_point(aes(group = Strategy, color = Strategy), position = position_dodge(width = 0.2)) +
  ylab(yLabel) + xlab("") + 
  theme_bw() + 
  ylim(c(0,5)) +
  scale_color_manual(values = labelData$colors[-1]) +
  theme(legend.position = "bottom") + coord_flip() + 
  ggtitle("Multimodal Simulation")

pMulti

ggsave(filename = "figures/SM_multimodalminESS.png", plot = pMulti,
        width = plot_width/3, height = plot_height , 
        dpi = 300, scale = 1.4, units = unit, device='png')
##-----------------------------------------#


allPlotsTMP <- plot_grid(
    plot_grid(
    pUni  + theme(legend.position = "none"),
    pBi  + theme(legend.position = "none"),
    pMulti  + theme(legend.position = "none"),
    nrow = 1, align = "h"
   ),
   get_legend(pUni + theme(legend.position = "bottom")), 
   rel_heights = c(1, .1), nrow=2)

# title <- ggdraw() + draw_label("Parametric 2PL model", 
#       fontface='bold')

# allPlots <- plot_grid( title, allPlotsTMP, rel_heights = c(.1, 1), nrow=2)


ggsave(filename = "figures/SM_allScenarioMinESS.png", 
        plot = allPlotsTMP,
        width = plot_width, height = plot_height , 
        dpi = 300, scale = 1.4, units = unit, device='png')

##-----------------------------------------###-----------------------------------------#
##-----------------------------------------#
## Unimodal - extraSimulation - bnp
##-----------------------------------------#

unimodalFiles <- fileList[grep("unimodal", fileList)]

unimodalList <- list()
for(i in 1:length(unimodalFiles)){
      unimodalList[[i]] <- read.table(paste0(unimodalFiles[i], "/", bnpFileName), header = TRUE )
      unimodalList[[i]]$simulation <- strsplit(unimodalFiles[i], "\\/")[[1]][3]
}     


unimodalDf <- as.data.frame(do.call(rbind, unimodalList))
unimodalDf$ESS_second <- unimodalDf$multiEssItemsAbility/unimodalDf$runningTime

unimodalDf$labels <- gsub("bnp_", "", unimodalDf$fileName)
## match R labels to plot labels
unimodalDf$labels <- labelData[match(unimodalDf$labels, labelData$R_label), ]$plot_label
unimodalDf <- droplevels(unimodalDf[!is.na(unimodalDf$labels), ])

## data frame for plotting
dfParametricEff <- data.frame(unimodalDf[, c("labels", "ESS_second", "simulation")])

colnames(dfParametricEff) <- c("Strategy", "ESS", "Simulation")
dfParametricEff$Strategy  <- factor(dfParametricEff$Strategy, levels = labelData$plot_label)
dfParametricEff$Simulation <- factor(dfParametricEff$Simulation)

levels(dfParametricEff$Simulation) <-  c("N = 2000, I = 15", 
                                        "N = 1000, I = 10", 
                                        "N = 5000, I = 10", 
                                        "N = 1000, I = 30", 
                                        "N = 5000, I = 30")
dfParametricEff$Simulation <- factor(dfParametricEff$Simulatio, levels = levels(dfParametricEff$Simulation)[c(5,3, 1,4,2)])

yLabel <- paste0("mESS/second (total time)")

pUni <- ggplot(dfParametricEff, 
      aes(x = Simulation, y= ESS)) +
  geom_line(aes(group = Strategy, color = Strategy), position = position_dodge(width = 0.2)) +
  geom_point(aes(group = Strategy, color = Strategy), position = position_dodge(width = 0.2)) +
  ylab(yLabel) + xlab("") + 
  theme_bw() + 
  ylim(c(0,40)) +
  scale_color_manual(values = labelData$colors[-1]) +
  theme(legend.position = "bottom") + coord_flip() + 
  ggtitle("Unimodal Simulation")

pUni

ggsave(filename = "figures/SM_BNP_unimodalMultiESS.png", plot = pUni,
        width = plot_width/3, height = plot_height , 
        dpi = 300, scale = 1.4, units = unit, device='png')

##-----------------------------------------#
## Bimodal
##-----------------------------------------#

bimodalFiles <- fileList[grep("bimodal", fileList)]

bimodalList <- list()
for(i in 1:length(bimodalFiles)){
      bimodalList[[i]] <- read.table(paste0(bimodalFiles[i], "/", bnpFileName), header = TRUE )
      bimodalList[[i]]$simulation <- strsplit(bimodalFiles[i], "\\/")[[1]][3]

}     

bimodalDf <- as.data.frame(do.call(rbind, bimodalList))

bimodalDf$ESS_second <- bimodalDf$multiEssItemsAbility/bimodalDf$runningTime

bimodalDf$labels <- gsub("bnp_", "", bimodalDf$fileName)
## match R labels to plot labels
bimodalDf$labels <- labelData[match(bimodalDf$labels, labelData$R_label), ]$plot_label
bimodalDf <- droplevels(bimodalDf[!is.na(bimodalDf$labels), ])

## data frame for plotting
dfParametricEff <- data.frame(bimodalDf[, c("labels", "ESS_second", "simulation")])

colnames(dfParametricEff) <- c("Strategy", "ESS", "Simulation")
dfParametricEff$Strategy  <- factor(dfParametricEff$Strategy, levels = labelData$plot_label)
dfParametricEff$Simulation <- factor(dfParametricEff$Simulation)

levels(dfParametricEff$Simulation) <-  c("N = 2000, I = 15", 
                                        "N = 1000, I = 10", 
                                        "N = 5000, I = 10", 
                                        "N = 1000, I = 30", 
                                        "N = 5000, I = 30")
dfParametricEff$Simulation <- factor(dfParametricEff$Simulation, levels = levels(dfParametricEff$Simulation)[c(5,3, 1,4,2)])


yLabel <- paste0("mESS/second (total time)")

pBi <- ggplot(dfParametricEff, 
      aes(x = Simulation, y= ESS)) +
  geom_line(aes(group = Strategy, color = Strategy), position = position_dodge(width = 0.2)) +
  geom_point(aes(group = Strategy, color = Strategy), position = position_dodge(width = 0.2)) +
  ylab(yLabel) + xlab("") + 
  theme_bw() + 
  ylim(c(0,40)) +
  scale_color_manual(values = labelData$colors[-1]) +
  theme(legend.position = "bottom") + coord_flip() + 
  ggtitle("Bimodal Simulation")

pBi

ggsave(filename = "figures/SM_BNP_bimodalMultiESS.png", plot = pBi,
        width = plot_width/3, height = plot_height , 
        dpi = 300, scale = 1.4, units = unit, device='png')

##-----------------------------------------#
## multimodal
##-----------------------------------------#

multimodalFiles <- fileList[grep("multimodal", fileList)]

multimodalList <- list()
for(i in 1:length(multimodalFiles)){
      multimodalList[[i]] <- read.table(paste0(multimodalFiles[i], "/", bnpFileName), header = TRUE )
      multimodalList[[i]]$simulation <- strsplit(multimodalFiles[i], "\\/")[[1]][3]

}     

multimodalDf <- as.data.frame(do.call(rbind, multimodalList))

unique(multimodalDf$fileName)
unique(labelData$R_label)

multimodalDf$ESS_second <- multimodalDf$multiEssItemsAbility/multimodalDf$runningTime

multimodalDf$labels <- gsub("bnp_", "", multimodalDf$fileName)
## match R labels to plot labels
multimodalDf$labels <- labelData[match(multimodalDf$labels, labelData$R_label), ]$plot_label
multimodalDf <- droplevels(multimodalDf[!is.na(multimodalDf$labels), ])


## data frame for plotting
dfParametricEff <- data.frame(multimodalDf[, c("labels", "ESS_second", "simulation")])

colnames(dfParametricEff) <- c("Strategy", "ESS", "Simulation")
dfParametricEff$Strategy  <- factor(dfParametricEff$Strategy, levels = labelData$plot_label)
dfParametricEff$Simulation <- factor(dfParametricEff$Simulation)

levels(dfParametricEff$Simulation) <-  c("N = 2000, I = 15", 
                                        "N = 1000, I = 10", 
                                        "N = 5000, I = 10", 
                                        "N = 1000, I = 30", 
                                        "N = 5000, I = 30")

dfParametricEff$Simulation <- factor(dfParametricEff$Simulation, levels = levels(dfParametricEff$Simulation)[c(5,3, 1,4,2)])


pMulti <- ggplot(dfParametricEff, 
      aes(x = Simulation, y= ESS)) +
  geom_line(aes(group = Strategy, color = Strategy), position = position_dodge(width = 0.2)) +
  geom_point(aes(group = Strategy, color = Strategy), position = position_dodge(width = 0.2)) +
  ylab(yLabel) + xlab("") + 
  theme_bw() + 
  ylim(c(0,40)) +
  scale_color_manual(values = labelData$colors[-1]) +
  theme(legend.position = "bottom") + coord_flip() + 
  ggtitle("Multimodal Simulation")

pMulti

ggsave(filename = "figures/SM_BNP_multimodalMultiESS.png", plot = pMulti,
        width = plot_width/3, height = plot_height , 
        dpi = 300, scale = 1.4, units = unit, device='png')
##-----------------------------------------#

allPlotsTMP <- plot_grid(
    plot_grid(
    pUni  + theme(legend.position = "none"),
    pBi  + theme(legend.position = "none"),
    pMulti  + theme(legend.position = "none"),
    nrow = 1, align = "h"
   ),
   get_legend(pUni + theme(legend.position = "bottom")), 
   rel_heights = c(1, .1), nrow=2)

# title <- ggdraw() + draw_label("Semiparametric 2PL model", 
#       fontface='bold')

# allPlots <-  plot_grid( title, allPlotsTMP, rel_heights = c(.1, 1), nrow=2)


ggsave(filename = "figures/SM_fig2_BNP_allScenarioMultiESS.png", 
        plot = allPlotsTMP,
        width = plot_width, height = plot_height , 
        dpi = 300, scale = 1.4, units = unit, device='png')

##-----------------------------------------#
## ESS distribution plots
##-----------------------------------------#
## Unimodal

xx <- readRDS("output/posterior_samples_elaborated/simulation_unimodal/ESS_parametric_IRT_constrainedItem.rds")
unimodalDF <- data.frame(cbind(names(xx), xx, "IRT constrained item"))

xx <- readRDS("output/posterior_samples_elaborated/simulation_unimodal/ESS_parametric_IRT_unconstrained.rds")
unimodalDF <- rbind(unimodalDF, (cbind(names(xx), xx, "IRT unconstrained")))

xx <- readRDS("output/posterior_samples_elaborated/simulation_unimodal/ESS_parametric_IRT_stan.rds")
unimodalDF <- rbind(unimodalDF, (cbind(names(xx), xx, "IRT HMC (stan)**")))

colnames(unimodalDF) <- c("parameter", "efficiency", "model")

unimodalDF$parameter[grep("lambda", unimodalDF$parameter)] <- "lambda"
unimodalDF$parameter[grep("beta", unimodalDF$parameter)] <- "beta"
unimodalDF$parameter[grep("^eta", unimodalDF$parameter)] <- "eta"

unimodalDF$efficiency <- as.numeric(unimodalDF$efficiency)
unimodalDF$model <- factor(unimodalDF$model, levels = unique(unimodalDF$model)[c(3, 2, 1)])

unimodalDF$parameter <- factor(unimodalDF$parameter, levels = c("beta", "lambda", "eta"))


plotCol <- labelData[c(1, 4, 2), ]$colors


# ggplot(unimodalDF, aes(x = efficiency, fill = model)) + 
#       geom_histogram(aes(y = ..density..), color = "white", binwidth=0.2) + 
#       facet_grid(model ~ parameter , scales = "free") + 
#       scale_fill_manual(values = plotCol) +
#       xlab("ESS/time") + ylab("Density")

plotUni <- ggplot(unimodalDF, aes(x = efficiency, fill = model)) + 
      geom_histogram(data=subset(unimodalDF, parameter == "eta"), aes(y = ..density..), color = "white", binwidth=0.5) + 
      geom_histogram(data=subset(unimodalDF, parameter == "beta"), aes(y = ..density..), color = "white", binwidth=0.2) + 
      geom_histogram(data=subset(unimodalDF, parameter == "lambda"), aes(y = ..density..), color = "white", binwidth=0.2) + 
      facet_grid(model ~ parameter , scales = "free") + 
      scale_fill_manual(values = plotCol) +
      xlab("ESS/second") + ylab("Density") + theme(legend.position = 'none')
plotUni


# infoUni <- read.table(paste0("output/mcmc_time/simulation_unimodal/", paraFileName), header = T)
# infoUni <- subset(infoUni, fileName %in% c("parametric_IRT_stan","parametric_IRT_unconstrained", "parametric_IRT_constrainedItem"))

# infoUni <- infoUni[, c("fileName", "multiEssItemsAbility", "runningTime")]
# colnames(infoUni) <- c("model", "mESS", "totalTime")
# infoUni$model  <- levels(unimodalDF$model)[c(3, 2, 1)]
# infoUni$parameters <- "eta"
# infoUni$x <- rep(17, 3)
# infoUni$y <- c(3, 0.8,1)
# infoUni$text <- paste0("mESS/second: ", round(infoUni$mESS/infoUni$totalTime))

# infoUni$model <- factor(infoUni$model, levels = levels(unimodalDF$model))

# infoUni$parameter <- factor(infoUni$parameter, levels = c("beta", "lambda", "eta"))


# plotUni <- plotUni + geom_text(data = infoUni, 
#                    aes(x = x, y = y, label = text), size = 4)

# plotUni

ggsave(filename = "figures/SM_figESS_unimodal.png", 
        plot = plotUni,
        width = plot_width, height = plot_height*1.5 , 
        dpi = 300, scale = 1.4, units = unit, device='png')
######


xx <- readRDS("output/posterior_samples_elaborated/simulation_bimodal/ESS_parametric_IRT_constrainedItem.rds")
bimodalDF <- data.frame(cbind(names(xx), xx, "IRT constrained item"))

xx <- readRDS("output/posterior_samples_elaborated/simulation_bimodal/ESS_parametric_IRT_unconstrained.rds")
bimodalDF <- rbind(bimodalDF, (cbind(names(xx), xx, "IRT unconstrained")))

xx <- readRDS("output/posterior_samples_elaborated/simulation_bimodal/ESS_parametric_IRT_stan.rds")
bimodalDF <- rbind(bimodalDF, (cbind(names(xx), xx, "IRT HMC (stan)**")))

colnames(bimodalDF) <- c("parameter", "efficiency", "model")

bimodalDF$parameter[grep("lambda", bimodalDF$parameter)] <- "lambda"
bimodalDF$parameter[grep("beta", bimodalDF$parameter)] <- "beta"
bimodalDF$parameter[grep("^eta", bimodalDF$parameter)] <- "eta"

bimodalDF$efficiency <- as.numeric(bimodalDF$efficiency)
bimodalDF$model <- factor(bimodalDF$model, levels = unique(bimodalDF$model)[c(3, 2, 1)])

bimodalDF$parameter <- factor(bimodalDF$parameter, levels = c("beta", "lambda", "eta"))


plotCol <- labelData[c(1, 4, 2), ]$colors


# ggplot(bimodalDF, aes(x = efficiency, fill = model)) + 
#       geom_histogram(aes(y = ..density..), color = "white", binwidth=0.2) + 
#       facet_grid(model ~ parameter , scales = "free") + 
#       scale_fill_manual(values = plotCol) +
#       xlab("ESS/second") + ylab("Density")

plotBi <- ggplot(bimodalDF, aes(x = efficiency, fill = model)) + 
      geom_histogram(data=subset(bimodalDF, parameter == "eta"), aes(y = ..density..), color = "white", binwidth=0.2) + 
      geom_histogram(data=subset(bimodalDF, parameter == "beta"), aes(y = ..density..), color = "white", binwidth=0.2) + 
      geom_histogram(data=subset(bimodalDF, parameter == "lambda"), aes(y = ..density..), color = "white", binwidth=0.2) + 
      facet_grid(model ~ parameter , scales = "free") + 
      scale_fill_manual(values = plotCol) +
      xlab("ESS/second") + ylab("Density") + theme(legend.position = 'none')
plotBi


# infoBi <- read.table(paste0("output/mcmc_time/simulation_bimodal/", paraFileName), header = T)
# infoBi <- subset(infoBi, fileName %in% c("parametric_IRT_stan","parametric_IRT_unconstrained", "parametric_IRT_constrainedItem"))

# infoBi <- infoBi[, c("fileName", "multiEssItemsAbility", "runningTime")]
# colnames(infoBi) <- c("model", "mESS", "totalTime")
# infoBi$model  <- levels(unimodalDF$model)[c(3, 2, 1)]
# infoBi$parameters <- "eta"
# infoBi$x <- rep(5, 3)
# infoBi$y <- c(3, 0.8, 1)
# infoBi$text <- paste0("mESS/second: ", infoBi$mESS/infoBi$totalTime)

# infoBi$model <- factor(infoBi$model, levels = levels(unimodalDF$model))

# infoBi$parameter <- factor(infoBi$parameter, levels = c("beta", "lambda", "eta"))


# plotBi <- plotBi + geom_text(data = infoBi, 
#                    aes(x = x, y = y, label = text), size = 4)

ggsave(filename = "figures/SM_figESS_bimodal.png", 
        plot = plotBi,
        width = plot_width, height = plot_height*1.5 , 
        dpi = 300, scale = 1.4, units = unit, device='png')


######


xx <- readRDS("output/posterior_samples_elaborated/simulation_multimodal/ESS_parametric_IRT_constrainedItem.rds")
multimodalDF <- data.frame(cbind(names(xx), xx, "IRT constrained item"))

xx <- readRDS("output/posterior_samples_elaborated/simulation_multimodal/ESS_parametric_IRT_unconstrained.rds")
multimodalDF <- rbind(multimodalDF, (cbind(names(xx), xx, "IRT unconstrained")))

xx <- readRDS("output/posterior_samples_elaborated/simulation_multimodal/ESS_parametric_IRT_stan.rds")
multimodalDF <- rbind(multimodalDF, (cbind(names(xx), xx, "IRT HMC (stan)**")))


multimodalDF$ESS_second

colnames(multimodalDF) <- c("parameter", "efficiency", "model")

multimodalDF$parameter[grep("lambda", multimodalDF$parameter)] <- "lambda"
multimodalDF$parameter[grep("beta", multimodalDF$parameter)] <- "beta"
multimodalDF$parameter[grep("^eta", multimodalDF$parameter)] <- "eta"

multimodalDF$efficiency <- as.numeric(multimodalDF$efficiency)
multimodalDF$model <- factor(multimodalDF$model, levels = unique(multimodalDF$model)[c(3, 2, 1)])

multimodalDF$parameter <- factor(multimodalDF$parameter, levels = c("beta", "lambda", "eta"))


plotCol <- labelData[c(1, 4, 2), ]$colors



plotMulti <- ggplot(multimodalDF, aes(x = efficiency, fill = model)) + 
      geom_histogram(data=subset(multimodalDF, parameter == "eta"), aes(y = ..density..), color = "white", binwidth=0.3) + 
      geom_histogram(data=subset(multimodalDF, parameter == "beta"), aes(y = ..density..), color = "white", binwidth=0.2) + 
      geom_histogram(data=subset(multimodalDF, parameter == "lambda"), aes(y = ..density..), color = "white", binwidth=0.2) + 
      facet_grid(model ~ parameter , scales = "free") + 
      scale_fill_manual(values = plotCol) +
      xlab("ESS/time") + ylab("Density") + theme(legend.position = 'none')
plotMulti

# infoMulti <- read.table(paste0("output/mcmc_time/simulation_multimodal/", paraFileName), header = T)
# infoMulti <- subset(infoMulti, fileName %in% c("parametric_IRT_stan","parametric_IRT_unconstrained", "parametric_IRT_constrainedItem"))

# infoMulti <- infoMulti[, c("fileName", "multiEssItemsAbility", "runningTime")]
# colnames(infoMulti) <- c("model", "mESS", "totalTime")
# infoMulti$model  <- levels(unimodalDF$model)[c(3, 2, 1)]
# infoMulti$parameters <- "eta"
# infoMulti$x <- rep(11, 3)
# infoMulti$y <- c(1, 1, 0.8)
# infoMulti$text <- paste0("mESS/second: ", infoMulti$mESS/infoMulti$totalTime)

# infoMulti$model <- factor(infoMulti$model, levels = levels(unimodalDF$model))

# infoMulti$parameter <- factor(infoMulti$parameter, levels = c("beta", "lambda", "eta"))


# plotMulti <- plotMulti + geom_text(data = infoMulti, 
#                    aes(x = x, y = y, label = text), size = 4)

ggsave(filename = "figures/SM_figESS_multimodal.png", 
        plot = plotMulti,
        width = plot_width, height = plot_height*1.5 , 
        dpi = 300, scale = 1.4, units = unit, device='png')



