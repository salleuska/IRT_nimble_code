##-----------------------------------------#
## Bayesian semiparametric Item Response Theory models using NIMBLE 
## Sally Paganin
## November 2020
##-----------------------------------------#
library(ggplot2)
library(cowplot)
source("R_functions/ggplot_settings.R")
##-----------------------------------------#
## Set dimensions
plot_width  <- 21*1.2 ## a4 paper
plot_height <- (plot_width/5*2)*1.2
unit   <- "cm" 

##-----------------------------------------#
## Unimodal
##-----------------------------------------#
dataName <- "simulation_unimodal"
unimodalRes <- readRDS("figures/dataForFigures/unimodal.rds")


##-----------------------------------------#
## Discriminations
##-----------------------------------------#

estimateDiscr <- data.frame(estimate =  c(unimodalRes$paraEstimates$lambda,
                                           unimodalRes$bnpEstimates$lambda), 
                             CI_low = c(unimodalRes$paraLow$lambda,
                                        unimodalRes$bnpLow$lambda), 
                             CI_upp = c(unimodalRes$paraUpper$lambda,
                                        unimodalRes$bnpUpper$lambda))

estimateDiscr$trueVal <- rep(unimodalRes$truValues$lambda, 2)
estimateDiscr$Model <- rep(c("Parametric", "Semiparametric"), each = length(unimodalRes$truValues$lambda))

pUniDiscr <- ggplot(estimateDiscr, aes(x = trueVal, y = estimate, color = Model)) + 
        geom_point(size = 1.5) + 
          scale_x_continuous(breaks = seq(0.2, 1.6, by = 0.2)) +
          scale_y_continuous(breaks = seq(0.2, 1.6, by = 0.2)) + 
        geom_errorbar(aes(ymin=CI_low, ymax=CI_upp), width = 0.05) + 
        scale_color_manual(values=c(paraColor, bnpColor),
            guide = guide_legend(override.aes = list(linetype = rep("blank", 2), 
                                                     shape    = c(16, 16)))) +
        labs(x = "True value", y = "Estimate", title = "Discrimination parameters") + 
        geom_abline(color = 'grey10', lty = 3)

pUniDiscr
ggsave(filename = "figures/unimodal_discriminations.png", plot = pUniDiscr,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')

##-----------------------------------------#
## Difficulties
##-----------------------------------------#

estimateDiff <- data.frame(estimate =  c(unimodalRes$paraEstimates$beta,
                                         unimodalRes$bnpEstimates$beta), 
                             CI_low = c(unimodalRes$paraLow$beta,
                                        unimodalRes$bnpLow$beta), 
                             CI_upp = c(unimodalRes$paraUpper$beta,
                                        unimodalRes$bnpUpper$beta))

estimateDiff$trueVal <- rep(unimodalRes$truValues$beta, 2)
estimateDiff$Model <- rep(c("Parametric", "Semiparametric"), each = length(unimodalRes$truValues$beta))

pUniDiff <- ggplot(estimateDiff, aes(x = trueVal, y = estimate, color = Model)) + 
        geom_point(size = 1.5) + 
          scale_x_continuous(breaks = seq(-3, 3, by = 1)) +
          scale_y_continuous(breaks = seq(-3, 3, by = 1)) +
        geom_errorbar(aes(ymin=CI_low, ymax=CI_upp), width = 0.2) + 
        scale_color_manual(values=c(paraColor, bnpColor),
            guide = guide_legend(override.aes = list(linetype = rep("blank", 2), 
                                                     shape    = c(16, 16)))) +
        labs(x = "True value", y = "Estimate", title = "Difficulty parameters") + 
        geom_abline(color = 'grey10', lty = 3)

pUniDiff

ggsave(filename = "figures/unimodal_difficulties.png", plot = pUniDiff,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')

##-----------------------------------------#
## Density - mean 
##-----------------------------------------#


##-----------------------------------------#
## Bimodal
##-----------------------------------------#
bimodalRes <- readRDS("figures/dataForFigures/bimodal.rds")

##-----------------------------------------#
## Discriminations
##-----------------------------------------#

estimateDiscr <- data.frame(estimate =  c(bimodalRes$paraEstimates$lambda,
                                           bimodalRes$bnpEstimates$lambda), 
                             CI_low = c(bimodalRes$paraLow$lambda,
                                        bimodalRes$bnpLow$lambda), 
                             CI_upp = c(bimodalRes$paraUpper$lambda,
                                        bimodalRes$bnpUpper$lambda))

estimateDiscr$trueVal <- rep(bimodalRes$truValues$lambda, 2)
estimateDiscr$Model <- rep(c("Parametric", "Semiparametric"), each = length(bimodalRes$truValues$lambda))

pBiDiscr <- ggplot(estimateDiscr, aes(x = trueVal, y = estimate, color = Model)) + 
        geom_point(size = 1.5) + 
          scale_x_continuous(breaks = seq(0.2, 1.6, by = 0.2)) +
          scale_y_continuous(breaks = seq(0.2, 1.6, by = 0.2)) + 
        geom_errorbar(aes(ymin=CI_low, ymax=CI_upp), width = 0.05) + 
        scale_color_manual(values=c(paraColor, bnpColor),
            guide = guide_legend(override.aes = list(linetype = rep("blank", 2), 
                                                     shape    = c(16, 16)))) +
        labs(x = "True value", y = "Estimate", title = "Discrimination parameters") + 
        geom_abline(color = 'grey10', lty = 3)

pBiDiscr
ggsave(filename = "figures/bimodal_discriminations.png", plot = pBiDiscr,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')

##-----------------------------------------#
## Difficulties
##-----------------------------------------#

estimateDiff <- data.frame(estimate =  c(bimodalRes$paraEstimates$beta,
                                         bimodalRes$bnpEstimates$beta), 
                             CI_low = c(bimodalRes$paraLow$beta,
                                        bimodalRes$bnpLow$beta), 
                             CI_upp = c(bimodalRes$paraUpper$beta,
                                        bimodalRes$bnpUpper$beta))

estimateDiff$trueVal <- rep(bimodalRes$truValues$beta, 2)
estimateDiff$Model <- rep(c("Parametric", "Semiparametric"), each = length(bimodalRes$truValues$beta))

pBiDiff <- ggplot(estimateDiff, aes(x = trueVal, y = estimate, color = Model)) + 
        geom_point(size = 1.5) + 
          scale_x_continuous(breaks = seq(-3, 3, by = 1)) +
          scale_y_continuous(breaks = seq(-3, 3, by = 1)) +
        geom_errorbar(aes(ymin=CI_low, ymax=CI_upp), width = 0.2) + 
        scale_color_manual(values=c(paraColor, bnpColor),
            guide = guide_legend(override.aes = list(linetype = rep("blank", 2), 
                                                     shape    = c(16, 16)))) +
        labs(x = "True value", y = "Estimate", title = "Difficulty parameters") + 
        geom_abline(color = 'grey10', lty = 3)

pBiDiff

ggsave(filename = "figures/bimodal_difficulties.png", plot = pBiDiff,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')


##-----------------------------------------#
## Figure 6
## Item estimates unimodal simulation
##-----------------------------------------#
unimodalItemPlot <- plot_grid(
    plot_grid(
    pUniDiff  + theme(legend.position = "none"),
    pUniDiscr + theme(legend.position = "none"),
    nrow = 1, align = "h"
   ),
   get_legend(pUniDiscr + theme(legend.position = "bottom")), rel_heights = c(1, .1), nrow=2)

ggsave(filename = "figures/fig6_unimodal_items.png", plot = unimodalItemPlot,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')

##-----------------------------------------#
## Figure 7
## Item estimates bimodal simulation
##-----------------------------------------#

bimodalItemPlot <- plot_grid(
    plot_grid(
    pBiDiff  + theme(legend.position = "none"),
    pBiDiscr + theme(legend.position = "none"),
    nrow = 1, align = "h"
   ),
   get_legend(pBiDiscr + theme(legend.position = "bottom")), rel_heights = c(1, .1), nrow=2)

ggsave(filename = "figures/fig7_bimodal_items.png", plot = bimodalItemPlot,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')
