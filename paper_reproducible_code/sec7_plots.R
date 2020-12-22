##-----------------------------------------#
## Bayesian semiparametric Item Response Theory models using NIMBLE 
## Sally Paganin
## November 2020
##-----------------------------------------#
library(ggplot2)
library(cowplot)
library(bayestestR)
source("R_functions/ggplot_settings.R")
##-----------------------------------------#
## Set dimensions
plot_width  <- 21*1.1 ## a4 paper
plot_height <- plot_width/5*2
unit   <- "cm" 

##-----------------------------------------#
## Unimodal
##-----------------------------------------#
dataName <- "simulation_unimodal"
unimodalRes <- readRDS("figures/dataForFigures/unimodal.rds")

##-----------------------------------------#
## Discriminations
##-----------------------------------------#

estimateDiscr <- data.frame(estimate = c(unimodalRes$paraEstimates$lambda,
                                           unimodalRes$bnpEstimates$lambda), 
                             CI_low  = c(unimodalRes$paraLow$lambda,
                                        unimodalRes$bnpLow$lambda), 
                             CI_upp  = c(unimodalRes$paraUpper$lambda,
                                        unimodalRes$bnpUpper$lambda))

estimateDiscr$trueVal <- rep(unimodalRes$truValues$lambda, 2)
estimateDiscr$Model <- rep(c("Parametric", "Semiparametric"), each = length(unimodalRes$truValues$lambda))

pUniDiscr <- ggplot(estimateDiscr, aes(x = trueVal, y = estimate, color = Model)) + 
        geom_point(size = 1.5) + 
          scale_x_continuous(breaks = seq(0.2, 2, by = 0.2)) +
          scale_y_continuous(breaks = seq(0.2, 2, by = 0.2)) + 
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

dfEtaMeansUni <- data.frame(mean = c(unimodalRes$paraEstimates$eta,
					   				 unimodalRes$bnpEstimates$eta))

dfEtaMeansUni$Model <- rep(c('Parametric', 'Semiparametric'), each = dim(dfEtaMeansUni)[1]/2)

title <- paste0("Unimodal simulation\nDistribution of individual posterior mean abilities")

pEtaMeanUni <- ggplot(dfEtaMeansUni, aes(x=mean, color = Model)) +
			    geom_histogram(aes(y=..density..), position="identity",
			      binwidth = 0.1, fill="white",key_glyph = "path")+
			    geom_line(stat = "density", lwd = 1, key_glyph = "path") + 
			    xlim(-8, 8) +
			      labs(title=title,x="Ability", y = "Density")+
			    scale_color_manual(values=c(paraColor, bnpColor, "black"),
			      guide=guide_legend(override.aes=list(linetype=c(1,1,3))))+
			    stat_function(
			            fun = function(x) dnorm(x, 0, sd = 1.25), 
			            aes(col = 'True density'),lwd = 1, lty = 3) 

pEtaMeanUni <- pEtaMeanUni + theme(legend.title = element_blank())

ggsave(filename = "figures/unimodal_posterior_means.png", plot = pEtaMeanUni,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')
##-----------------------------------------#
## Density - DP Measure
##-----------------------------------------#

dfEtaDensity <- data.frame(grid            = unimodalRes$grid, 
		        		   semiparametric  = apply(unimodalRes$densityDPMeasure, 2, mean),
    		    		   parametric 	   = apply(unimodalRes$densitySamplesPara, 2, mean))

dfEtaDensity$trueDensity <- dnorm(dfEtaDensity$grid, mean = 0, sd = 1.25)

dfEtaDensityPlot <- reshape2::melt(dfEtaDensity, id.vars = "grid")
colnames(dfEtaDensityPlot) <- c("grid", "Model", "value")
levels(dfEtaDensityPlot$Model) <- c("Semiparametric", "Parametric", "True density")
dfEtaDensityPlot$Model <- factor(dfEtaDensityPlot$Model, levels = c("Parametric", "Semiparametric", "True density"))

dfIntervalPara <- data.frame(grid  = unimodalRes$grid, 
                             lower = apply(unimodalRes$densitySamplesPara, 2, quantile, 0.025), 
                             upper = apply(unimodalRes$densitySamplesPara, 2, quantile, 0.975))
dfIntervalPara$Model <- "Parametric"

dfIntervalBNP <- data.frame(grid  = unimodalRes$grid, 
                             lower = apply(unimodalRes$densityDPMeasure, 2, quantile, 0.025), 
                             upper = apply(unimodalRes$densityDPMeasure, 2, quantile, 0.975))

dfIntervalBNP$Model <- "Semiparametric"

title <- paste0("Unimodal simulation\nEstimate of distribution of abilities")

pEtaDensityUni <- ggplot(dfEtaDensityPlot, aes(x=grid, y = value, color = Model)) +
		        geom_line(lwd = 1, aes(linetype=Model)) +
		        geom_line(data = dfIntervalBNP, aes(x = grid, y=lower), lwd = 1,  linetype="dashed") + 
		        geom_line(data = dfIntervalBNP, aes(x = grid, y=upper), lwd = 1,  linetype="dashed") +
		        geom_line(data = dfIntervalPara, aes(x = grid, y=lower), lwd = 1, linetype="dashed") + 
		        geom_line(data = dfIntervalPara, aes(x = grid, y=upper), lwd = 1, linetype="dashed") +
		        xlim(-8, 8) +
		        scale_linetype_manual(values=c("solid", "solid", "dotted")) +
		        labs(title=title, x="Ability", y = "Density")+
		        scale_color_manual(values=c(paraColor, bnpColor, "black"),
		            guide=guide_legend(override.aes=list(linetype=c("solid", "solid", "dotted"))))

pEtaDensityUni <- pEtaDensityUni + theme(legend.title = element_blank())
pEtaDensityUni
ggsave(filename = "figures/unimodal_posterior_density.png", plot = pEtaDensityUni,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')

##-----------------------------------------#
## Percentiles
##-----------------------------------------#
dfPercentile <- data.frame(ind       = rep(1:50, 2),
                           estimate  = c(apply(unimodalRes$paraPerc, 2, mean),
                                         apply(unimodalRes$bnpPerc, 2, mean)),
                           CI_low    = c(apply(unimodalRes$paraPerc, 2, quantile, 0.025),
                                         apply(unimodalRes$bnpPerc, 2, quantile, 0.025)),
                           CI_upp    = c(apply(unimodalRes$paraPerc, 2, quantile, 0.975),
                                         apply(unimodalRes$bnpPerc, 2, quantile, 0.975)))



dfPercentile$trueVal <- rep(unimodalRes$truePerc, 2)
dfPercentile$Model <- rep(c("Parametric", "Semiparametric"), each = dim(dfPercentile)[1]/2)

pUniPerc <- ggplot(dfPercentile, aes(x = ind, y = estimate*100, color = Model)) + 
        geom_point(size = 1.5) + 
          scale_y_continuous(breaks = seq(0, 100, by = 10)) +
          scale_x_continuous(breaks = seq(1, 50, by = 2)) +
            theme(axis.text.x = element_blank(), legend.title=element_blank()) + 
        geom_errorbar(aes(ymin=CI_low*100, ymax=CI_upp*100), width = 0.8) + 
        geom_point(aes(x = ind, y = trueVal*100, fill = "True value"), color = "black", size = 1.5) + 
        scale_color_manual(values=c(paraColor, bnpColor),
            guide = guide_legend(override.aes = list(linetype = rep("blank"), 
                                                     shape    = c(16, 16),
                                                     color    = c(paraColor, bnpColor)))) + 
        labs(y = "Percentile", x = "Individual", 
                title = "Unimodal simulation") 
pUniPerc

ggsave(filename = "figures/unimodal_percentiles.png", plot = pUniPerc,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')

##-----------------------------------------#
##-----------------------------------------#
## Bimodal
##-----------------------------------------#
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
          scale_x_continuous(breaks = seq(0.2, 1.8, by = 0.2)) +
          scale_y_continuous(breaks = seq(0.2, 1.8, by = 0.2)) + 
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
## Density - mean 
##-----------------------------------------#

dfEtaMeansBi <- data.frame(mean = c(bimodalRes$paraEstimates$eta,
					   				bimodalRes$bnpEstimates$eta))

dfEtaMeansBi$Model <- rep(c('Parametric', 'Semiparametric'), each = dim(dfEtaMeansBi)[1]/2)

title <- paste0("Bimodal simulation\nDistribution of individual posterior mean abilities")

pEtaMeanBi <- ggplot(dfEtaMeansBi, aes(x=mean, color = Model)) +
			    geom_histogram(aes(y=..density..), position="identity",
			      binwidth = 0.1, fill="white",key_glyph = "path")+
			    geom_line(stat = "density", lwd = 1, key_glyph = "path") + 
			    xlim(-8, 8) +
			      labs(title=title,x="Ability", y = "Density")+
			    scale_color_manual(values=c(paraColor, bnpColor, "black"),
			      guide=guide_legend(override.aes=list(linetype=c(1,1,3))))+
			    stat_function(
			            fun = function(x) 0.5*dnorm(x, -2, sd = 1.25) + 0.5*dnorm(x, 2, sd = 1.25), 
			            aes(col = 'True density'),lwd = 1, lty = 3) 

pEtaMeanBi <- pEtaMeanBi + theme(legend.title = element_blank())

ggsave(filename = "figures/bimodal_posterior_means.png", plot = pEtaMeanBi,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')
##-----------------------------------------#
## Density - DP Measure
##-----------------------------------------#

dfEtaDensity <- data.frame(grid            = bimodalRes$grid, 
		        		   semiparametric  = apply(bimodalRes$densityDPMeasure, 2, mean),
    		    		   parametric 	   = apply(bimodalRes$densitySamplesPara, 2, mean))

dfEtaDensity$trueDensity <- 0.5*dnorm(dfEtaDensity$grid, -2, sd = 1.25) + 0.5*dnorm(dfEtaDensity$grid, 2, sd = 1.25)

dfEtaDensityPlot <- reshape2::melt(dfEtaDensity, id.vars = "grid")
colnames(dfEtaDensityPlot) <- c("grid", "Model", "value")
levels(dfEtaDensityPlot$Model) <- c("Semiparametric", "Parametric", "True density")
dfEtaDensityPlot$Model <- factor(dfEtaDensityPlot$Model, levels = c("Parametric", "Semiparametric", "True density"))

dfIntervalPara <- data.frame(grid  = bimodalRes$grid, 
                             lower = apply(bimodalRes$densitySamplesPara, 2, quantile, 0.025), 
                             upper = apply(bimodalRes$densitySamplesPara, 2, quantile, 0.975))

dfIntervalPara$Model <- "Parametric"

dfIntervalBNP <- data.frame(grid  = bimodalRes$grid, 
                             lower = apply(bimodalRes$densityDPMeasure, 2, quantile, 0.025), 
                             upper = apply(bimodalRes$densityDPMeasure, 2, quantile, 0.975))

dfIntervalBNP$Model <- "Semiparametric"

title <- paste0("Bimodal simulation\nEstimate of distribution of abilities")

pEtaDensityBi <- ggplot(dfEtaDensityPlot, aes(x=grid, y = value, color = Model)) +
		        geom_line(lwd = 1, aes(linetype=Model)) +
		        geom_line(data = dfIntervalBNP, aes(x = grid, y=lower), lwd = 1,  linetype="dashed") + 
		        geom_line(data = dfIntervalBNP, aes(x = grid, y=upper), lwd = 1,  linetype="dashed") +
		        geom_line(data = dfIntervalPara, aes(x = grid, y=lower), lwd = 1, linetype="dashed") + 
		        geom_line(data = dfIntervalPara, aes(x = grid, y=upper), lwd = 1, linetype="dashed") +
		        xlim(-8, 8) +
		        scale_linetype_manual(values=c("solid", "solid", "dotted")) +
		        labs(title=title, x="Ability", y = "Density")+
		        scale_color_manual(values=c(paraColor, bnpColor, "black"),
		            guide=guide_legend(override.aes=list(linetype=c("solid", "solid", "dotted"))))

pEtaDensityBi <- pEtaDensityBi + theme(legend.title = element_blank())

ggsave(filename = "figures/bimodal_posterior_density.png", plot = pEtaDensityBi,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')

##-----------------------------------------#
## Percentiles
##-----------------------------------------#
dfPercentile <- data.frame(ind       = rep(1:50, 2),
                           estimate  = c(apply(bimodalRes$paraPerc, 2, mean),
                                         apply(bimodalRes$bnpPerc, 2, mean)),
                           CI_low    = c(apply(bimodalRes$paraPerc, 2, quantile, 0.025),
                                         apply(bimodalRes$bnpPerc, 2, quantile, 0.025)),
                           CI_upp    = c(apply(bimodalRes$paraPerc, 2, quantile, 0.975),
                                         apply(bimodalRes$bnpPerc, 2, quantile, 0.975)))


dfPercentile$trueVal <- rep(bimodalRes$truePerc, 2)
dfPercentile$Model <- rep(c("Parametric", "Semiparametric"), each = dim(dfPercentile)[1]/2)

pBiPerc <- ggplot(dfPercentile, aes(x = ind, y = estimate*100, color = Model)) + 
        geom_point(size = 1.5) + 
          scale_y_continuous(breaks = seq(0, 100, by = 10)) +
          scale_x_continuous(breaks = seq(1, 50, by = 2)) +
            theme(axis.text.x = element_blank(), legend.title = element_blank()) + 
        geom_errorbar(aes(ymin=CI_low*100, ymax=CI_upp*100), width = 0.8) + 
        geom_point(aes(x = ind, y = trueVal*100, fill = "True value"), color = "black", size = 1.5) + 
        scale_color_manual(values=c(paraColor, bnpColor),
            guide = guide_legend(override.aes = list(linetype = rep("blank", 2), 
                                                     shape    = c(16, 16)))) +
        labs(y = "Percentile", x = "Individual", 
                title = "Bimodal simulation") 
pBiPerc

ggsave(filename = "figures/bimodal_percentiles.png", plot = pBiPerc,
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
   get_legend(pUniDiscr + theme(legend.position = "bottom", legend.title = element_blank())), rel_heights = c(1, .1), nrow=2)

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
   get_legend(pBiDiscr + theme(legend.position = "bottom", legend.title = element_blank())), 
              rel_heights = c(1, .1), nrow=2)

ggsave(filename = "figures/fig7_bimodal_items.png", plot = bimodalItemPlot,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')


##-----------------------------------------#
## Figure 8
## Distribution of posterior means - unimodal & bimodal
##-----------------------------------------#

densityPosteriorMeans <- plot_grid(
    plot_grid(
    pEtaMeanUni  + theme(legend.position = "none"),
    pEtaMeanBi + theme(legend.position = "none"),
    nrow = 1, align = "h"
   ),
   get_legend(pEtaDensityBi + theme(legend.position = "bottom")), rel_heights = c(1, .1), nrow=2)

ggsave(filename = "figures/fig8_density_means.png", plot = densityPosteriorMeans,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')

##-----------------------------------------#
## Figure 9
## density- unimodal & bimodal
##-----------------------------------------#

densityPosteriorPred <- plot_grid(
    plot_grid(
    pEtaDensityUni  + theme(legend.position = "none"),
    pEtaDensityBi + theme(legend.position = "none"),
    nrow = 1, align = "h"
   ),
   get_legend(pEtaDensityBi + theme(legend.position = "bottom")), rel_heights = c(1, .1), nrow=2)

ggsave(filename = "figures/fig9_density_estimate.png", plot = densityPosteriorPred,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')


##-----------------------------------------#
## Figure 10
## Percentile estimation - unimodal & bimodal
##-----------------------------------------#

simuPercPlot <- plot_grid(
    plot_grid(
    pUniPerc  + theme(legend.position = "none") + theme(axis.text = element_text(size= 15),
                axis.title=element_text(size=15),
                strip.text = element_text(size=15),
                plot.title = element_text(size=17)),
    pBiPerc + theme(legend.position = "none")+ 
     theme(axis.text = element_text(size= 15),
                axis.title=element_text(size=15),
                strip.text = element_text(size=15),
                plot.title = element_text(size=17)),
    nrow = 1, align = "h"
   ),
   get_legend(pBiPerc + theme(legend.position = "bottom", legend.text=element_text(size=15), legend.title = element_blank())), 
              rel_heights = c(1, .1), nrow=2)

simuPercPlot

ggsave(filename = "figures/fig10_simulation_percentiles.png", plot = simuPercPlot,
        width = plot_width*1.6, height = plot_height*1.6 , dpi = 300, units = unit, device='png')

#########################################################################################################################
#########################################################################################################################
##-----------------------------------------#
## Data Health
##-----------------------------------------#
healthRes <- readRDS("figures/dataForFigures/health.rds")
##-----------------------------------------#
## Discriminations
##-----------------------------------------#

estimateDiscr <- data.frame(estimate =  c(healthRes$paraEstimates$lambda,
                                           healthRes$bnpEstimates$lambda), 
                             CI_low = c(healthRes$paraLow$lambda,
                                        healthRes$bnpLow$lambda), 
                             CI_upp = c(healthRes$paraUpper$lambda,
                                        healthRes$bnpUpper$lambda))

estimateDiscr$items <- itemDataHealth
estimateDiscr$Model <- rep(c("Parametric", "Semiparametric"), each = dim(estimateDiscr)[1]/2)

pHealthDiscr <- ggplot(estimateDiscr, aes(y = reorder(items, estimate), x = estimate, color = Model)) + 
        geom_point(size = 1.5) + 
          scale_x_continuous(breaks = seq(0.2, 2, by = 0.2)) +
        geom_errorbar(aes(xmin=CI_low, xmax=CI_upp), width= 0.5) + 
        scale_color_manual(values=c(paraColor, bnpColor),
            guide = guide_legend(override.aes = list(linetype = rep("blank", 2), 
                                                     shape    = c(16, 16)))) +
        labs(x = "Estimate", y = "", title = "Discrimination parameters") 

pHealthDiscr

ggsave(filename = "figures/data_health_discriminations.png", plot = pHealthDiscr,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')

##-----------------------------------------#
## Difficulties
##-----------------------------------------#

estimateDiff <- data.frame(estimate =  c(healthRes$paraEstimates$beta,
                                         healthRes$bnpEstimates$beta), 
                             CI_low = c(healthRes$paraLow$beta,
                                        healthRes$bnpLow$beta), 
                             CI_upp = c(healthRes$paraUpper$beta,
                                        healthRes$bnpUpper$beta))

estimateDiff$items <- itemDataHealth
estimateDiff$Model <- rep(c("Parametric", "Semiparametric"), each = dim(estimateDiff)[1]/2)

pHealthDiff <- ggplot(estimateDiff, aes(y = reorder(items, estimate), x = estimate, color = Model)) + 
        geom_point(size = 1.5) + 
		scale_x_continuous(breaks = seq(-3, 6, by = 1)) +
        geom_errorbar(aes(xmin=CI_low, xmax=CI_upp), width= 0.5) + 
        scale_color_manual(values=c(paraColor, bnpColor),
            guide = guide_legend(override.aes = list(linetype = rep("blank", 2), 
                                                     shape    = c(16, 16)))) +
        labs(x = "Estimate", y = "", title = "Difficulty parameters") 

pHealthDiff

ggsave(filename = "figures/data_health_difficulties.png", plot = pHealthDiff,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')

##-----------------------------------------#
## Density - mean 
##-----------------------------------------#

dfEtaMeansHealth <- data.frame(mean = c(healthRes$paraEstimates$eta,
					   				    healthRes$bnpEstimates$eta))

dfEtaMeansHealth$Model <- rep(c('Parametric', 'Semiparametric'), each = dim(dfEtaMeansHealth)[1]/2)

title <- paste0("Health data\nDistribution of individual posterior mean abilities")

pEtaMeanHealth <- ggplot(dfEtaMeansHealth, aes(x=mean, color = Model)) +
			    geom_histogram(aes(y=..density..), position="identity",
			      binwidth = 0.1, fill="white",key_glyph = "path")+
			    geom_line(stat = "density", lwd = 1, key_glyph = "path") + 
			    xlim(-10, 20) +
			      labs(title=title,x="Ability", y = "Density")+
			    scale_color_manual(values=c(paraColor, bnpColor),
			      guide=guide_legend(override.aes=list(linetype=c(1,1))))


pEtaMeanHealth <- pEtaMeanHealth + theme(legend.title = element_blank())
pEtaMeanHealth
ggsave(filename = "figures/data_health_posterior_means.png", plot = pEtaMeanHealth,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')
##-----------------------------------------#
## Density - DP Measure
##-----------------------------------------#

dfEtaDensity <- data.frame(grid            = healthRes$grid, 
		        		   semiparametric  = apply(healthRes$densityDPMeasure, 2, mean),
    		    		   parametric 	   = apply(healthRes$densitySamplesPara, 2, mean))


dfEtaDensityPlot <- reshape2::melt(dfEtaDensity, id.vars = "grid")
colnames(dfEtaDensityPlot) <- c("grid", "Model", "value")
levels(dfEtaDensityPlot$Model) <- c("Semiparametric", "Parametric")
dfEtaDensityPlot$Model <- factor(dfEtaDensityPlot$Model, levels = c("Parametric", "Semiparametric"))

dfIntervalPara <- data.frame(grid  = healthRes$grid, 
                             lower = apply(healthRes$densitySamplesPara, 2, quantile, 0.025), 
                             upper = apply(healthRes$densitySamplesPara, 2, quantile, 0.975))
dfIntervalPara$Model <- "Parametric"

dfIntervalBNP <- data.frame(grid  = healthRes$grid, 
                             lower = apply(healthRes$densityDPMeasure, 2, quantile, 0.025), 
                             upper = apply(healthRes$densityDPMeasure, 2, quantile, 0.975))

dfIntervalBNP$Model <- "Semiparametric"

title <- paste0("Health data\nEstimate of distribution of abilities")

pEtaDensityHealth <- ggplot(dfEtaDensityPlot, aes(x=grid, y = value, color = Model)) +
		        geom_line(lwd = 1, aes(linetype=Model)) +
		        geom_line(data = dfIntervalBNP, aes(x = grid, y=lower), lwd = 1,  linetype="dashed") + 
		        geom_line(data = dfIntervalBNP, aes(x = grid, y=upper), lwd = 1,  linetype="dashed") +
		        geom_line(data = dfIntervalPara, aes(x = grid, y=lower), lwd = 1, linetype="dashed") + 
		        geom_line(data = dfIntervalPara, aes(x = grid, y=upper), lwd = 1, linetype="dashed") +
		        xlim(-15, 30) +
		        scale_linetype_manual(values=c("solid", "solid")) +
		        labs(title=title, x="Ability", y = "Density")+
		        scale_color_manual(values=c(paraColor, bnpColor),
		            guide=guide_legend(override.aes=list(linetype=c("solid", "solid"))))

pEtaDensityHealth <- pEtaDensityHealth + theme(legend.title = element_blank())
pEtaDensityHealth

ggsave(filename = "figures/data_health_posterior_density.png", plot = pEtaDensityHealth,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')

##-----------------------------------------#
## Percentiles
##-----------------------------------------#
dfPercentile <- data.frame(ind       = rep(1:50, 2),
                           estimate  = c(apply(healthRes$paraPerc, 2, mean),
                                         apply(healthRes$bnpPerc, 2, mean)),
                           CI_low    = c(apply(healthRes$paraPerc, 2, quantile, 0.025),
                                         apply(healthRes$bnpPerc, 2, quantile, 0.025)),
                           CI_upp    = c(apply(healthRes$paraPerc, 2, quantile, 0.975),
                                         apply(healthRes$bnpPerc, 2, quantile, 0.975)))


dfPercentile$Model <- rep(c("Parametric", "Semiparametric"), each = dim(dfPercentile)[1]/2)

pHealthPerc <- ggplot(dfPercentile, aes(x = ind, y = estimate*100, color = Model)) + 
        geom_point(size = 1.5) + 
          scale_y_continuous(breaks = seq(0, 100, by = 10)) +
          scale_x_continuous(breaks = seq(1, 50, by = 2)) +
            theme(axis.text.x = element_blank()) + 
        geom_errorbar(aes(ymin=CI_low*100, ymax=CI_upp*100), width = 0.8) + 
        scale_color_manual(values=c(paraColor, bnpColor),
            guide = guide_legend(override.aes = list(linetype = rep("blank", 2), 
                                                     shape    = c(16, 16)))) +
        labs(y = "Percentile", x = "Individual", 
                title = "Health data") 
pHealthPerc

ggsave(filename = "figures/health_percentiles.png", plot = pHealthPerc,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')

##-----------------------------------------#
## Data Timss
##-----------------------------------------#
timssRes <- readRDS("figures/dataForFigures/timss.rds")
##-----------------------------------------#
## Discriminations
##-----------------------------------------#

estimateDiscr <- data.frame(Parametric     = timssRes$paraEstimates$lambda,
						    Semiparametric = timssRes$bnpEstimates$lambda)


pTimssDiscr <- ggplot(estimateDiscr, aes(x = Parametric, y =  Semiparametric)) + 
				geom_point(size = 1.5) + 
				scale_x_continuous(breaks = seq(0.5, 3, by = 0.5)) +
				scale_y_continuous(breaks = seq(0.5, 3, by = 0.5)) +
    			labs(x = "Parametric", y = "Semiparametric", title = "Discrimination parameters") + 
				geom_abline(color = 'grey10', lty = 3)

pTimssDiscr

ggsave(filename = "figures/data_timss_discriminations.png", plot = pTimssDiscr,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')


##-----------------------------------------#
## Difficulties
##-----------------------------------------#

estimateDiff <- data.frame(Parametric     = timssRes$paraEstimates$beta,
						   Semiparametric = timssRes$bnpEstimates$beta)


pTimssDiff <- ggplot(estimateDiff, aes(x = Parametric, y =  Semiparametric)) + 
				geom_point(size = 1.5) + 
				scale_x_continuous(breaks = seq(-3, 3, by = 1)) +
				scale_y_continuous(breaks = seq(-3, 3, by = 1)) +
    			labs(x = "Parametric", y = "Semiparametric", title = "Difficulty parameters") + 
				geom_abline(color = 'grey10', lty = 3)

pTimssDiff

ggsave(filename = "figures/data_timss_difficulties.png", plot = pTimssDiff,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')

##-----------------------------------------#
## Density - mean 
##-----------------------------------------#

dfEtaMeansTimss <- data.frame(mean = c(timssRes$paraEstimates$eta,
					   				   timssRes$bnpEstimates$eta))

dfEtaMeansTimss$Model <- rep(c('Parametric', 'Semiparametric'), each = dim(dfEtaMeansTimss)[1]/2)

title <- paste0("TIMSS data\nDistribution of individual posterior mean abilities")

pEtaMeanTimss <- ggplot(dfEtaMeansTimss, aes(x=mean, color = Model)) +
			    geom_histogram(aes(y=..density..), position="identity",
			      binwidth = 0.1, fill="white",key_glyph = "path")+
			    geom_line(stat = "density", lwd = 1, key_glyph = "path") + 
			    xlim(-6, 6) +
			      labs(title=title,x="Ability", y = "Density")+
			    scale_color_manual(values=c(paraColor, bnpColor),
			      guide=guide_legend(override.aes=list(linetype=c(1,1))))


pEtaMeanTimss <- pEtaMeanTimss + theme(legend.title = element_blank())

ggsave(filename = "figures/data_timss_posterior_means.png", plot = pEtaMeanTimss,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')
##-----------------------------------------#
## Density - DP Measure
##-----------------------------------------#

dfEtaDensity <- data.frame(grid            = timssRes$grid, 
		        		   semiparametric  = apply(timssRes$densityDPMeasure, 2, mean),
    		    		   parametric 	   = apply(timssRes$densitySamplesPara, 2, mean))


dfEtaDensityPlot <- reshape2::melt(dfEtaDensity, id.vars = "grid")
colnames(dfEtaDensityPlot) <- c("grid", "Model", "value")
levels(dfEtaDensityPlot$Model) <- c("Semiparametric", "Parametric")
dfEtaDensityPlot$Model <- factor(dfEtaDensityPlot$Model, levels = c("Parametric", "Semiparametric"))

dfIntervalPara <- data.frame(grid  = timssRes$grid, 
                             lower = apply(timssRes$densitySamplesPara, 2, quantile, 0.025), 
                             upper = apply(timssRes$densitySamplesPara, 2, quantile, 0.975))
dfIntervalPara$Model <- "Parametric"

dfIntervalBNP <- data.frame(grid  = timssRes$grid, 
                             lower = apply(timssRes$densityDPMeasure, 2, quantile, 0.025), 
                             upper = apply(timssRes$densityDPMeasure, 2, quantile, 0.975))

dfIntervalBNP$Model <- "Semiparametric"

title <- paste0("TIMSS data\nEstimate of distribution of abilities")

pEtaDensityTimss <- ggplot(dfEtaDensityPlot, aes(x=grid, y = value, color = Model)) +
		        geom_line(lwd = 1, aes(linetype=Model)) +
		        geom_line(data = dfIntervalBNP, aes(x = grid, y=lower), lwd = 1,  linetype="dashed") + 
		        geom_line(data = dfIntervalBNP, aes(x = grid, y=upper), lwd = 1,  linetype="dashed") +
		        geom_line(data = dfIntervalPara, aes(x = grid, y=lower), lwd = 1, linetype="dashed") + 
		        geom_line(data = dfIntervalPara, aes(x = grid, y=upper), lwd = 1, linetype="dashed") +
		        xlim(-6, 6) +
		        scale_linetype_manual(values=c("solid", "solid")) +
		        labs(title=title, x="Ability", y = "Density")+
		        scale_color_manual(values=c(paraColor, bnpColor),
		            guide=guide_legend(override.aes=list(linetype=c("solid", "solid"))))

pEtaDensityTimss <- pEtaDensityTimss + theme(legend.title = element_blank())

ggsave(filename = "figures/data_timss_posterior_density.png", plot = pEtaDensityTimss,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')


##-----------------------------------------#
## Percentiles
##-----------------------------------------#
dfPercentile <- data.frame(ind       = rep(1:50, 2),
                           estimate  = c(apply(timssRes$paraPerc, 2, mean),
                                         apply(timssRes$bnpPerc, 2, mean)),
                           CI_low    = c(apply(timssRes$paraPerc, 2, quantile, 0.025),
                                         apply(timssRes$bnpPerc, 2, quantile, 0.025)),
                           CI_upp    = c(apply(timssRes$paraPerc, 2, quantile, 0.975),
                                         apply(timssRes$bnpPerc, 2, quantile, 0.975)))

dfPercentile$Model <- rep(c("Parametric", "Semiparametric"), each = dim(dfPercentile)[1]/2)

pTimssPerc <- ggplot(dfPercentile, aes(x = ind, y = estimate*100, color = Model)) + 
        geom_point(size = 1.5) + 
          scale_y_continuous(breaks = seq(0, 100, by = 10)) +
          scale_x_continuous(breaks = seq(1, 50, by = 2)) +
            theme(axis.text.x = element_blank()) + 
        geom_errorbar(aes(ymin=CI_low*100, ymax=CI_upp*100), width = 0.8) + 
        scale_color_manual(values=c(paraColor, bnpColor),
            guide = guide_legend(override.aes = list(linetype = rep("blank", 2), 
                                                     shape    = c(16, 16)))) +
        labs(y = "Percentile", x = "Individual", 
                title = "TIMSS data") 

pTimssPerc

ggsave(filename = "figures/timss_percentiles.png", plot = pTimssPerc,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')

##-----------------------------------------#
## Figure 11
## Item estimates health
##-----------------------------------------#

healthItemPlot <- plot_grid(
    plot_grid(
    pHealthDiff  + theme(legend.position = "none"),
    pHealthDiscr + theme(legend.position = "none"),
    nrow = 1, align = "h",rel_widths = c(1, 0.9)
   ),
   get_legend(pHealthDiscr + theme(legend.position = "bottom",legend.title = element_blank())), rel_heights = c(1, .1), nrow=2)


healthItemPlot
ggsave(filename = "figures/fig11_health_item.png", plot = healthItemPlot,
        width = plot_width*1.2, height = plot_height*0.9, dpi = 300, units = unit, device='png')

##-----------------------------------------#
## Figure 12
## Distribution of posterior means & density - health data
##-----------------------------------------#

densityHealth <- plot_grid(
    plot_grid(
    pEtaMeanHealth  + theme(legend.position = "none"),
    pEtaDensityHealth + theme(legend.position = "none"),
    nrow = 1, align = "h"
   ),
   get_legend(pEtaDensityHealth + theme(legend.position = "bottom", legend.text=element_text(size=15))), rel_heights = c(1, .1), nrow=2)
densityHealth

ggsave(filename = "figures/fig12_data_health_density_estimate.png", plot = densityHealth,
        width = plot_width*1.5, height = plot_height*1.5 , dpi = 300, units = unit, device='png')


##-----------------------------------------#
## Figure 13
##-----------------------------------------#

timssItemPlot <- plot_grid(
    pTimssDiff + theme(legend.position = "none"),
    pTimssDiscr  + theme(legend.position = "none"),
    nrow = 1, align = "h")

timssItemPlot
ggsave(filename = "figures/fig13_timss_item.png", plot = timssItemPlot,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')

##-----------------------------------------#
## Figure 14
## Distribution of posterior means & density - timssdata
##-----------------------------------------#

densityTimss <- plot_grid(
    plot_grid(
    pEtaMeanTimss  + theme(legend.position = "none"),
    pEtaDensityTimss + theme(legend.position = "none"),
    nrow = 1, align = "h"
   ),
   get_legend(pEtaDensityTimss + theme(legend.position = "bottom")), rel_heights = c(1, .1), nrow=2)
densityTimss


ggsave(filename = "figures/fig14_data_timss_density_estimate.png", plot = densityTimss,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')


##-----------------------------------------#
## Figure 15
## Percentile estimation - health & timss
##-----------------------------------------#

dataPercPlot <- plot_grid(
    plot_grid(
    pHealthPerc  + theme(legend.position = "none") + theme(axis.text = element_text(size= 15),
                axis.title=element_text(size=15),
                strip.text = element_text(size=15),
                plot.title = element_text(size=17)),
    pTimssPerc + theme(legend.position = "none")+ 
     theme(axis.text = element_text(size= 15),
                axis.title=element_text(size=15),
                strip.text = element_text(size=15),
                plot.title = element_text(size=17)),
    nrow = 1, align = "h"
   ),
   get_legend(pHealthPerc + theme(legend.position = "bottom", legend.text=element_text(size=15), legend.title = element_blank())), 
              rel_heights = c(1, .1), nrow=2)

dataPercPlot

ggsave(filename = "figures/fig15_data_percentiles.png", plot = dataPercPlot,
        width = plot_width*1.6, height = plot_height*1.6 , dpi = 300, units = unit, device='png')
