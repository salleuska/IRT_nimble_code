##-----------------------------------------#
## Computational strategies and estimation performance with Bayesian semiparametric Item Response Theory model
## Sally Paganin
## last update: August 2022
## R version 4.2.0 (2022-04-22) -- "Vigorous Calisthenics"
## nimble version 0.12.2
##-----------------------------------------#
library(ggplot2)
library(cowplot)
library(bayestestR)
source("R_functions/ggplot_settings.R")
source("R_functions/multimodalDensity.R")
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
        geom_errorbar(aes(ymin=CI_low, ymax=CI_upp), 
        	width = 0.05) + 
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
pEtaMeanUni

ggsave(filename = "figures/unimodal_posterior_means.png", plot = pEtaMeanUni,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')
##-----------------------------------------#
## Density - DP Measure
##-----------------------------------------#

dfEtaDensity <- data.frame(grid  = unimodalRes$grid, 
		              semiparametric  = apply(unimodalRes$densityDPMeasure, 2, mean),
    		    	      parametric      = apply(unimodalRes$densitySamplesPara, 2, mean))

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

title <- paste0("Unimodal simulation\nEstimate of the distribution of ability")

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
        geom_point(size = 1.5, position = position_dodge(width = 0.6)) + 
          scale_y_continuous(breaks = seq(0, 100, by = 10)) +
          scale_x_continuous(breaks = seq(1, 50, by = 2)) +
            theme(axis.text.x = element_blank(), legend.title=element_blank()) + 
        geom_errorbar(aes(ymin=CI_low*100, ymax=CI_upp*100), 
        	width = 0.8, position = position_dodge(width = 0.6)) + 
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

dfEtaDensity <- data.frame(grid = bimodalRes$grid, 
        semiparametric  = apply(bimodalRes$densityDPMeasure, 2, mean),
        parametric 	= apply(bimodalRes$densitySamplesPara, 2, mean))

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

title <- paste0("Bimodal simulation\nEstimate of the distribution of ability")

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
        geom_point(size = 1.5, position = position_dodge(width = 0.6)) + 
          scale_y_continuous(breaks = seq(0, 100, by = 10)) +
          scale_x_continuous(breaks = seq(1, 50, by = 2)) +
            theme(axis.text.x = element_blank(), legend.title = element_blank()) + 
        geom_errorbar(aes(ymin=CI_low*100, ymax=CI_upp*100), 
        	width = 0.8, position = position_dodge(width = 0.6)) + 
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
## Multimodal
##-----------------------------------------#
dataName <- "simulation_multimodal2"
multimodalRes <- readRDS("figures/dataForFigures/multimodal2.rds")

##-----------------------------------------#
## Discriminations
##-----------------------------------------#
estimateDiscr <- data.frame(estimate = c(multimodalRes$paraEstimates$lambda,
                                           multimodalRes$bnpEstimates$lambda), 
                             CI_low  = c(multimodalRes$paraLow$lambda,
                                        multimodalRes$bnpLow$lambda), 
                             CI_upp  = c(multimodalRes$paraUpper$lambda,
                                        multimodalRes$bnpUpper$lambda))

estimateDiscr$trueVal <- rep(multimodalRes$truValues$lambda, 2)
estimateDiscr$Model <- rep(c("Parametric", "Semiparametric"), each = length(multimodalRes$truValues$lambda))

pMultiDiscr <- ggplot(estimateDiscr, aes(x = trueVal, y = estimate, color = Model)) + 
        geom_point(size = 1.5) + 
          scale_x_continuous(breaks = seq(0.2, 2, by = 0.2)) +
          scale_y_continuous(breaks = seq(0.2, 2, by = 0.2)) + 
        geom_errorbar(aes(ymin=CI_low, ymax=CI_upp), 
                width = 0.05) + 
        scale_color_manual(values=c(paraColor, bnpColor),
            guide = guide_legend(override.aes = list(linetype = rep("blank", 2), 
                                                     shape    = c(16, 16)))) +
        labs(x = "True value", y = "Estimate", title = "Discrimination parameters") + 
        geom_abline(color = 'grey10', lty = 3)

pMultiDiscr
ggsave(filename = "figures/multimodal_discriminations.png", plot = pMultiDiscr,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')

##-----------------------------------------#
## Difficulties
##-----------------------------------------#

estimateDiff <- data.frame(estimate =  c(multimodalRes$paraEstimates$beta,
                                         multimodalRes$bnpEstimates$beta), 
                             CI_low = c(multimodalRes$paraLow$beta,
                                        multimodalRes$bnpLow$beta), 
                             CI_upp = c(multimodalRes$paraUpper$beta,
                                        multimodalRes$bnpUpper$beta))

estimateDiff$trueVal <- rep(multimodalRes$truValues$beta, 2)
estimateDiff$Model <- rep(c("Parametric", "Semiparametric"), each = length(multimodalRes$truValues$beta))

pMultiDiff <- ggplot(estimateDiff, aes(x = trueVal, y = estimate, color = Model)) + 
        geom_point(size = 1.5) + 
          scale_x_continuous(breaks = seq(-3, 3, by = 1)) +
          scale_y_continuous(breaks = seq(-3, 3, by = 1)) +
        geom_errorbar(aes(ymin=CI_low, ymax=CI_upp), width = 0.2) + 
        scale_color_manual(values=c(paraColor, bnpColor),
            guide = guide_legend(override.aes = list(linetype = rep("blank", 2), 
                                                     shape    = c(16, 16)))) +
        labs(x = "True value", y = "Estimate", title = "Difficulty parameters") + 
        geom_abline(color = 'grey10', lty = 3)

pMultiDiff

ggsave(filename = "figures/multimodal_difficulties.png", plot = pMultiDiff,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')

##-----------------------------------------#
## Density - mean 
##-----------------------------------------#

dfEtaMeansMulti <- data.frame(mean = c(multimodalRes$paraEstimates$eta,
                                                                         multimodalRes$bnpEstimates$eta))

dfEtaMeansMulti$Model <- rep(c('Parametric', 'Semiparametric'), each = dim(dfEtaMeansMulti)[1]/2)

title <- paste0("Multimodal simulation\nDistribution of individual posterior mean abilities")

pEtaMeanMulti <- ggplot(dfEtaMeansMulti, aes(x=mean, color = Model)) +
                            geom_histogram(aes(y=..density..), position="identity",
                              binwidth = 0.1, fill="white",key_glyph = "path")+
                            geom_line(stat = "density", lwd = 1, key_glyph = "path") + 
                            xlim(-8, 8) +
                              labs(title=title,x="Ability", y = "Density")+
                            scale_color_manual(values=c(paraColor, bnpColor, "black"),
                              guide=guide_legend(override.aes=list(linetype=c(1,1,3))))+
                            stat_function(
                                    fun = function(x) dMultiModal(x,
                                        weights = c(2,4,4), 
                                        means= c(-2, 0, 3)), 
                                    aes(col = 'True density'),lwd = 1, lty = 3) 

pEtaMeanMulti <- pEtaMeanMulti + theme(legend.title = element_blank())
pEtaMeanMulti

ggsave(filename = "figures/multimodal_posterior_means.png", plot = pEtaMeanMulti,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')
##-----------------------------------------#
## Density - DP Measure
##-----------------------------------------#

dfEtaDensity <- data.frame(grid  = multimodalRes$grid, 
                              semiparametric  = apply(multimodalRes$densityDPMeasure, 2, mean),
                              parametric      = apply(multimodalRes$densitySamplesPara, 2, mean))


dfEtaDensity$trueDensity <- dMultiModal(dfEtaDensity$grid,
                                        weights = c(2,4,4), 
                                        means= c(-2, 0, 3))

dfEtaDensityPlot <- reshape2::melt(dfEtaDensity, id.vars = "grid")
colnames(dfEtaDensityPlot) <- c("grid", "Model", "value")
levels(dfEtaDensityPlot$Model) <- c("Semiparametric", "Parametric", "True density")
dfEtaDensityPlot$Model <- factor(dfEtaDensityPlot$Model, levels = c("Parametric", "Semiparametric", "True density"))

dfIntervalPara <- data.frame(grid  = multimodalRes$grid, 
                             lower = apply(multimodalRes$densitySamplesPara, 2, quantile, 0.025), 
                             upper = apply(multimodalRes$densitySamplesPara, 2, quantile, 0.975))
dfIntervalPara$Model <- "Parametric"

dfIntervalBNP <- data.frame(grid  = multimodalRes$grid, 
                             lower = apply(multimodalRes$densityDPMeasure, 2, quantile, 0.025), 
                             upper = apply(multimodalRes$densityDPMeasure, 2, quantile, 0.975))

dfIntervalBNP$Model <- "Semiparametric"


title <- paste0("Multimodal simulation\nEstimate of the distribution of ability")

pEtaDensityMulti <- ggplot(dfEtaDensityPlot, aes(x=grid, y = value, color = Model)) +
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

pEtaDensityMulti <- pEtaDensityMulti + theme(legend.title = element_blank())
pEtaDensityMulti
ggsave(filename = "figures/multimodal_posterior_density.png", plot = pEtaDensityMulti,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')

##-----------------------------------------#
## Percentiles
##-----------------------------------------#
dfPercentile <- data.frame(ind       = rep(1:50, 2),
                           estimate  = c(apply(multimodalRes$paraPerc, 2, mean),
                                         apply(multimodalRes$bnpPerc, 2, mean)),
                           CI_low    = c(apply(multimodalRes$paraPerc, 2, quantile, 0.025),
                                         apply(multimodalRes$bnpPerc, 2, quantile, 0.025)),
                           CI_upp    = c(apply(multimodalRes$paraPerc, 2, quantile, 0.975),
                                         apply(multimodalRes$bnpPerc, 2, quantile, 0.975)))


dfPercentile$trueVal <- rep(multimodalRes$truePerc, 2)
dfPercentile$Model <- rep(c("Parametric", "Semiparametric"), each = dim(dfPercentile)[1]/2)

pMultiPerc <- ggplot(dfPercentile, aes(x = ind, y = estimate*100, color = Model)) + 
        geom_point(size = 1.5, position = position_dodge(width = 0.6)) + 
          scale_y_continuous(breaks = seq(0, 100, by = 10)) +
          scale_x_continuous(breaks = seq(1, 50, by = 2)) +
            theme(axis.text.x = element_blank(), legend.title=element_blank()) + 
        geom_errorbar(aes(ymin=CI_low*100, ymax=CI_upp*100), 
                width = 0.8, position = position_dodge(width = 0.6)) + 
        geom_point(aes(x = ind, y = trueVal*100, fill = "True value"), color = "black", size = 1.5) + 
        scale_color_manual(values=c(paraColor, bnpColor),
            guide = guide_legend(override.aes = list(linetype = rep("blank"), 
                                                     shape    = c(16, 16),
                                                     color    = c(paraColor, bnpColor)))) + 
        labs(y = "Percentile", x = "Individual", 
                title = "Multimodal simulation")  
pMultiPerc

ggsave(filename = "figures/multimodal_percentiles.png", plot = pMultiPerc,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')

##-----------------------------------------#
## Figure 6
## Item estimates unimodal simulation
##-----------------------------------------#

#   allItemPlot <- plot_grid(
#     plot_grid(
#     pUniDiff  + theme(legend.position = "none"),
#     pUniDiscr + theme(legend.position = "none"),
#     pBiDiff  + theme(legend.position = "none"),
#     pBiDiscr + theme(legend.position = "none"),
#     pMultiDiff  + theme(legend.position = "none"),
#     pMultiDiscr + theme(legend.position = "none"),
#     nrow = 3, align = "h",
#     labels= c("sim ", "quo", "qua")),
#    get_legend(pUniDiscr + theme(legend.position = "bottom", legend.title = element_blank())), rel_heights = c(1, .1), nrow=2)

# ggsave(filename = "figures/fig6_unimodal_items.png", plot = unimodalItemPlot,
#         width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')




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


multiItemPlot <- plot_grid(
    plot_grid(
    pMultiDiff  + theme(legend.position = "none"),
    pMultiDiscr + theme(legend.position = "none"),
    nrow = 1, align = "h"
   ),
   get_legend(pMultiDiscr + theme(legend.position = "bottom", legend.title = element_blank())), 
              rel_heights = c(1, .1), nrow=2)


title1 <- ggdraw() + draw_label("Bimodal Simulation", 
      fontface='bold')
title2 <- ggdraw() + draw_label("Multimodal Simulation", 
      fontface='bold')

allPlots <-  plot_grid(title1, bimodalItemPlot, title2,multiItemPlot,
                  rel_heights = c(.1, 1, .1, 1), ncol =1)

ggsave(filename = "figures/fig7_bimodal_multimodal_items.png", plot = allPlots,
        width = plot_width, height = plot_height*1.8 , dpi = 300, units = unit, device='png')
##-----------------------------------------#
## Figure 8
## Item estimates multimodal simulation
##-----------------------------------------#


# ggsave(filename = "figures/fig8_multimodal_items.png", plot = multiItemPlot,
#         width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')


##-----------------------------------------#
## Figure 9
## Distribution of posterior means - unimodal & bimodal & multi
##-----------------------------------------#

densityPosteriorMeans <- plot_grid(
    plot_grid(
    pEtaMeanUni  + theme(legend.position = "none"),
    pEtaMeanBi + theme(legend.position = "none"),
    pEtaMeanMulti + theme(legend.position = "none"),
    nrow = 1, align = "h"
   ),
   get_legend(pEtaDensityBi + theme(legend.position = "bottom")), rel_heights = c(1, .1), nrow=2)

ggsave(filename = "figures/fig8_density_means.png", plot = densityPosteriorMeans,
        width = plot_width*1.2, height = plot_height, 
        scale = 1.3,
        dpi = 300, units = unit, device='png')
##-----------------------------------------#
## Figure 10
## density- unimodal & bimodal & multi
##-----------------------------------------#

densityPosteriorPred <- plot_grid(
    plot_grid(
    pEtaDensityUni  + theme(legend.position = "none"),
    pEtaDensityBi + theme(legend.position = "none"),
    pEtaDensityMulti + theme(legend.position = "none"),
    nrow = 1, align = "h"
   ),
   get_legend(pEtaDensityBi + theme(legend.position = "bottom")), rel_heights = c(1, .1), nrow=2)

ggsave(filename = "figures/fig9_density_estimate.png", plot = densityPosteriorPred,
        width = plot_width*1.2, height = plot_height ,
        scale = 1.3,
         dpi = 300, units = unit, device='png')


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
    pMultiPerc + theme(legend.position = "none")+ 
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
        width = plot_width*1.2, height = plot_height,
        scale = 1.4,
         dpi = 300, units = unit, device='png')







