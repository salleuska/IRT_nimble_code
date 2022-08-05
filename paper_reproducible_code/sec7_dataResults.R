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

title <- paste0("Health data\nEstimate of the distribution of ability")

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
        geom_point(size = 1.5, position = position_dodge(width = 0.6)) + 
          scale_y_continuous(breaks = seq(0, 100, by = 10)) +
          scale_x_continuous(breaks = seq(1, 50, by = 2)) +
            theme(axis.text.x = element_blank()) + 
        geom_errorbar(aes(ymin=CI_low*100, ymax=CI_upp*100),
         width = 0.8, position = position_dodge(width = 0.6)) + 
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

title <- paste0("TIMSS data\nEstimate of the distribution of ability")

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
        geom_point(size = 1.5, position = position_dodge(width = 0.6)) + 
          scale_y_continuous(breaks = seq(0, 100, by = 10)) +
          scale_x_continuous(breaks = seq(1, 50, by = 2)) +
            theme(axis.text.x = element_blank()) + 
        geom_errorbar(aes(ymin=CI_low*100, ymax=CI_upp*100), 
        	width = 0.8, position = position_dodge(width = 0.6)) + 
        scale_color_manual(values=c(paraColor, bnpColor),
            guide = guide_legend(override.aes = list(linetype = rep("blank", 2), 
                                                     shape    = c(16, 16)))) +
        labs(y = "Percentile", x = "Individual", 
                title = "TIMSS data") 

pTimssPerc

ggsave(filename = "figures/timss_percentiles.png", plot = pTimssPerc,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')

##-----------------------------------------#
## Data Timss 3PL
##-----------------------------------------#
timssRes <- readRDS("figures/dataForFigures/timss3PL.rds")
##-----------------------------------------#
## Discriminations
##-----------------------------------------#

estimateDiscr3PL <- data.frame(Parametric     = timssRes$paraEstimates$lambda,
                Semiparametric = timssRes$bnpEstimates$lambda)


pTimssDiscr3PL <- ggplot(estimateDiscr3PL, aes(x = Parametric, y =  Semiparametric)) + 
        geom_point(size = 1.5) + 
        scale_x_continuous(breaks = seq(0.5, 3, by = 0.5)) +
        scale_y_continuous(breaks = seq(0.5, 3, by = 0.5)) +
          labs(x = "Parametric", y = "Semiparametric", title = "Discrimination parameters") + 
        geom_abline(color = 'grey10', lty = 3)

pTimssDiscr3PL

ggsave(filename = "figures/3PL_data_timss_discriminations.png", plot = pTimssDiscr3PL,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')


##-----------------------------------------#
## Difficulties
##-----------------------------------------#

estimateDiff3PL <- data.frame(Parametric     = timssRes$paraEstimates$beta,
               Semiparametric = timssRes$bnpEstimates$beta)


pTimssDiff3PL <- ggplot(estimateDiff3PL, aes(x = Parametric, y =  Semiparametric)) + 
        geom_point(size = 1.5) + 
        scale_x_continuous(breaks = seq(-3, 3, by = 1)) +
        scale_y_continuous(breaks = seq(-3, 3, by = 1)) +
          labs(x = "Parametric", y = "Semiparametric", title = "Difficulty parameters") + 
        geom_abline(color = 'grey10', lty = 3)

pTimssDiff3PL

ggsave(filename = "figures/3PL_data_timss_difficulties.png", plot = pTimssDiff3PL,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')

##-----------------------------------------#
## Guessing
##-----------------------------------------#

estimateGuess3PL <- data.frame(Parametric     = timssRes$paraEstimates$delta,
                           Semiparametric = timssRes$bnpEstimates$delta)


pTimssGuess3PL <- ggplot(estimateGuess3PL, aes(x = Parametric, y =  Semiparametric)) + 
        geom_point(size = 1.5) + 
        scale_x_continuous(breaks = seq(0, 0.5, by = 0.1)) +
        scale_y_continuous(breaks = seq(0, 0.5, by = 0.1)) +
          labs(x = "Parametric", y = "Semiparametric", title = "Guessing parameters") + 
        geom_abline(color = 'grey10', lty = 3)

pTimssGuess3PL

ggsave(filename = "figures/3PL_data_timss_guessing.png", plot = pTimssGuess3PL,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')

##-----------------------------------------#
## Density - mean 
##-----------------------------------------#

dfEtaMeansTimss3PL <- data.frame(mean = c(timssRes$paraEstimates$eta,
                       timssRes$bnpEstimates$eta))

dfEtaMeansTimss3PL$Model <- rep(c('Parametric', 'Semiparametric'), each = dim(dfEtaMeansTimss3PL)[1]/2)

title <- paste0("TIMSS data- 3PL\nDistribution of individual posterior mean abilities")

pEtaMeanTimss3PL <- ggplot(dfEtaMeansTimss3PL, aes(x=mean, color = Model)) +
          geom_histogram(aes(y=..density..), position="identity",
            binwidth = 0.1, fill="white",key_glyph = "path")+
          geom_line(stat = "density", lwd = 1, key_glyph = "path") + 
          xlim(-10, 6) +
            labs(title=title,x="Ability", y = "Density")+
          scale_color_manual(values=c(paraColor, bnpColor),
            guide=guide_legend(override.aes=list(linetype=c(1,1))))


pEtaMeanTimss3PL <- pEtaMeanTimss3PL + theme(legend.title = element_blank())
pEtaMeanTimss3PL

ggsave(filename = "figures/3PL_data_timss_posterior_means.png", plot = pEtaMeanTimss3PL,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')
##-----------------------------------------#
## Density - DP Measure
##-----------------------------------------#

dfEtaDensity3PL <- data.frame(grid            = timssRes$grid, 
                   semiparametric  = apply(timssRes$densityDPMeasure, 2, mean),
                   parametric      = apply(timssRes$densitySamplesPara, 2, mean))


dfEtaDensityPlot3PL <- reshape2::melt(dfEtaDensity3PL, id.vars = "grid")
colnames(dfEtaDensityPlot3PL) <- c("grid", "Model", "value")
levels(dfEtaDensityPlot3PL$Model) <- c("Semiparametric", "Parametric")
dfEtaDensityPlot3PL$Model <- factor(dfEtaDensityPlot3PL$Model, levels = c("Parametric", "Semiparametric"))

dfIntervalPara <- data.frame(grid  = timssRes$grid, 
                             lower = apply(timssRes$densitySamplesPara, 2, quantile, 0.025), 
                             upper = apply(timssRes$densitySamplesPara, 2, quantile, 0.975))
dfIntervalPara$Model <- "Parametric"

dfIntervalBNP <- data.frame(grid  = timssRes$grid, 
                             lower = apply(timssRes$densityDPMeasure, 2, quantile, 0.025), 
                             upper = apply(timssRes$densityDPMeasure, 2, quantile, 0.975))

dfIntervalBNP$Model <- "Semiparametric"

title <- paste0("TIMSS data - 3PL\nEstimate of the distribution of ability")

pEtaDensityTimss3PL <- ggplot(dfEtaDensityPlot3PL, aes(x=grid, y = value, color = Model)) +
            geom_line(lwd = 1, aes(linetype=Model)) +
            geom_line(data = dfIntervalBNP, aes(x = grid, y=lower), lwd = 1,  linetype="dashed") + 
            geom_line(data = dfIntervalBNP, aes(x = grid, y=upper), lwd = 1,  linetype="dashed") +
            geom_line(data = dfIntervalPara, aes(x = grid, y=lower), lwd = 1, linetype="dashed") + 
            geom_line(data = dfIntervalPara, aes(x = grid, y=upper), lwd = 1, linetype="dashed") +
            xlim(-10, 6) +
            scale_linetype_manual(values=c("solid", "solid")) +
            labs(title=title, x="Ability", y = "Density")+
            scale_color_manual(values=c(paraColor, bnpColor),
                guide=guide_legend(override.aes=list(linetype=c("solid", "solid"))))

pEtaDensityTimss3PL <- pEtaDensityTimss3PL + theme(legend.title = element_blank())
pEtaDensityTimss3PL

ggsave(filename = "figures/3PL_data_timss_posterior_density.png", plot = pEtaDensityTimss3PL,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')


##-----------------------------------------#
## Percentiles
##-----------------------------------------#
dfPercentile3PL <- data.frame(ind       = rep(1:50, 2),
                           estimate  = c(apply(timssRes$paraPerc, 2, mean),
                                         apply(timssRes$bnpPerc, 2, mean)),
                           CI_low    = c(apply(timssRes$paraPerc, 2, quantile, 0.025),
                                         apply(timssRes$bnpPerc, 2, quantile, 0.025)),
                           CI_upp    = c(apply(timssRes$paraPerc, 2, quantile, 0.975),
                                         apply(timssRes$bnpPerc, 2, quantile, 0.975)))

dfPercentile3PL$Model <- rep(c("Parametric", "Semiparametric"), each = dim(dfPercentile3PL)[1]/2)

pTimssPerc3PL <- ggplot(dfPercentile3PL, aes(x = ind, y = estimate*100, color = Model)) + 
        geom_point(size = 1.5, position = position_dodge(width = 0.6)) + 
          scale_y_continuous(breaks = seq(0, 100, by = 10)) +
          scale_x_continuous(breaks = seq(1, 50, by = 2)) +
            theme(axis.text.x = element_blank()) + 
        geom_errorbar(aes(ymin=CI_low*100, ymax=CI_upp*100), 
          width = 0.8, position = position_dodge(width = 0.6)) + 
        scale_color_manual(values=c(paraColor, bnpColor),
            guide = guide_legend(override.aes = list(linetype = rep("blank", 2), 
                                                     shape    = c(16, 16)))) +
        labs(y = "Percentile", x = "Individual", 
                title = "TIMSS data") 

pTimssPerc3PL

ggsave(filename = "figures/3PL_data_timss_percentiles.png", plot = pTimssPerc3PL,
        width = plot_width, height = plot_height , dpi = 300, units = unit, device='png')

# ##-----------------------------------------#
# ## Figure 11
# ## Item estimates health
# ##-----------------------------------------#

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


# ##-----------------------------------------#
# ## Figure 13
# ##-----------------------------------------#

timssItemPlot <- plot_grid(
    pTimssDiff3PL + theme(legend.position = "none"),
    pTimssDiscr3PL  + theme(legend.position = "none"),
    pTimssGuess3PL  + theme(legend.position = "none"),
    nrow = 1, align = "h")

timssItemPlot

ggsave(filename = "figures/fig13_timss3PL_item.png", plot = timssItemPlot,
        width = plot_width, height = plot_height*0.9 ,
        scale = 1.3, 
        dpi = 300, units = unit, device='png')

# ##-----------------------------------------#
# ## Figure 14
# ## Distribution of posterior means & density - timssdata
# ##-----------------------------------------#

densityTimss <- plot_grid(
    plot_grid(
    pEtaMeanTimss3PL  + theme(legend.position = "none"),
    pEtaDensityTimss3PL + theme(legend.position = "none"),
    nrow = 1, align = "h"
   ),
   get_legend(pEtaDensityTimss + theme(legend.position = "bottom")), rel_heights = c(1, .1), nrow=2)
densityTimss


ggsave(filename = "figures/fig14_data_timss3PL_density_estimate.png", plot = densityTimss,
        width = plot_width, height = plot_height ,
        scale = 1.2, dpi = 300, units = unit, device='png')


# ##-----------------------------------------#
# ## Figure 15
# ## Percentile estimation - health & timss
# ##-----------------------------------------#

dataPercPlot <- plot_grid(
    plot_grid(
    pHealthPerc  + theme(legend.position = "none") + theme(axis.text = element_text(size= 15),
                axis.title=element_text(size=15),
                strip.text = element_text(size=15),
                plot.title = element_text(size=17)),
    pTimssPerc3PL + theme(legend.position = "none")+ 
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



# ##-----------------------------------------#
# ## Supplementary materials - timss 2PL
# ##-----------------------------------------#

timssItemPlot <- plot_grid(
    pTimssDiff + theme(legend.position = "none"),
    pTimssDiscr  + theme(legend.position = "none"),
    nrow = 1, align = "h")

timssItemPlot

ggsave(filename = "figures/SM_timss2PL_item.png", plot = timssItemPlot,
        width = plot_width, height = plot_height , scale = 1.2,
        dpi = 300, units = unit, device='png')

# ##-----------------------------------------#
# ## Figure SM
# ## Distribution of posterior means & density - timssdata
# ##-----------------------------------------#

densityTimss <- plot_grid(
    plot_grid(
    pEtaMeanTimss  + theme(legend.position = "none"),
    pEtaDensityTimss + theme(legend.position = "none"),
    nrow = 1, align = "h"
   ),
   get_legend(pEtaDensityTimss + theme(legend.position = "bottom")), rel_heights = c(1, .1), nrow=2)
densityTimss


ggsave(filename = "figures/SM_timss2PL_density_estimate.png", plot = densityTimss,
        width = plot_width, height = plot_height ,
        scale = 1.2, dpi = 300, units = unit, device='png')