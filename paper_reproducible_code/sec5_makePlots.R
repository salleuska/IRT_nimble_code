library(ggplot2)
library(cowplot)
## Global ggplot theme
theme_set(theme_bw() +
          theme(axis.text = element_text(size= 16),
                axis.title=element_text(size=18),
                strip.text = element_text(size=18),
                strip.background =element_rect(fill="white"),
                plot.title = element_text(size=16), 
                legend.text=element_text(size=16)))

outdir <- "figures/"
########################################
## samples form  parametric models
########################################
dir <- "output/prior_samples/parametric/"
modelType <- "parametric"
listFiles <- list.files(dir)

priorSamples  <- matrix(0, ncol = length(listFiles), nrow = 10*100*100)

for(i in 1:length(listFiles)) {
	ff <- readRDS(paste0(dir, listFiles[i]))
	priorSamples[,i] <- as.numeric(ff)
}

colnames(priorSamples) <- gsub("priorSamples_|\\.rds", "", listFiles)

dd <- reshape2::melt(priorSamples)
str(dd)

levels(dd$Var2) <- gsub("\\_", " ", levels(dd$Var2))
levels(dd$Var2) <- gsub("([[:lower:]])([[:upper:]][[:lower:]])", "\\1 \\2", levels(dd$Var2))
levels(dd$Var2) <- gsub("Abilities", "abilities", levels(dd$Var2))
levels(dd$Var2) <- gsub("Item", "item", levels(dd$Var2))


p_hist <- ggplot(dd, aes(x= value)) + 
               geom_histogram(aes(y = ..density..),  
               bins = 50, fill = "white", col = "gray50") +
               facet_wrap(.~ Var2)  +
               theme_bw()  + theme(legend.position = 'none') +
               xlab("") + ylab("") +
               ggtitle("Parametric models") + 
               theme(strip.background = element_rect(color="black", 
                fill="white", linetype="solid"))

p_hist


p_hist <- p_hist  + ylim(0, 8) + stat_function(fun = dbeta, colour="black", lty = 2,  n = 200, xlim = c(0.0001, 0.9999),
                      args = list(shape1 = 0.5, shape2 = 0.5)) 
p_hist

ggsave(filename = paste0(outdir, modelType, "_dist.png"), plot = p_hist,
        width = 8, height = 6, dpi = 400, units = "in", device='png')

################
## qqplots 
################
qq_df <- apply(priorSamples, 2,  sort)
qq_df <- as.data.frame(qq_df)

qq_df$y <- qbeta(ppoints(length(qq_df[,1])), 0.5, 0.5) 

vv <- round(seq.int(1, dim(qq_df)[1], length = 5000))
qq_df <- qq_df[vv, ]
str(qq_df)

dd_qq <- reshape2::melt(qq_df, id = "y")
str(dd_qq)
# levels(dd_qq$variable) <- paste0("alpha = ", levels(dd_qq$variable))

p <- ggplot(dd_qq, aes(sample= value)) +
               scale_y_continuous(breaks = seq(0.0, 1.0, length = 5)) + 
               scale_x_continuous(breaks = seq(0.0, 1.0, length = 5)) + 
               stat_qq(distribution = qbeta, dparams = c(0.5, 0.5), col = "gray") +
               stat_qq_line(distribution = qbeta, dparams = c(0.5, 0.5), col = 1) +
               xlab("Theoretical quantiles - Beta(0.5, 0.5)") + ylab("Empirical quantiles") +
               facet_wrap(.~ variable) +   
               theme_bw()  + theme(legend.position = 'none')+
               theme(strip.background = element_rect(color="black", 
                fill="white", linetype="solid"))

p  

ggsave(filename = paste0(outdir, modelType, "_qqplot.png"), plot = p,
        width = 8, height = 6, dpi = 400, units = "in", device='png')

########################################
## samples from bnp models
########################################
dir <- "output/prior_samples/bnp/"
modelType <- "bnp"

listFiles <- list.files(dir)

priorSamples  <- matrix(0, ncol = length(listFiles), nrow = 10*100*100)

for(i in 1:length(listFiles)) {
	ff <- readRDS(paste0(dir, listFiles[i]))
	priorSamples[,i] <- as.numeric(ff)
}

colnames(priorSamples) <- gsub("priorSamples_|\\.rds", "", listFiles)
str(priorSamples)

##########################################
dd <- reshape2::melt(priorSamples)
str(dd)

if(length(levels(dd$Var2)) == 4) levels(dd$Var2) <- c("IRT_constrainedItem", "IRT_unconstrained" ,"SI_constrainedItem","SI_unconstrained")

levels(dd$Var2) <- gsub("\\_", " ", levels(dd$Var2))
levels(dd$Var2) <- gsub("([[:lower:]])([[:upper:]][[:lower:]])", "\\1 \\2", levels(dd$Var2))
levels(dd$Var2) <- gsub("Abilities", "abilities", levels(dd$Var2))
levels(dd$Var2) <- gsub("Item", "item", levels(dd$Var2))

## binwidth=0.03,
p_hist_bnp <- ggplot(dd, aes(x= value)) + 
               geom_histogram(aes(y = ..density..),  
               bins = 50, col = "gray50", fill = "white") +
               facet_wrap(.~ Var2)  +
               theme_bw()  + theme(legend.position = 'none') +
               xlab("") + ylab("") +
               ggtitle("semiparametric models") + 
               theme(strip.background = element_rect(color="black", 
               	fill="white", linetype="solid"))

p_hist_bnp


p_hist_bnp <- p_hist_bnp  + ylim(0, 8) + stat_function(fun = dbeta, colour="black", lty = 2,  n = 200, xlim = c(0.0001, 0.9999),
                      args = list(shape1 = 0.5, shape2 = 0.5)) 
p_hist_bnp


ggsave(filename = paste0(outdir, modelType, "_dist.png"), plot = p_hist_bnp,
        width = 8/3*2, height = 6, dpi = 400, units = "in", device='png')

###############################
### Subset for paper plot 

dd_subset <- droplevels(subset(dd, dd$Var2 %in% levels(dd$Var2)[grep("a 2 b 4", levels(dd$Var2))]))
levels(dd_subset$Var2) <- gsub(" a 2 b 4", "", levels(dd_subset$Var2))

p_hist_subset <- ggplot(dd_subset, aes(x= value)) + 
               geom_histogram(aes(y = ..density..),  
               bins = 50, col = "gray50", fill = "white") +
               facet_wrap(.~ Var2)  +
               theme_bw()  + theme(legend.position = 'none') +
               xlab("") + ylab("") +
               ggtitle("Semiparametric models") + 
               theme(strip.background = element_rect(color="black", 
                fill="white", linetype="solid"))

p_hist_subset

p_hist_subset <- p_hist_subset  + ylim(0, 8) + stat_function(fun = dbeta, colour="black", lty = 2,  n = 200, xlim = c(0.0001, 0.9999),
                      args = list(shape1 = 0.5, shape2 = 0.5)) 
p_hist_subset

plot_para_bnp <- plot_grid(p_hist, p_hist_subset, rel_widths = c(3/5, 2/5))
plot_para_bnp

ggsave(filename = paste0(outdir,"fig2_prior_simulations.png"), plot = plot_para_bnp,
        width = 12, height = 6, dpi = 400, units = "in", device='png')

################
## qqplots 
################
qq_df <- apply(priorSamples, 2,  sort)
qq_df <- as.data.frame(qq_df)

qq_df$y <- qbeta(ppoints(length(qq_df[,1])), 0.5, 0.5) 

vv <- round(seq.int(1, dim(qq_df)[1], length = 5000))
qq_df <- qq_df[vv, ]
str(qq_df)

dd_qq <- reshape2::melt(qq_df, id = "y")
str(dd_qq)
# levels(dd_qq$variable) <- paste0("alpha = ", levels(dd_qq$variable))

p <- ggplot(dd_qq, aes(sample= value)) + ylim(-0.1,1.1)+
               scale_y_continuous(breaks = seq(0.0, 1.0, length = 5)) + 
               scale_x_continuous(breaks = seq(0.0, 1.0, length = 5)) + 
               stat_qq(distribution = qbeta, dparams = c(0.5, 0.5), col = "gray") +
               stat_qq_line(distribution = qbeta, dparams = c(0.5, 0.5), col = 1) +
               xlab("Theoretical quantiles - Beta(0.5, 0.5)") + ylab("Empirical quantiles") +
               facet_wrap(.~ variable) +   
               theme_bw()  + theme(legend.position = 'none')+
               theme(strip.background = element_rect(color="black", 
                fill="white", linetype="solid"))

p  

ggsave(filename = paste0(outdir, modelType, "_qqplot.png"), plot = p,
        width = 8, height = 6, dpi = 400, units = "in", device='png')


########################################
## samples from bnp models - fixed alpha
########################################
dir <- "output/prior_samples/bnpfixedAlpha/"
modelType <- "bnpFixedAlpha"

listFiles <- list.files(dir)
# listFiles <- listFiles[c(3, 6, 9, 12)]

priorSamples  <- matrix(0, ncol = length(listFiles), nrow = 10*100*100)

for(i in 1:length(listFiles)) {
  ff <- readRDS(paste0(dir, listFiles[i]))
  priorSamples[,i] <- as.numeric(ff)
}

colnames(priorSamples) <- gsub("priorSamples_|\\.rds", "", listFiles)
str(priorSamples)

##########################################
dd <- reshape2::melt(priorSamples)
str(dd)

p_hist <- ggplot(dd, aes(x= value)) + 
               geom_histogram(aes(y = ..density..),  
               bins = 50, col = "gray50", fill = "white") +
               facet_wrap(.~ Var2, ncol = 6)  +
               theme_bw()  + theme(legend.position = 'none') +
               xlab("") + ylab("") +
               ggtitle("Distribution of probabilities samples") + 
               theme(strip.background = element_rect(color="black", 
                fill="white", linetype="solid"))

p_hist


p_hist <- p_hist  + ylim(0, 8) + stat_function(fun = dbeta, colour="black", lty = 2,  n = 200, xlim = c(0.0001, 0.9999),
                      args = list(shape1 = 0.5, shape2 = 0.5)) 
p_hist

ggsave(filename = paste0(outdir, modelType, "_dist.png"), plot = p_hist,
        width = 14, height = 10, dpi = 400, units = "in", device='png')

################
## qqplots 
################
qq_df <- apply(priorSamples, 2,  sort)
qq_df <- as.data.frame(qq_df)

qq_df$y <- qbeta(ppoints(length(qq_df[,1])), 0.5, 0.5) 

vv <- round(seq.int(1, dim(qq_df)[1], length = 5000))
qq_df <- qq_df[vv, ]
str(qq_df)

dd_qq <- reshape2::melt(qq_df, id = "y")
str(dd_qq)
# levels(dd_qq$variable) <- paste0("alpha = ", levels(dd_qq$variable))

p <- ggplot(dd_qq, aes(sample= value)) + 
               scale_y_continuous(breaks = seq(0.0, 1.0, length = 5)) + 
               scale_x_continuous(breaks = seq(0.0, 1.0, length = 5)) + 
               stat_qq(distribution = qbeta, dparams = c(0.5, 0.5), col = "gray") +
               stat_qq_line(distribution = qbeta, dparams = c(0.5, 0.5), col = 1) +
               xlab("Theoretical quantiles - Beta(0.5, 0.5)") + ylab("Empirical quantiles") +
               facet_wrap(.~ variable, ncol = 6) +   
               theme_bw()  + theme(legend.position = 'none')+
               theme(strip.background = element_rect(color="black", 
                fill="white", linetype="solid"))

p  

ggsave(filename = paste0(outdir, modelType, "_qqplot.png"), plot = p,
        width = 14, height = 10, dpi = 400, units = "in", device='png')
