library(ggplot2)

outdir <- "graphs/"
########################################
## samples form  parametric models
########################################
dir <- "simulated_samples/parametric/"
modelType <- "parametric"
listFiles <- list.files(dir)

priorSamples  <- matrix(0, ncol = length(listFiles), nrow = 10*100*100)

for(i in 1:length(listFiles)) {
	ff <- readRDS(paste0(dir, listFiles[i]))
	priorSamples[,i] <- as.numeric(ff)
}

colnames(priorSamples) <- gsub("priorSamples_|\\.rds", "", listFiles)

########################################
## samples from bnp models
########################################
dir <- "simulated_samples/bnp/"

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

apply(priorSamples, 2, function(x) density(x)$bw)

p_hist <- ggplot(dd, aes(x= value, col = Var2)) + 
               geom_histogram(aes(y = ..density..),  
               binwidth=0.03, fill = "white") +
               facet_wrap(.~ Var2)  +
               theme_bw()  + theme(legend.position = 'none') +
               xlab("Probability") + ylab("") +
               theme(strip.background = element_rect(color="black", 
               	fill="white", size=1, linetype="solid"))

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

p <- ggplot(dd_qq, aes(sample= value, col = variable)) + 
               stat_qq(distribution = qbeta, dparams = c(0.5, 0.5)) +
               stat_qq_line(distribution = qbeta, dparams = c(0.5, 0.5), col = 1) +
               xlab("theoretical - Beta(0.5, 0.5)") + ylab("Empirical") +
               facet_wrap(.~ variable) +   
               theme_bw()  + theme(legend.position = 'none')

p  

ggsave(filename = paste0(outdir, modelType, "_qqplot.png"), plot = p,
        width = 8, height = 6, dpi = 400, units = "in", device='png')